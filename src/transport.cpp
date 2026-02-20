#include "../include/transport.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include "../include/math_constants.h"  // MSVC compatibility for M_PI
#include <stdexcept>

using combaero::thermo::R_GAS;
using combaero::thermo::BOLTZMANN;
using combaero::thermo::AVOGADRO;

// Transport-related helpers factored out of thermo_transport.cpp.

// Linear interpolation helper for collision integrals
double linear_interp(double x,
                     const std::vector<double>& x_values,
                     const std::vector<double>& y_values)
{
    if (x_values.size() != y_values.size()) {
        throw std::invalid_argument("linear_interp: x_values and y_values must have the same size");
    }
    if (x_values.size() < 2) {
        throw std::invalid_argument("linear_interp: at least two points are required");
    }

    auto it = std::lower_bound(x_values.begin(), x_values.end(), x);

    if (it == x_values.begin()) {
        return y_values.front();
    }
    if (it == x_values.end()) {
        return y_values.back();
    }

    std::size_t idx = static_cast<std::size_t>(std::distance(x_values.begin(), it) - 1);

    double x0 = x_values[idx];
    double x1 = x_values[idx + 1];
    double y0 = y_values[idx];
    double y1 = y_values[idx + 1];

    if (x1 == x0) {
        throw std::invalid_argument("linear_interp: x_values must be strictly increasing");
    }

    double t = (x - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}

// Collision integral omega(2,2)
double omega22(double T, double well_depth)
{
    if (T <= 0.0) {
        throw std::invalid_argument("omega22: temperature must be positive");
    }
    if (well_depth <= 0.0) {
        throw std::invalid_argument("omega22: well_depth must be positive");
    }

    static const std::vector<double> T_star_values = {
        0.1, 0.2, 0.5, 1.0, 2.0,
        5.0, 10.0, 20.0, 50.0, 100.0
    };

    static const std::vector<double> omega_values = {
        4.008, 2.995, 2.313, 1.710, 1.276,
        0.922, 0.711, 0.567, 0.432, 0.364
    };

    double T_star = (BOLTZMANN * T) / well_depth;
    double T_star_clamped = std::max(T_star_values.front(),
                                     std::min(T_star, T_star_values.back()));
    return linear_interp(T_star_clamped, T_star_values, omega_values);
}

// -------------------------------------------------------------
// File-local helper: compute pure-species μᵢ and kᵢ in one pass.
// omega22 and the kinetic-theory formula are called once per species,
// not once per public function.
// -------------------------------------------------------------

namespace {

struct PureSpeciesTransport {
    std::vector<double> mu;  // pure-species dynamic viscosity [Pa·s]
    std::vector<double> k;   // pure-species thermal conductivity [W/(m·K)]
};

PureSpeciesTransport pure_species_transport(double T, const std::vector<double>& X)
{
    const std::size_t N = X.size();
    PureSpeciesTransport out;
    out.mu.resize(N);
    out.k.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        const double mw_kg  = molar_masses[i] / 1000.0;           // kg/mol
        const double diam   = transport_props[i].diameter * 1.0e-10; // m
        const double eps    = transport_props[i].well_depth * BOLTZMANN; // J
        const double omega  = omega22(T, eps);
        const double m_mol  = mw_kg / AVOGADRO;                   // kg/molecule

        // Kinetic theory viscosity: μ = (5/16) √(π m k T) / (π σ² Ω)
        out.mu[i] = (5.0 / 16.0) * std::sqrt(M_PI * m_mol * BOLTZMANN * T) /
                    (M_PI * diam * diam * omega);

        // Eucken thermal conductivity: k = μ (cp_mass + f_rot R_specific)
        const double cp_i       = cp_R(i, T) * R_GAS;             // J/(mol·K)
        const double cp_mass_i  = cp_i / mw_kg;                   // J/(kg·K)
        const double R_spec_i   = R_GAS / mw_kg;                  // J/(kg·K)

        double f_rot = 0.0;
        switch (transport_props[i].geometry) {
            case MolecularGeometry::Atom:      f_rot = 0.0; break;
            case MolecularGeometry::Linear:    f_rot = 1.0; break;
            case MolecularGeometry::Nonlinear: f_rot = 1.5; break;
        }

        out.k[i] = out.mu[i] * (cp_mass_i + f_rot * R_spec_i);
    }

    return out;
}

// Wilke mixing rule for viscosity given pre-computed pure-species μᵢ.
double wilke_viscosity(const std::vector<double>& X,
                       const std::vector<double>& pure_mu)
{
    const std::size_t N = X.size();
    double mix_visc = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        if (X[i] <= 0.0) continue;
        double denom = 0.0;
        for (std::size_t j = 0; j < N; ++j) {
            if (X[j] <= 0.0) continue;
            double phi_ij;
            if (i == j) {
                phi_ij = 1.0;
            } else {
                const double M_i   = molar_masses[i];
                const double M_j   = molar_masses[j];
                const double term1 = 1.0 / std::sqrt(8.0) / std::sqrt(1.0 + M_i / M_j);
                const double term2 = std::pow(1.0 + std::sqrt(pure_mu[i] / pure_mu[j]) *
                                     std::pow(M_j / M_i, 0.25), 2);
                phi_ij = term1 * term2;
            }
            denom += X[j] * phi_ij;
        }
        mix_visc += X[i] * pure_mu[i] / denom;
    }
    return mix_visc;
}

// Mixture thermal conductivity (average of upper/lower bounds) given pre-computed kᵢ.
double mixture_conductivity(const std::vector<double>& X,
                            const std::vector<double>& pure_k)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;
        sum1 += X[i] * pure_k[i];
        sum2 += X[i] / pure_k[i];
    }
    return 0.5 * (sum1 + 1.0 / sum2);
}

}  // namespace

// Dynamic viscosity [Pa·s]
double viscosity(double T, double P, const std::vector<double>& X)
{
    (void)P;
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("viscosity: Mole fraction vector size does not match number of species");
    }
    auto ps = pure_species_transport(T, X);
    return wilke_viscosity(X, ps.mu);
}

// Thermal conductivity [W/(m·K)]
double thermal_conductivity(double T, double P, const std::vector<double>& X)
{
    (void)P;
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("thermal_conductivity: Mole fraction vector size does not match number of species");
    }
    auto ps = pure_species_transport(T, X);
    return mixture_conductivity(X, ps.k);
}

// Prandtl number (dimensionless)
double prandtl(double T, double P, const std::vector<double>& X)
{
    (void)P;
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("prandtl: Mole fraction vector size does not match number of species");
    }
    auto ps  = pure_species_transport(T, X);
    double mu_mix = wilke_viscosity(X, ps.mu);
    double k_mix  = mixture_conductivity(X, ps.k);
    if (k_mix <= 0.0) {
        throw std::runtime_error("prandtl: thermal conductivity must be positive");
    }
    const double MW      = mwmix(X) / 1000.0;   // kg/mol
    const double cp_mass = cp(T, X) / MW;        // J/(kg·K)
    return (cp_mass * mu_mix) / k_mix;
}

// Kinematic viscosity [m²/s]
double kinematic_viscosity(double T, double P, const std::vector<double>& X)
{
    return viscosity(T, P, X) / density(T, P, X);
}

// Thermal diffusivity [m²/s] — α = k / (ρ cp)
double thermal_diffusivity(double T, double P, const std::vector<double>& X)
{
    const double k_mix   = thermal_conductivity(T, P, X);
    const double rho     = density(T, P, X);
    const double MW      = mwmix(X) / 1000.0;
    const double cp_mass = cp(T, X) / MW;
    return k_mix / (rho * cp_mass);
}

// Reynolds number — Re = ρ V L / μ
double reynolds(double T, double P, const std::vector<double>& X, double V, double L)
{
    const double rho_val = density(T, P, X);
    const double mu      = viscosity(T, P, X);
    if (mu <= 0.0) {
        throw std::runtime_error("reynolds: viscosity must be positive");
    }
    return (rho_val * V * L) / mu;
}

// Peclet number (thermal) — Pe = V L / α
double peclet(double T, double P, const std::vector<double>& X, double V, double L)
{
    const double alpha = thermal_diffusivity(T, P, X);
    if (alpha <= 0.0) {
        throw std::runtime_error("peclet: thermal diffusivity must be positive");
    }
    return (V * L) / alpha;
}

// -------------------------------------------------------------
// Transport State Bundle — single-pass implementation
// -------------------------------------------------------------

TransportState transport_state(double T, double P, const std::vector<double>& X) {
    if (T <= 0) {
        throw std::invalid_argument("transport_state: temperature must be positive");
    }
    if (P <= 0) {
        throw std::invalid_argument("transport_state: pressure must be positive");
    }
    if (X.empty()) {
        throw std::invalid_argument("transport_state: mole fractions vector cannot be empty");
    }
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("transport_state: mole fraction vector size mismatch");
    }

    // Single species-loop pass — omega22 called N times total
    auto ps = pure_species_transport(T, X);

    const double mu  = wilke_viscosity(X, ps.mu);
    const double k   = mixture_conductivity(X, ps.k);
    const double rho = density(T, P, X);
    const double MW  = mwmix(X) / 1000.0;   // kg/mol
    const double cp_m = cp(T, X) / MW;      // J/(kg·K)
    const double cv_m = cv(T, X) / MW;      // J/(kg·K)

    TransportState state;
    state.T     = T;
    state.P     = P;
    state.rho   = rho;
    state.mu    = mu;
    state.k     = k;
    state.nu    = mu / rho;
    state.alpha = k / (rho * cp_m);
    state.Pr    = (cp_m * mu) / k;
    state.cp    = cp_m;
    state.cv    = cv_m;
    state.gamma = cp_m / cv_m;
    state.a     = speed_of_sound(T, X);
    return state;
}
