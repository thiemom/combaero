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

    double k_B = BOLTZMANN;

    double T_star = (k_B * T) / well_depth;

    double T_star_clamped = std::max(T_star_values.front(),
                                     std::min(T_star, T_star_values.back()));

    return linear_interp(T_star_clamped, T_star_values, omega_values);
}

// Dynamic viscosity [Pa·s] using kinetic theory with per-species transport data
double viscosity(double T, double P, const std::vector<double>& X)
{
    (void)P; // Not used in the current simple model

    if (X.size() != species_names.size()) {
        throw std::invalid_argument("viscosity: Mole fraction vector size does not match number of species");
    }

    // Calculate pure species viscosities using kinetic theory
    std::vector<double> pure_visc(X.size());
    for (std::size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;

        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;

        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;

        // Calculate collision integral
        double omega = omega22(T, eps);

        // Calculate viscosity using kinetic theory formula
        // μ = (5/16) * √(π*m*k*T) / (π*σ²*Ω(2,2))
        // where m is mass of a molecule = MW/NA
        double mass_molecule = mw_kg / AVOGADRO;
        pure_visc[i] = (5.0 / 16.0) * std::sqrt(M_PI * mass_molecule * BOLTZMANN * T) /
                       (M_PI * diam * diam * omega);
    }

    // Wilke's mixing rule for viscosity
    double mix_visc = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;

        double denom = 0.0;
        for (std::size_t j = 0; j < X.size(); ++j) {
            if (X[j] <= 0.0) continue;

            double phi_ij;
            if (i == j) {
                phi_ij = 1.0;
            } else {
                double M_i = molar_masses[i];
                double M_j = molar_masses[j];

                // Wilke's formula
                double term1 = 1.0 / std::sqrt(8.0) / std::sqrt(1.0 + M_i / M_j);
                double term2 = std::pow(1.0 + std::sqrt(pure_visc[i] / pure_visc[j]) *
                               std::pow(M_j / M_i, 0.25), 2);
                phi_ij = term1 * term2;
            }

            denom += X[j] * phi_ij;
        }

        mix_visc += X[i] * pure_visc[i] / denom;
    }

    return mix_visc;
}

// Thermal conductivity [W/(m·K)] using Eucken formula with per-species data
double thermal_conductivity(double T, double P, const std::vector<double>& X)
{
    (void)P; // Not used in the current simple model

    if (X.size() != species_names.size()) {
        throw std::invalid_argument("thermal_conductivity: Mole fraction vector size does not match number of species");
    }

    // Calculate pure species thermal conductivities using Eucken formula
    std::vector<double> pure_cond(X.size());
    for (std::size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;

        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;

        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;

        // Calculate pure species viscosity
        double omega = omega22(T, eps);
        double mass_molecule = mw_kg / AVOGADRO;
        double visc = (5.0 / 16.0) * std::sqrt(M_PI * mass_molecule * BOLTZMANN * T) /
                      (M_PI * diam * diam * omega);

        // Get specific heat capacity for the species in J/(mol·K)
        double cp_species = cp_R(i, T) * R_GAS;

        // Convert cp from J/(mol·K) to J/(kg·K)
        // mw_kg is in kg/mol, so cp_species / mw_kg = J/(kg·K)
        double cp_mass = cp_species / mw_kg;

        // Specific gas constant R_specific = R / M [J/(kg·K)]
        double R_specific = R_GAS / mw_kg;

        // Determine rotational contribution based on geometry
        double f_rot = 0.0;
        switch (transport_props[i].geometry) {
            case MolecularGeometry::Atom:     f_rot = 0.0; break;
            case MolecularGeometry::Linear:   f_rot = 1.0; break;
            case MolecularGeometry::Nonlinear: f_rot = 1.5; break;
        }

        // Modified Eucken formula for thermal conductivity
        // k = mu * (Cv + 9/4 * R) for monatomic
        // k = mu * (Cp + f_rot * R) is a simplified form
        pure_cond[i] = visc * (cp_mass + f_rot * R_specific);
    }

    // Mixture rule (average of upper and lower bounds)
    double sum1 = 0.0;
    double sum2 = 0.0;

    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;

        sum1 += X[i] * pure_cond[i];
        sum2 += X[i] / pure_cond[i];
    }

    return 0.5 * (sum1 + 1.0 / sum2);
}

// Prandtl number (dimensionless)
double prandtl(double T, double P, const std::vector<double>& X)
{
    double mu = viscosity(T, P, X);
    double k = thermal_conductivity(T, P, X);

    double Cp = cp(T, X);
    double MW = mwmix(X) / 1000.0; // kg/mol
    double Cp_mass = Cp / MW;      // J/(kg·K)

    if (k <= 0.0) {
        throw std::runtime_error("prandtl: thermal conductivity must be positive");
    }

    return (Cp_mass * mu) / k;
}

// Kinematic viscosity [m^2/s]
double kinematic_viscosity(double T, double P, const std::vector<double>& X)
{
    double visc_val = viscosity(T, P, X);
    double rho_val = density(T, P, X);

    return visc_val / rho_val;
}

// Thermal diffusivity [m^2/s]
// α = k / (ρ * cp)
double thermal_diffusivity(double T, double P, const std::vector<double>& X)
{
    double k = thermal_conductivity(T, P, X);
    double rho_val = density(T, P, X);

    double Cp = cp(T, X);
    double MW = mwmix(X) / 1000.0;  // kg/mol
    double Cp_mass = Cp / MW;       // J/(kg·K)

    return k / (rho_val * Cp_mass);
}

// Reynolds number (dimensionless)
// Re = ρ * V * L / μ
double reynolds(double T, double P, const std::vector<double>& X, double V, double L)
{
    double rho_val = density(T, P, X);
    double mu = viscosity(T, P, X);

    if (mu <= 0.0) {
        throw std::runtime_error("reynolds: viscosity must be positive");
    }

    return (rho_val * V * L) / mu;
}

// Peclet number (thermal, dimensionless)
// Pe = V * L / α  (ratio of convective to conductive heat transfer)
double peclet(double T, double P, const std::vector<double>& X, double V, double L)
{
    double alpha = thermal_diffusivity(T, P, X);

    if (alpha <= 0.0) {
        throw std::runtime_error("peclet: thermal diffusivity must be positive");
    }

    return (V * L) / alpha;
}

// -------------------------------------------------------------
// Transport State Bundle
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

    TransportState state;

    // Echo back inputs
    state.T = T;
    state.P = P;

    // Compute all transport properties
    state.rho = density(T, P, X);
    state.mu = viscosity(T, P, X);
    state.k = thermal_conductivity(T, P, X);
    state.nu = kinematic_viscosity(T, P, X);
    state.alpha = thermal_diffusivity(T, P, X);
    state.Pr = prandtl(T, P, X);
    state.cp = cp_mass(T, X);

    return state;
}
