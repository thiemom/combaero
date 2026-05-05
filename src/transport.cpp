#include "../include/transport.h"
#include "../include/thermo.h"
#include "thermo.h"
#include "composition.h"
#include "thermo_transport_data.h"
#include "../include/math_constants.h"  // MSVC compatibility for M_PI
#include <stdexcept>

namespace combaero {

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

// Collision integral omega(2,2) -- Monchick-Mason table (37 T* x 8 delta* points).
// Table values from Cantera MMCollisionInt.cpp (authoritative transcription).
// T* axis (37 points):
static const double MM_T_STAR[37] = {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0
};
// delta* axis (8 columns): 0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5
static const double MM_DELTA[8] = {0.0, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0, 2.5};
// omega22 table [37 rows x 8 cols], row-major
static const double MM_OMEGA22[37 * 8] = {
    4.1005, 4.266,  4.833,  5.742,  6.729,  8.624,  10.34,  11.89,
    3.2626, 3.305,  3.516,  3.914,  4.433,  5.57,   6.637,  7.618,
    2.8399, 2.836,  2.936,  3.168,  3.511,  4.329,  5.126,  5.874,
    2.531,  2.522,  2.586,  2.749,  3.004,  3.64,   4.282,  4.895,
    2.2837, 2.277,  2.329,  2.46,   2.665,  3.187,  3.727,  4.249,
    2.0838, 2.081,  2.13,   2.243,  2.417,  2.862,  3.329,  3.786,
    1.922,  1.924,  1.97,   2.072,  2.225,  2.614,  3.028,  3.435,
    1.7902, 1.795,  1.84,   1.934,  2.07,   2.417,  2.788,  3.156,
    1.6823, 1.689,  1.733,  1.82,   1.944,  2.258,  2.596,  2.933,
    1.5929, 1.601,  1.644,  1.725,  1.838,  2.124,  2.435,  2.746,
    1.4551, 1.465,  1.504,  1.574,  1.67,   1.913,  2.181,  2.451,
    1.3551, 1.365,  1.4,    1.461,  1.544,  1.754,  1.989,  2.228,
    1.28,   1.289,  1.321,  1.374,  1.447,  1.63,   1.838,  2.053,
    1.2219, 1.231,  1.259,  1.306,  1.37,   1.532,  1.718,  1.912,
    1.1757, 1.184,  1.209,  1.251,  1.307,  1.451,  1.618,  1.795,
    1.0933, 1.1,    1.119,  1.15,   1.193,  1.304,  1.435,  1.578,
    1.0388, 1.044,  1.059,  1.083,  1.117,  1.204,  1.31,   1.428,
    0.99963, 1.004, 1.016,  1.035,  1.062,  1.133,  1.22,   1.319,
    0.96988, 0.9732, 0.983, 0.9991, 1.021,  1.079,  1.153,  1.236,
    0.92676, 0.9291, 0.936, 0.9473, 0.9628, 1.005,  1.058,  1.121,
    0.89616, 0.8979, 0.903, 0.9114, 0.923,  0.9545, 0.9955, 1.044,
    0.87272, 0.8741, 0.878, 0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
    0.85379, 0.8549, 0.858, 0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
    0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
    0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
    0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
    0.78363, 0.784,  0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
    0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
    0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
    0.74364, 0.7438, 0.7445, 0.7455, 0.747,  0.7512, 0.7569, 0.7642,
    0.71982, 0.72,   0.7204, 0.7211, 0.7221, 0.725,  0.7289, 0.7339,
    0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
    0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
    0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
    0.65099, 0.651,  0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
    0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
    0.5887,  0.5889, 0.5894, 0.59,   0.5903, 0.5901, 0.5895, 0.5885
};

// Collision integral Omega*(2,2) via bilinear interpolation in Monchick-Mason table.
// T_star = T / (epsilon/k_B); delta_star = 0.0 for non-polar species.
double omega22(double T_star, double delta_star)
{
    // Clamp to table bounds
    const double ts = std::max(MM_T_STAR[0], std::min(T_star, MM_T_STAR[36]));
    const double ds = std::max(0.0, std::min(delta_star, MM_DELTA[7]));

    // Find bracketing T* row indices
    int iT = 35;
    for (int k = 0; k < 36; ++k) {
        if (ts <= MM_T_STAR[k + 1]) { iT = k; break; }
    }
    const double fT = (ts - MM_T_STAR[iT]) / (MM_T_STAR[iT + 1] - MM_T_STAR[iT]);

    // Find bracketing delta* column indices
    int iD = 6;
    for (int k = 0; k < 7; ++k) {
        if (ds <= MM_DELTA[k + 1]) { iD = k; break; }
    }
    const double fD = (ds - MM_DELTA[iD]) / (MM_DELTA[iD + 1] - MM_DELTA[iD]);

    // Bilinear interpolation
    const double q00 = MM_OMEGA22[iT       * 8 + iD];
    const double q01 = MM_OMEGA22[iT       * 8 + iD + 1];
    const double q10 = MM_OMEGA22[(iT + 1) * 8 + iD];
    const double q11 = MM_OMEGA22[(iT + 1) * 8 + iD + 1];
    return (1.0 - fT) * ((1.0 - fD) * q00 + fD * q01)
               + fT   * ((1.0 - fD) * q10 + fD * q11);
}

// -------------------------------------------------------------
// File-local helpers
// -------------------------------------------------------------

namespace {

// Reduced dipole moment (Stockmayer parameter) delta*.
// well_depth_K in K, diameter_A in Angstrom, dipole_D in Debye.
double compute_delta_star(double well_depth_K, double diameter_A, double dipole_D)
{
    if (dipole_D == 0.0) return 0.0;
    constexpr double K_E     = 8.98755e9;    // N*m^2/C^2  (1/4*pi*eps0)
    constexpr double D_TO_CM = 3.33564e-30;  // C*m per Debye
    const double mu_SI  = dipole_D * D_TO_CM;
    const double eps_SI = well_depth_K * thermo::BOLTZMANN;
    const double sig_m  = diameter_A * 1.0e-10;
    return K_E * mu_SI * mu_SI / (2.0 * eps_SI * sig_m * sig_m * sig_m);
}

struct PureSpeciesTransport {
    std::vector<double> mu;  // pure-species dynamic viscosity [Pa·s]
    std::vector<double> k;   // pure-species thermal conductivity [W/(m·K)]
};

PureSpeciesTransport pure_species_transport(double T_in, const std::vector<double>& X)
{
    double T = std::max(T_in, 10.0);

    const std::size_t N = X.size();
    PureSpeciesTransport out;
    out.mu.resize(N);
    out.k.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        const double mw_kg  = molar_masses[i] / 1000.0;              // kg/mol
        const double diam   = transport_props[i].diameter * 1.0e-10; // m
        const double T_star = T / transport_props[i].well_depth;
        const double d_star = compute_delta_star(transport_props[i].well_depth,
                                                 transport_props[i].diameter,
                                                 transport_props[i].dipole_moment);
        const double omega  = omega22(T_star, d_star);
        const double m_mol  = mw_kg / thermo::AVOGADRO;              // kg/molecule

        // Kinetic theory viscosity: mu = (5/16) sqrt(pi*m*k_B*T) / (pi*sigma^2*Omega)
        out.mu[i] = (5.0 / 16.0) * std::sqrt(M_PI * m_mol * thermo::BOLTZMANN * T) /
                    (M_PI * diam * diam * omega);

        // Mason-Monchick modified Eucken thermal conductivity
        // Formula from Kee, Coltrin & Glarborg (2003) as implemented in Cantera.
        // All cv_*_R are dimensionless (in units of R per mole).

        // cv_rot in units of R: nrot/2 (Atom=0, Linear=1, Nonlinear=1.5)
        double cv_rot_R = 0.0;
        switch (transport_props[i].geometry) {
            case MolecularGeometry::Atom:      cv_rot_R = 0.0; break;
            case MolecularGeometry::Linear:    cv_rot_R = 1.0; break;
            case MolecularGeometry::Nonlinear: cv_rot_R = 1.5; break;
        }

        // cv_vib from NASA: cp/R - 5/2 - cv_rot/R (clamp to 0 at low T)
        const double cp_R_val = cp_R(i, T);
        const double cv_vib_R = std::max(0.0, cp_R_val - 2.5 - cv_rot_R);

        // rho*D/mu ~ 6/5 (kinetic theory approximation)
        const double f_int_val = 1.2;
        const double A_factor  = 2.5 - f_int_val;  // ~ 1.3

        // Parker temperature correction to Z_rot (Kee 2003 eq. 12.112)
        const double z_rot = transport_props[i].z_rot;
        double z_rot_eff = z_rot;
        if (z_rot > 0.0) {
            auto parker_fz = [](double ts) -> double {
                return 1.0 + std::pow(M_PI, 1.5) / std::sqrt(ts) * (0.5 + 1.0 / ts)
                       + (M_PI * M_PI / 4.0 + 2.0) / ts;
            };
            const double tstar_298 = 298.0 / transport_props[i].well_depth;
            z_rot_eff = z_rot * parker_fz(tstar_298) / parker_fz(T_star);
        }

        const double B_factor = z_rot_eff + (2.0 / M_PI) * (5.0 / 3.0 * cv_rot_R + f_int_val);
        const double c1       = (2.0 / M_PI) * A_factor / B_factor;
        const double f_rot    = f_int_val * (1.0 + c1);
        const double f_trans  = 2.5 * (1.0 - c1 * cv_rot_R / 1.5);

        out.k[i] = (out.mu[i] / mw_kg) * thermo::R_GAS
                   * (f_trans * 1.5 + f_rot * cv_rot_R + f_int_val * cv_vib_R);
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

} // namespace combaero
