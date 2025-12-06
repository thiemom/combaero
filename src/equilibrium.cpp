#include "../include/equilibrium.h"
#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"
#include <algorithm>
#include <cmath>

// -------------------------------------------------------------
// Internal utilities
// -------------------------------------------------------------

namespace {

// Internal WGS configuration (species indices)
struct WgsConfig {
    std::size_t i_CO;
    std::size_t i_H2O;
    std::size_t i_CO2;
    std::size_t i_H2;
};

static inline double clamp(double v, double lo, double hi)
{
    return std::max(lo, std::min(v, hi));
}

static inline void composition_from_xi(const std::vector<double>& n0,
                                       double xi,
                                       const WgsConfig& cfg,
                                       std::vector<double>& n)
{
    n = n0;
    n[cfg.i_CO]  -= xi;
    n[cfg.i_H2O] -= xi;
    n[cfg.i_CO2] += xi;
    n[cfg.i_H2]  += xi;
}

static inline void molefractions(const std::vector<double>& n, std::vector<double>& X)
{
    double nt = 0.0;
    for (double v : n) nt += v;

    X.resize(n.size());
    if (nt <= 0.0) {
        double inv = 1.0 / static_cast<double>(n.size());
        for (std::size_t k = 0; k < n.size(); ++k) X[k] = inv;
    } else {
        for (std::size_t k = 0; k < n.size(); ++k)
            X[k] = n[k] / nt;
    }
}

static inline double Kp_WGS(double T, const WgsConfig& cfg)
{
    double g_CO2 = ::g_over_RT(cfg.i_CO2, T);
    double g_H2  = ::g_over_RT(cfg.i_H2,  T);
    double g_CO  = ::g_over_RT(cfg.i_CO,  T);
    double g_H2O = ::g_over_RT(cfg.i_H2O, T);

    double d_gRT = (g_CO2 + g_H2) - (g_CO + g_H2O);
    return std::exp(-d_gRT);
}

// -------------------------------------------------------------
// Newton with damping
// -------------------------------------------------------------

template<typename Fun>
static double newton_damped(Fun&& F, double x0,
                            double xmin, double xmax,
                            std::size_t maxIt = 40)
{
    double x = x0;

    for (std::size_t it = 0; it < maxIt; ++it) {
        double f0 = F.f(x);

        // Converged?
        if (std::abs(f0) < 1e-12)
            return clamp(x, xmin, xmax);

        // Finite diff derivative
        double eps = 1e-6 * std::max(1.0, std::abs(x));
        double f1  = F.f(x + eps);
        double df  = (f1 - f0) / eps;

        if (std::abs(df) < 1e-14)
            return clamp(x, xmin, xmax);

        double dx = -f0 / df;
        double x_trial = x + dx;

        // Damping if outside bounds
        double lambda = 1.0;
        while ( (x_trial < xmin || x_trial > xmax) && lambda > 1e-6 ) {
            lambda *= 0.5;
            x_trial = x + lambda * dx;
        }

        x = clamp(x_trial, xmin, xmax);
    }

    return clamp(x, xmin, xmax);
}

// -------------------------------------------------------------
// Solve WGS equilibrium at fixed T
// -------------------------------------------------------------

struct WgsIsothermalFun {
    const std::vector<double>& n0;
    WgsConfig cfg;
    double lnKp;
    mutable std::vector<double> n;

    WgsIsothermalFun(const std::vector<double>& n0_, double T,
                     const WgsConfig& cfg_)
        : n0(n0_), cfg(cfg_), n(n0_.size())
    {
        lnKp = std::log(Kp_WGS(T, cfg));
    }

    double f(double xi) const
    {
        composition_from_xi(n0, xi, cfg, n);

        double nt = 0.0;
        for (double v : n) nt += v;

        // Avoid divisions by zero
        const double tiny = 1e-300;
        double yCO2 = std::max(n[cfg.i_CO2] / nt, tiny);
        double yH2  = std::max(n[cfg.i_H2]  / nt, tiny);
        double yCO  = std::max(n[cfg.i_CO]  / nt, tiny);
        double yH2O = std::max(n[cfg.i_H2O] / nt, tiny);

        double Q = (yCO2 * yH2) / (yCO * yH2O);
        return std::log(Q) - lnKp;
    }
};

static double solve_wgs_isothermal(const std::vector<double>& n0,
                                   double T,
                                   const WgsConfig& cfg)
{
    WgsIsothermalFun F{n0, T, cfg};

    // Determine ξ bounds from non-negativity constraints
    // Reaction: CO + H2O -> CO2 + H2 (forward: xi > 0)
    //
    // n[CO]  = n0[CO]  - xi  >= 0  =>  xi <= n0[CO]
    // n[H2O] = n0[H2O] - xi  >= 0  =>  xi <= n0[H2O]
    // n[CO2] = n0[CO2] + xi  >= 0  =>  xi >= -n0[CO2]
    // n[H2]  = n0[H2]  + xi  >= 0  =>  xi >= -n0[H2]

    double xi_min = std::max(-n0[cfg.i_CO2], -n0[cfg.i_H2]);
    double xi_max = std::min(n0[cfg.i_CO], n0[cfg.i_H2O]);

    return newton_damped(F, 0.0, xi_min, xi_max, 40);
}

// -------------------------------------------------------------
// Adiabatic solver
// -------------------------------------------------------------

struct AdiabaticFun {
    const std::vector<double>& n_in;
    double H_target;
    WgsConfig cfg;

    mutable std::vector<double> n;
    mutable std::vector<double> X;

    AdiabaticFun(const std::vector<double>& n_in_, double H_target_, const WgsConfig& cfg_)
        : n_in(n_in_), H_target(H_target_), cfg(cfg_),
          n(n_in_.size()), X(n_in_.size()) {}

    double f(double T) const
    {
        // Inner: solve ξ(T)
        double xi = solve_wgs_isothermal(n_in, T, cfg);

        composition_from_xi(n_in, xi, cfg, n);
        molefractions(n, X);

        double H = h(T, X);
        return H - H_target;
    }
};

// Build WgsConfig from species names
static WgsConfig make_wgs_config()
{
    WgsConfig cfg;
    cfg.i_CO  = species_index.at("CO");
    cfg.i_H2O = species_index.at("H2O");
    cfg.i_CO2 = species_index.at("CO2");
    cfg.i_H2  = species_index.at("H2");
    return cfg;
}

// Solve for adiabatic T with WGS partial equilibrium (internal)
static double solve_adiabatic_T_wgs(const std::vector<double>& n_in,
                                    double H_target,
                                    double T_guess,
                                    const WgsConfig& cfg)
{
    AdiabaticFun F{n_in, H_target, cfg};

    // Bracket temperature physically
    double Tmin = 300.0;
    double Tmax = 4500.0;

    return newton_damped(F, T_guess, Tmin, Tmax, 50);
}

// -------------------------------------------------------------
// SMR + WGS configuration (species indices)
// -------------------------------------------------------------

struct SmrWgsConfig {
    std::size_t i_CH4;
    std::size_t i_H2O;
    std::size_t i_CO;
    std::size_t i_CO2;
    std::size_t i_H2;
};

static SmrWgsConfig make_smr_wgs_config()
{
    SmrWgsConfig cfg;
    cfg.i_CH4 = species_index.at("CH4");
    cfg.i_H2O = species_index.at("H2O");
    cfg.i_CO  = species_index.at("CO");
    cfg.i_CO2 = species_index.at("CO2");
    cfg.i_H2  = species_index.at("H2");
    return cfg;
}

// Apply SMR and WGS extents to composition
// SMR: CH4 + H2O -> CO + 3H2  (extent xi1)
// WGS: CO + H2O -> CO2 + H2   (extent xi2)
static inline void composition_from_xi_smr_wgs(const std::vector<double>& n0,
                                                double xi1, double xi2,
                                                const SmrWgsConfig& cfg,
                                                std::vector<double>& n)
{
    n = n0;
    n[cfg.i_CH4] -= xi1;
    n[cfg.i_H2O] -= xi1 + xi2;
    n[cfg.i_CO]  += xi1 - xi2;
    n[cfg.i_CO2] += xi2;
    n[cfg.i_H2]  += 3.0 * xi1 + xi2;
}

// Equilibrium constant for SMR: CH4 + H2O <-> CO + 3H2
// Kp = (P_CO * P_H2^3) / (P_CH4 * P_H2O) * (P0/P)^2
// For mole fractions at pressure P: Kp = (y_CO * y_H2^3) / (y_CH4 * y_H2O) * (P/P0)^2
static inline double Kp_SMR(double T, const SmrWgsConfig& cfg)
{
    double g_CO  = ::g_over_RT(cfg.i_CO,  T);
    double g_H2  = ::g_over_RT(cfg.i_H2,  T);
    double g_CH4 = ::g_over_RT(cfg.i_CH4, T);
    double g_H2O = ::g_over_RT(cfg.i_H2O, T);
    
    // ΔG/RT = (g_CO + 3*g_H2) - (g_CH4 + g_H2O)
    double d_gRT = (g_CO + 3.0 * g_H2) - (g_CH4 + g_H2O);
    return std::exp(-d_gRT);
}

// Equilibrium constant for WGS (reuse existing)
static inline double Kp_WGS_smr(double T, const SmrWgsConfig& cfg)
{
    double g_CO2 = ::g_over_RT(cfg.i_CO2, T);
    double g_H2  = ::g_over_RT(cfg.i_H2,  T);
    double g_CO  = ::g_over_RT(cfg.i_CO,  T);
    double g_H2O = ::g_over_RT(cfg.i_H2O, T);
    
    double d_gRT = (g_CO2 + g_H2) - (g_CO + g_H2O);
    return std::exp(-d_gRT);
}

// -------------------------------------------------------------
// Solve SMR+WGS equilibrium at fixed T using 2D Newton
// -------------------------------------------------------------

// Compute residuals for SMR and WGS equilibrium
// F1 = ln(Q_SMR) - ln(Kp_SMR) = 0
// F2 = ln(Q_WGS) - ln(Kp_WGS) = 0
static void smr_wgs_residuals(const std::vector<double>& n0,
                              double xi1, double xi2,
                              double T, double P,
                              const SmrWgsConfig& cfg,
                              double& F1, double& F2)
{
    std::vector<double> n;
    composition_from_xi_smr_wgs(n0, xi1, xi2, cfg, n);
    
    double nt = 0.0;
    for (double v : n) nt += v;
    
    // Use a small but not tiny value to avoid numerical issues
    const double tiny = 1e-15;
    double y_CH4 = std::max(n[cfg.i_CH4] / nt, tiny);
    double y_H2O = std::max(n[cfg.i_H2O] / nt, tiny);
    double y_CO  = std::max(n[cfg.i_CO]  / nt, tiny);
    double y_CO2 = std::max(n[cfg.i_CO2] / nt, tiny);
    double y_H2  = std::max(n[cfg.i_H2]  / nt, tiny);
    
    // Reference pressure (1 atm = 101325 Pa)
    const double P0 = 101325.0;
    
    // SMR: Kp = (y_CO * y_H2^3) / (y_CH4 * y_H2O) * (P/P0)^2
    // Note: Δn = (1 + 3) - (1 + 1) = 2, so pressure factor is (P/P0)^2
    double Q_SMR = (y_CO * y_H2 * y_H2 * y_H2) / (y_CH4 * y_H2O) * (P / P0) * (P / P0);
    double Kp_smr = Kp_SMR(T, cfg);
    F1 = std::log(Q_SMR) - std::log(Kp_smr);
    
    // WGS: Kp = (y_CO2 * y_H2) / (y_CO * y_H2O)
    // Note: Δn = 0, no pressure dependence
    double Q_WGS = (y_CO2 * y_H2) / (y_CO * y_H2O);
    double Kp_wgs = Kp_WGS_smr(T, cfg);
    F2 = std::log(Q_WGS) - std::log(Kp_wgs);
}

// 2D Newton solver for SMR+WGS equilibrium at fixed T
static bool solve_smr_wgs_isothermal(const std::vector<double>& n0,
                                     double T, double P,
                                     const SmrWgsConfig& cfg,
                                     double& xi1_out, double& xi2_out,
                                     std::size_t maxIt = 100)
{
    // Check if SMR can proceed (need CH4 and H2O)
    const double min_species = 1e-12;
    bool can_smr = (n0[cfg.i_CH4] > min_species && n0[cfg.i_H2O] > min_species);
    
    // Check if WGS can proceed (need either CO+H2O or CO2+H2)
    bool can_wgs_fwd = (n0[cfg.i_CO] > min_species && n0[cfg.i_H2O] > min_species);
    bool can_wgs_rev = (n0[cfg.i_CO2] > min_species && n0[cfg.i_H2] > min_species);
    bool can_wgs = can_wgs_fwd || can_wgs_rev;
    
    // If neither reaction can proceed, return zero extents
    if (!can_smr && !can_wgs) {
        xi1_out = 0.0;
        xi2_out = 0.0;
        return true;
    }
    
    // Initial guess: small extents
    double xi1 = can_smr ? std::min(0.001, n0[cfg.i_CH4] * 0.1) : 0.0;
    double xi2 = 0.0;
    
    // Determine bounds from non-negativity
    // n[CH4] = n0[CH4] - xi1 >= 0           => xi1 <= n0[CH4]
    // n[H2O] = n0[H2O] - xi1 - xi2 >= 0     => xi1 + xi2 <= n0[H2O]
    // n[CO]  = n0[CO] + xi1 - xi2 >= 0      => xi2 <= n0[CO] + xi1
    // n[CO2] = n0[CO2] + xi2 >= 0           => xi2 >= -n0[CO2]
    // n[H2]  = n0[H2] + 3*xi1 + xi2 >= 0    => 3*xi1 + xi2 >= -n0[H2]
    
    const double tol = 1e-8;
    const double eps = 1e-8;
    
    for (std::size_t it = 0; it < maxIt; ++it) {
        double F1, F2;
        smr_wgs_residuals(n0, xi1, xi2, T, P, cfg, F1, F2);
        
        // Check convergence
        if (std::abs(F1) < tol && std::abs(F2) < tol) {
            xi1_out = xi1;
            xi2_out = xi2;
            return true;
        }
        
        // Compute Jacobian by finite differences
        double F1_dxi1, F2_dxi1, F1_dxi2, F2_dxi2;
        
        double F1p, F2p;
        smr_wgs_residuals(n0, xi1 + eps, xi2, T, P, cfg, F1p, F2p);
        F1_dxi1 = (F1p - F1) / eps;
        F2_dxi1 = (F2p - F2) / eps;
        
        smr_wgs_residuals(n0, xi1, xi2 + eps, T, P, cfg, F1p, F2p);
        F1_dxi2 = (F1p - F1) / eps;
        F2_dxi2 = (F2p - F2) / eps;
        
        // Solve 2x2 system: J * dx = -F
        double det = F1_dxi1 * F2_dxi2 - F1_dxi2 * F2_dxi1;
        if (std::abs(det) < 1e-20) {
            // Singular Jacobian - try small perturbation
            xi1 += 1e-6;
            continue;
        }
        
        double dxi1 = (-F1 * F2_dxi2 + F2 * F1_dxi2) / det;
        double dxi2 = (-F2 * F1_dxi1 + F1 * F2_dxi1) / det;
        
        // Line search with bounds enforcement
        double lambda = 1.0;
        for (int ls = 0; ls < 20; ++ls) {
            double xi1_new = xi1 + lambda * dxi1;
            double xi2_new = xi2 + lambda * dxi2;
            
            // Check bounds
            bool valid = true;
            if (xi1_new < -1e-12) valid = false;  // xi1 >= 0 (can't un-reform)
            if (xi1_new > n0[cfg.i_CH4] + 1e-12) valid = false;
            if (xi1_new + xi2_new > n0[cfg.i_H2O] + 1e-12) valid = false;
            if (xi2_new < -n0[cfg.i_CO2] - 1e-12) valid = false;
            if (n0[cfg.i_CO] + xi1_new - xi2_new < -1e-12) valid = false;
            if (n0[cfg.i_H2] + 3.0 * xi1_new + xi2_new < -1e-12) valid = false;
            
            if (valid) {
                // Check if residual decreased
                double F1_new, F2_new;
                smr_wgs_residuals(n0, xi1_new, xi2_new, T, P, cfg, F1_new, F2_new);
                double norm_old = F1 * F1 + F2 * F2;
                double norm_new = F1_new * F1_new + F2_new * F2_new;
                
                if (norm_new < norm_old || lambda < 1e-4) {
                    xi1 = std::max(0.0, xi1_new);
                    xi2 = xi2_new;
                    break;
                }
            }
            lambda *= 0.5;
        }
    }
    
    // Return best guess even if not fully converged
    xi1_out = xi1;
    xi2_out = xi2;
    return false;
}

// -------------------------------------------------------------
// Adiabatic SMR+WGS solver
// -------------------------------------------------------------

struct SmrWgsAdiabaticFun {
    const std::vector<double>& n_in;
    double H_target;
    double P;
    SmrWgsConfig cfg;
    
    mutable std::vector<double> n;
    mutable std::vector<double> X;
    
    SmrWgsAdiabaticFun(const std::vector<double>& n_in_, double H_target_, 
                       double P_, const SmrWgsConfig& cfg_)
        : n_in(n_in_), H_target(H_target_), P(P_), cfg(cfg_),
          n(n_in_.size()), X(n_in_.size()) {}
    
    double f(double T) const
    {
        // Inner: solve equilibrium composition at T
        double xi1, xi2;
        solve_smr_wgs_isothermal(n_in, T, P, cfg, xi1, xi2);
        
        composition_from_xi_smr_wgs(n_in, xi1, xi2, cfg, n);
        molefractions(n, X);
        
        double H = h(T, X);
        return H - H_target;
    }
};

static double solve_adiabatic_T_smr_wgs(const std::vector<double>& n_in,
                                        double H_target,
                                        double P,
                                        double T_guess,
                                        const SmrWgsConfig& cfg)
{
    SmrWgsAdiabaticFun F{n_in, H_target, P, cfg};
    
    // Bracket temperature physically
    double Tmin = 300.0;
    double Tmax = 4500.0;
    
    return newton_damped(F, T_guess, Tmin, Tmax, 50);
}

} // anonymous namespace

// -------------------------------------------------------------
// State-based WGS equilibrium functions (Public API)
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
State wgs_equilibrium(const State& in)
{
    WgsConfig cfg = make_wgs_config();
    
    // Solve for extent of reaction at fixed T
    WgsIsothermalFun F{in.X, in.T, cfg};
    
    // Bounds from non-negativity (same logic as solve_wgs_isothermal)
    double xi_min = std::max(-in.X[cfg.i_CO2], -in.X[cfg.i_H2]);
    double xi_max = std::min(in.X[cfg.i_CO], in.X[cfg.i_H2O]);
    
    double xi = newton_damped(F, 0.0, xi_min, xi_max, 40);
    
    // Build output composition
    std::vector<double> n_out;
    composition_from_xi(in.X, xi, cfg, n_out);
    
    std::vector<double> X_out;
    molefractions(n_out, X_out);
    
    State out;
    out.T = in.T;
    out.P = in.P;
    out.X = X_out;
    return out;
}

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
State wgs_equilibrium_adiabatic(const State& in)
{
    WgsConfig cfg = make_wgs_config();
    
    // Target enthalpy (conserved)
    double H_target = h(in.T, in.X);
    
    // Solve for adiabatic temperature
    double T_ad = solve_adiabatic_T_wgs(in.X, H_target, in.T, cfg);
    
    // Get equilibrium composition at T_ad
    State at_T_ad;
    at_T_ad.T = T_ad;
    at_T_ad.P = in.P;
    at_T_ad.X = in.X;
    
    State out = wgs_equilibrium(at_T_ad);
    return out;
}

// -------------------------------------------------------------
// State-based SMR+WGS equilibrium functions (Public API)
// -------------------------------------------------------------

// SMR+WGS equilibrium (isothermal) - equilibrate at input temperature
State smr_wgs_equilibrium(const State& in)
{
    SmrWgsConfig cfg = make_smr_wgs_config();
    
    // Check if SMR can proceed (need CH4)
    // If no CH4, SMR+WGS reduces to just WGS
    const double min_ch4 = 1e-10;
    if (in.X[cfg.i_CH4] < min_ch4) {
        // No CH4 to reform - fall back to WGS-only equilibrium
        return wgs_equilibrium(in);
    }
    
    // Solve for extents of reaction at fixed T
    double xi1, xi2;
    solve_smr_wgs_isothermal(in.X, in.T, in.P, cfg, xi1, xi2);
    
    // Build output composition
    std::vector<double> n_out;
    composition_from_xi_smr_wgs(in.X, xi1, xi2, cfg, n_out);
    
    std::vector<double> X_out;
    molefractions(n_out, X_out);
    
    State out;
    out.T = in.T;
    out.P = in.P;
    out.X = X_out;
    return out;
}

// SMR+WGS equilibrium (adiabatic) - solve for equilibrium T and composition
State smr_wgs_equilibrium_adiabatic(const State& in)
{
    SmrWgsConfig cfg = make_smr_wgs_config();
    
    // Check if SMR can proceed (need CH4)
    // If no CH4, SMR+WGS reduces to just WGS
    const double min_ch4 = 1e-10;
    if (in.X[cfg.i_CH4] < min_ch4) {
        // No CH4 to reform - fall back to WGS-only equilibrium
        return wgs_equilibrium_adiabatic(in);
    }
    
    // Target enthalpy (conserved)
    double H_target = h(in.T, in.X);
    
    // Solve for adiabatic temperature
    double T_ad = solve_adiabatic_T_smr_wgs(in.X, H_target, in.P, in.T, cfg);
    
    // Get equilibrium composition at T_ad
    double xi1, xi2;
    solve_smr_wgs_isothermal(in.X, T_ad, in.P, cfg, xi1, xi2);
    
    std::vector<double> n_out;
    composition_from_xi_smr_wgs(in.X, xi1, xi2, cfg, n_out);
    
    std::vector<double> X_out;
    molefractions(n_out, X_out);
    
    State out;
    out.T = T_ad;
    out.P = in.P;
    out.X = X_out;
    return out;
}