#include "../include/equilibrium.h"

// -------------------------------------------------------------
// Internal utilities
// -------------------------------------------------------------

namespace {

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

    // Determine ξ bounds
    double xi_min = -1e30;
    double xi_max =  1e30;

    // From n[k] >= 0 → xi >= -n0[k]/nu
    xi_min = std::max(xi_min, -n0[cfg.i_CO]);
    xi_min = std::max(xi_min, -n0[cfg.i_H2O]);
    // From n[k] >= 0 → xi <= n0[k] for products
    xi_max = std::min(xi_max, n0[cfg.i_CO2]);
    xi_max = std::min(xi_max, n0[cfg.i_H2]);

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

} // anonymous namespace

// -------------------------------------------------------------
// Public API (legacy)
// -------------------------------------------------------------

double solve_adiabatic_T_wgs(const std::vector<double>& n_in,
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
// State-based WGS equilibrium functions
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
State wgs_equilibrium(const State& in)
{
    WgsConfig cfg = make_wgs_config();
    
    // Solve for extent of reaction at fixed T
    WgsIsothermalFun F{in.X, in.T, cfg};
    
    double xi_min = -1e30;
    double xi_max =  1e30;
    xi_min = std::max(xi_min, -in.X[cfg.i_CO]);
    xi_min = std::max(xi_min, -in.X[cfg.i_H2O]);
    xi_max = std::min(xi_max, in.X[cfg.i_CO2]);
    xi_max = std::min(xi_max, in.X[cfg.i_H2]);
    
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