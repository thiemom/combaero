#pragma once

#include "thermo_transport.h"
#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace combaero::equilibrium {

struct WgsConfig {
    std::size_t i_CO;
    std::size_t i_H2O;
    std::size_t i_CO2;
    std::size_t i_H2;
};

// must be provided by you, using thermo_transport_data.h + NASA
double g_over_RT(std::size_t k, double T);

// ----- small helpers -----

inline void composition_from_xi(const std::vector<double>& n0,
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

inline void normalize_to_X(const std::vector<double>& n,
                           std::vector<double>& X)
{
    double nt = 0.0;
    for (double v : n) nt += v;

    X.resize(n.size());
    if (nt <= 0.0) {
        double inv = 1.0 / static_cast<double>(n.size());
        for (std::size_t k = 0; k < n.size(); ++k) X[k] = inv;
    } else {
        for (std::size_t k = 0; k < n.size(); ++k) X[k] = n[k] / nt;
    }
}

// WGS Kp(T) from Gibbs functions
inline double Kp_WGS(double T, const WgsConfig& cfg)
{
    double gRT_CO2 = g_over_RT(cfg.i_CO2, T);
    double gRT_H2  = g_over_RT(cfg.i_H2,  T);
    double gRT_CO  = g_over_RT(cfg.i_CO,  T);
    double gRT_H2O = g_over_RT(cfg.i_H2O, T);

    double d_gRT = (gRT_CO2 + gRT_H2) - (gRT_CO + gRT_H2O);
    return std::exp(-d_gRT);
}

// ----- generic 1D Newton with finite-diff derivative -----

template <class Fun>
double newton1d(Fun&& F,
                double x0,
                double tol_f   = 1e-10,
                double tol_dx  = 1e-8,
                int    max_it  = 50)
{
    double x = x0;

    for (int it = 0; it < max_it; ++it) {
        double f0 = F.f(x);
        if (std::abs(f0) < tol_f) return x;

        double eps = 1e-6 * std::max(1.0, std::abs(x));
        double f1  = F.f(x + eps);
        double df  = (f1 - f0) / eps;
        if (df == 0.0) break;

        double dx = -f0 / df;
        x += dx;

        if (std::abs(dx) < tol_dx * std::max(1.0, std::abs(x))) return x;
    }
    return x; // let caller decide if convergence is good enough
}

// ----- inner: WGS equilibrium at fixed T -----

struct WgsIsothermalFun {
    const std::vector<double>& n0;
    WgsConfig cfg;
    double lnKp;
    mutable std::vector<double> n;

    WgsIsothermalFun(const std::vector<double>& n0_,
                     double T,
                     const WgsConfig& cfg_)
        : n0(n0_), cfg(cfg_), n(n0_.size())
    {
        lnKp = std::log(Kp_WGS(T, cfg));
    }

    double f(double xi) const {
        composition_from_xi(n0, xi, cfg, n);

        double nt = 0.0;
        for (double v : n) nt += v;

        double yCO2 = n[cfg.i_CO2] / nt;
        double yH2  = n[cfg.i_H2]  / nt;
        double yCO  = n[cfg.i_CO]  / nt;
        double yH2O = n[cfg.i_H2O] / nt;

        double Q = (yCO2 * yH2) / (yCO * yH2O);
        return std::log(Q) - lnKp;
    }
};

inline double solve_wgs_isothermal(const std::vector<double>& n0,
                                   double T,
                                   const WgsConfig& cfg)
{
    WgsIsothermalFun F{n0, T, cfg};
    return newton1d(F, 0.0);
}

// ----- outer: adiabatic constraint H(T, ξ(T)) = H_target -----

struct AdiabaticWgsFun {
    const std::vector<double>& n_in;
    double H_target;
    WgsConfig cfg;
    mutable std::vector<double> n;
    mutable std::vector<double> X;

    AdiabaticWgsFun(const std::vector<double>& n_in_,
                    double H_target_,
                    const WgsConfig& cfg_)
        : n_in(n_in_),
          H_target(H_target_),
          cfg(cfg_),
          n(n_in_.size()),
          X(n_in_.size())
    {}

    double f(double T) const {
        double xi = solve_wgs_isothermal(n_in, T, cfg);
        composition_from_xi(n_in, xi, cfg, n);
        normalize_to_X(n, X);

        double H = h(T, X);      // uses your existing mixture h(T, X)
        return H - H_target;
    }
};

// main entry: solve for adiabatic T with WGS partial equilibrium
inline double solve_adiabatic_T_wgs(const std::vector<double>& n_in,
                                    double H_target,
                                    double T_guess,
                                    const WgsConfig& cfg)
{
    AdiabaticWgsFun F{n_in, H_target, cfg};
    return newton1d(F, T_guess);
}

} // namespace combaero::equilibrium