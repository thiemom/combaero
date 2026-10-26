#include "../include/equilibrium.h"
#include "../include/combustion.h"
#include "thermo.h"
#include "composition.h"
#include "thermo_transport_data.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace combaero {

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

static inline double clamp(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

static inline void composition_from_xi(const std::vector<double> &n0, double xi,
                                       const WgsConfig &cfg,
                                       std::vector<double> &n) {
  n = n0;
  n[cfg.i_CO] -= xi;
  n[cfg.i_H2O] -= xi;
  n[cfg.i_CO2] += xi;
  n[cfg.i_H2] += xi;
}

static inline void molefractions(const std::vector<double> &n,
                                 std::vector<double> &X) {
  double nt = 0.0;
  for (double v : n)
    nt += v;

  X.resize(n.size());
  if (nt <= 0.0) {
    double inv = 1.0 / static_cast<double>(n.size());
    for (std::size_t k = 0; k < n.size(); ++k)
      X[k] = inv;
  } else {
    for (std::size_t k = 0; k < n.size(); ++k)
      X[k] = n[k] / nt;
  }
}

static inline double kp_wgs_from_indices(double T, std::size_t i_CO2,
                                         std::size_t i_H2, std::size_t i_CO,
                                         std::size_t i_H2O) {
  double d_gRT = (g_over_RT(i_CO2, T) + g_over_RT(i_H2, T)) -
                 (g_over_RT(i_CO, T) + g_over_RT(i_H2O, T));
  return std::exp(-d_gRT);
}

static inline double Kp_WGS(double T, const WgsConfig &cfg) {
  return kp_wgs_from_indices(T, cfg.i_CO2, cfg.i_H2, cfg.i_CO, cfg.i_H2O);
}

// -------------------------------------------------------------
// Newton with damping
// -------------------------------------------------------------

template <typename Fun>
static double newton_damped(Fun &&F, double x0, double xmin, double xmax,
                            std::size_t maxIt = 40) {
  double x = x0;

  for (std::size_t it = 0; it < maxIt; ++it) {
    double f0 = F.f(x);

    // Converged?
    if (std::abs(f0) < 1e-12)
      return clamp(x, xmin, xmax);

    // Finite diff derivative
    double eps = 1e-6 * std::max(1.0, std::abs(x));
    double f1 = F.f(x + eps);
    double df = (f1 - f0) / eps;

    if (std::abs(df) < 1e-14)
      return clamp(x, xmin, xmax);

    double dx = -f0 / df;
    double x_trial = x + dx;

    // Damping if outside bounds
    double lambda = 1.0;
    while ((x_trial < xmin || x_trial > xmax) && lambda > 1e-6) {
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
  const std::vector<double> &n0;
  WgsConfig cfg;
  double lnKp;
  mutable std::vector<double> n;

  WgsIsothermalFun(const std::vector<double> &n0_, double T,
                   const WgsConfig &cfg_)
      : n0(n0_), cfg(cfg_), n(n0_.size()) {
    lnKp = std::log(Kp_WGS(T, cfg));
  }

  double f(double xi) const {
    composition_from_xi(n0, xi, cfg, n);

    double nt = 0.0;
    for (double v : n)
      nt += v;

    // Avoid divisions by zero
    const double tiny = 1e-300;
    double yCO2 = std::max(n[cfg.i_CO2] / nt, tiny);
    double yH2 = std::max(n[cfg.i_H2] / nt, tiny);
    double yCO = std::max(n[cfg.i_CO] / nt, tiny);
    double yH2O = std::max(n[cfg.i_H2O] / nt, tiny);

    double Q = (yCO2 * yH2) / (yCO * yH2O);
    return std::log(Q) - lnKp;
  }
};

static double solve_wgs_isothermal(const std::vector<double> &n0, double T,
                                   const WgsConfig &cfg) {
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
  const std::vector<double> &n_in;
  double H_mass_target;
  WgsConfig cfg;

  mutable std::vector<double> n;
  mutable std::vector<double> X;

  AdiabaticFun(const std::vector<double> &n_in_, double H_mass_target_,
               const WgsConfig &cfg_)
      : n_in(n_in_), H_mass_target(H_mass_target_), cfg(cfg_), n(n_in_.size()),
        X(n_in_.size()) {}

  double f(double T) const {
    // Inner: solve ξ(T)
    double xi = solve_wgs_isothermal(n_in, T, cfg);

    composition_from_xi(n_in, xi, cfg, n);
    molefractions(n, X);

    double H_mass = h(T, X) / mwmix(X);
    return H_mass - H_mass_target;
  }
};

// Build WgsConfig from species names
static WgsConfig make_wgs_config() {
  WgsConfig cfg;
  cfg.i_CO = species_index.at("CO");
  cfg.i_H2O = species_index.at("H2O");
  cfg.i_CO2 = species_index.at("CO2");
  cfg.i_H2 = species_index.at("H2");
  return cfg;
}

// Solve for adiabatic T with WGS partial equilibrium (internal)
static double solve_adiabatic_T_wgs(const std::vector<double> &n_in,
                                    double H_mass_target, double T_guess,
                                    const WgsConfig &cfg) {
  AdiabaticFun F{n_in, H_mass_target, cfg};

  // Bracket temperature physically
  double Tmin = 300.0;
  double Tmax = 4500.0;

  return newton_damped(F, T_guess, Tmin, Tmax, 50);
}

// -------------------------------------------------------------
// General Steam Reforming + WGS configuration
// -------------------------------------------------------------
// Handles all hydrocarbons CnHm via:
//   CnHm + n*H2O -> n*CO + (n + m/2)*H2
// Plus WGS:
//   CO + H2O <-> CO2 + H2

struct ReformingConfig {
  std::size_t i_H2O;
  std::size_t i_CO;
  std::size_t i_CO2;
  std::size_t i_H2;

  // Hydrocarbon species indices and their C, H counts
  std::vector<std::size_t> hc_indices;
  std::vector<int> hc_C; // number of C atoms
  std::vector<int> hc_H; // number of H atoms
};

static ReformingConfig make_reforming_config() {
  ReformingConfig cfg;
  cfg.i_H2O = species_index.at("H2O");
  cfg.i_CO = species_index.at("CO");
  cfg.i_CO2 = species_index.at("CO2");
  cfg.i_H2 = species_index.at("H2");

  // Find all hydrocarbon species (C > 0, H > 0, O = 0, N = 0)
  for (std::size_t i = 0; i < species_names.size(); ++i) {
    const auto &ms = molecular_structures[i];
    if (ms.C > 0 && ms.H > 0 && ms.O == 0 && ms.N == 0) {
      cfg.hc_indices.push_back(i);
      cfg.hc_C.push_back(static_cast<int>(ms.C));
      cfg.hc_H.push_back(static_cast<int>(ms.H));
    }
  }
  return cfg;
}

// Legacy config for backward compatibility
struct SmrWgsConfig {
  std::size_t i_CH4;
  std::size_t i_H2O;
  std::size_t i_CO;
  std::size_t i_CO2;
  std::size_t i_H2;
};

static SmrWgsConfig make_smr_wgs_config() {
  SmrWgsConfig cfg;
  cfg.i_CH4 = species_index.at("CH4");
  cfg.i_H2O = species_index.at("H2O");
  cfg.i_CO = species_index.at("CO");
  cfg.i_CO2 = species_index.at("CO2");
  cfg.i_H2 = species_index.at("H2");
  return cfg;
}

// Apply SMR and WGS extents to composition
// SMR: CH4 + H2O -> CO + 3H2  (extent xi1)
// WGS: CO + H2O -> CO2 + H2   (extent xi2)
static inline void composition_from_xi_smr_wgs(const std::vector<double> &n0,
                                               double xi1, double xi2,
                                               const SmrWgsConfig &cfg,
                                               std::vector<double> &n) {
  n = n0;
  n[cfg.i_CH4] -= xi1;
  n[cfg.i_H2O] -= xi1 + xi2;
  n[cfg.i_CO] += xi1 - xi2;
  n[cfg.i_CO2] += xi2;
  n[cfg.i_H2] += 3.0 * xi1 + xi2;
}

// Equilibrium constant for SMR: CH4 + H2O <-> CO + 3H2
// Kp = (P_CO * P_H2^3) / (P_CH4 * P_H2O) * (P0/P)^2
// For mole fractions at pressure P: Kp = (y_CO * y_H2^3) / (y_CH4 * y_H2O) *
// (P/P0)^2
static inline double Kp_SMR(double T, const SmrWgsConfig &cfg) {
  double g_CO = g_over_RT(cfg.i_CO, T);
  double g_H2 = g_over_RT(cfg.i_H2, T);
  double g_CH4 = g_over_RT(cfg.i_CH4, T);
  double g_H2O = g_over_RT(cfg.i_H2O, T);

  // ΔG/RT = (g_CO + 3*g_H2) - (g_CH4 + g_H2O)
  double d_gRT = (g_CO + 3.0 * g_H2) - (g_CH4 + g_H2O);
  return std::exp(-d_gRT);
}

// Equilibrium constant for WGS — delegates to shared helper
static inline double Kp_WGS_smr(double T, const SmrWgsConfig &cfg) {
  return kp_wgs_from_indices(T, cfg.i_CO2, cfg.i_H2, cfg.i_CO, cfg.i_H2O);
}

// -------------------------------------------------------------
// Solve SMR+WGS equilibrium at fixed T using 2D Newton
// -------------------------------------------------------------

// Compute residuals for SMR and WGS equilibrium
// F1 = ln(Q_SMR) - ln(Kp_SMR) = 0
// F2 = ln(Q_WGS) - ln(Kp_WGS) = 0
static void smr_wgs_residuals(const std::vector<double> &n0, double xi1,
                              double xi2, double T, double P,
                              const SmrWgsConfig &cfg, double &F1, double &F2) {
  std::vector<double> n;
  composition_from_xi_smr_wgs(n0, xi1, xi2, cfg, n);

  double nt = 0.0;
  for (double v : n)
    nt += v;

  // Use a small but not tiny value to avoid numerical issues
  const double tiny = 1e-15;
  double y_CH4 = std::max(n[cfg.i_CH4] / nt, tiny);
  double y_H2O = std::max(n[cfg.i_H2O] / nt, tiny);
  double y_CO = std::max(n[cfg.i_CO] / nt, tiny);
  double y_CO2 = std::max(n[cfg.i_CO2] / nt, tiny);
  double y_H2 = std::max(n[cfg.i_H2] / nt, tiny);

  // Reference pressure (1 atm = 101325 Pa)
  const double P0 = 101325.0;

  // SMR: Kp = (y_CO * y_H2^3) / (y_CH4 * y_H2O) * (P/P0)^2
  // Note: Δn = (1 + 3) - (1 + 1) = 2, so pressure factor is (P/P0)^2
  double Q_SMR =
      (y_CO * y_H2 * y_H2 * y_H2) / (y_CH4 * y_H2O) * (P / P0) * (P / P0);
  double Kp_smr = Kp_SMR(T, cfg);
  F1 = std::log(Q_SMR) - std::log(Kp_smr);

  // WGS: Kp = (y_CO2 * y_H2) / (y_CO * y_H2O)
  // Note: Δn = 0, no pressure dependence
  double Q_WGS = (y_CO2 * y_H2) / (y_CO * y_H2O);
  double Kp_wgs = Kp_WGS_smr(T, cfg);
  F2 = std::log(Q_WGS) - std::log(Kp_wgs);
}

// 2D Newton solver for SMR+WGS equilibrium at fixed T
static bool solve_smr_wgs_isothermal(const std::vector<double> &n0, double T,
                                     double P, const SmrWgsConfig &cfg,
                                     double &xi1_out, double &xi2_out,
                                     std::size_t maxIt = 100) {
  // Check if SMR can proceed (need CH4 and H2O)
  const double min_species = 1e-12;
  bool can_smr = (n0[cfg.i_CH4] > min_species && n0[cfg.i_H2O] > min_species);

  // Check if WGS can proceed (need either CO+H2O or CO2+H2)
  bool can_wgs_fwd =
      (n0[cfg.i_CO] > min_species && n0[cfg.i_H2O] > min_species);
  bool can_wgs_rev =
      (n0[cfg.i_CO2] > min_species && n0[cfg.i_H2] > min_species);
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
      if (xi1_new < -1e-12)
        valid = false; // xi1 >= 0 (can't un-reform)
      if (xi1_new > n0[cfg.i_CH4] + 1e-12)
        valid = false;
      if (xi1_new + xi2_new > n0[cfg.i_H2O] + 1e-12)
        valid = false;
      if (xi2_new < -n0[cfg.i_CO2] - 1e-12)
        valid = false;
      if (n0[cfg.i_CO] + xi1_new - xi2_new < -1e-12)
        valid = false;
      if (n0[cfg.i_H2] + 3.0 * xi1_new + xi2_new < -1e-12)
        valid = false;

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
  const std::vector<double> &n_in;
  double H_mass_target;
  double P;
  SmrWgsConfig cfg;

  mutable std::vector<double> n;
  mutable std::vector<double> X;

  SmrWgsAdiabaticFun(const std::vector<double> &n_in_, double H_mass_target_,
                     double P_, const SmrWgsConfig &cfg_)
      : n_in(n_in_), H_mass_target(H_mass_target_), P(P_), cfg(cfg_),
        n(n_in_.size()), X(n_in_.size()) {}

  double f(double T) const {
    // Inner: solve equilibrium composition at T
    double xi1, xi2;
    solve_smr_wgs_isothermal(n_in, T, P, cfg, xi1, xi2);

    composition_from_xi_smr_wgs(n_in, xi1, xi2, cfg, n);
    molefractions(n, X);

    double H_mass = h(T, X) / mwmix(X);
    return H_mass - H_mass_target;
  }
};

static double solve_adiabatic_T_smr_wgs(const std::vector<double> &n_in,
                                        double H_mass_target, double P,
                                        double T_guess,
                                        const SmrWgsConfig &cfg) {
  SmrWgsAdiabaticFun F{n_in, H_mass_target, P, cfg};

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
EquilibriumResult wgs_equilibrium(const State &in) {
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
  out.set_TPX(in.T, in.P, X_out);
  return EquilibriumResult{out, in.T, 0.0, true};
}

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
EquilibriumResult wgs_equilibrium_adiabatic(const State &in) {
  WgsConfig cfg = make_wgs_config();

  // Target enthalpy (conserved)
  double H_mass_target = h(in.T, in.X) / mwmix(in.X);

  // Solve for adiabatic temperature
  double T_ad = solve_adiabatic_T_wgs(in.X, H_mass_target, in.T, cfg);

  // Get equilibrium composition at T_ad
  State at_T_ad;
  at_T_ad.set_TPX(T_ad, in.P, in.X);

  EquilibriumResult inner = wgs_equilibrium(at_T_ad);
  return EquilibriumResult{inner.state, in.T, inner.state.T - in.T, true};
}

// -------------------------------------------------------------
// State-based SMR+WGS equilibrium functions (Public API)
// -------------------------------------------------------------

// SMR+WGS equilibrium (isothermal) - equilibrate at input temperature
EquilibriumResult smr_wgs_equilibrium(const State &in) {
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
  bool conv = solve_smr_wgs_isothermal(in.X, in.T, in.P, cfg, xi1, xi2);

  // Build output composition
  std::vector<double> n_out;
  composition_from_xi_smr_wgs(in.X, xi1, xi2, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(in.T, in.P, X_out);
  return EquilibriumResult{out, in.T, 0.0, conv};
}

// SMR+WGS equilibrium (adiabatic) - solve for equilibrium T and composition
EquilibriumResult smr_wgs_equilibrium_adiabatic(const State &in) {
  SmrWgsConfig cfg = make_smr_wgs_config();

  // Check if SMR can proceed (need CH4)
  // If no CH4, SMR+WGS reduces to just WGS
  const double min_ch4 = 1e-10;
  if (in.X[cfg.i_CH4] < min_ch4) {
    // No CH4 to reform - fall back to WGS-only equilibrium
    return wgs_equilibrium_adiabatic(in);
  }

  // Target enthalpy (conserved)
  double H_mass_target = h(in.T, in.X) / mwmix(in.X);

  // Solve for adiabatic temperature
  double T_ad = solve_adiabatic_T_smr_wgs(in.X, H_mass_target, in.P, in.T, cfg);

  // Get equilibrium composition at T_ad
  double xi1, xi2;
  bool conv = solve_smr_wgs_isothermal(in.X, T_ad, in.P, cfg, xi1, xi2);

  std::vector<double> n_out;
  composition_from_xi_smr_wgs(in.X, xi1, xi2, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(T_ad, in.P, X_out);
  return EquilibriumResult{out, in.T, T_ad - in.T, conv};
}

// -------------------------------------------------------------
// General Reforming + WGS equilibrium (handles all hydrocarbons)
// -------------------------------------------------------------

// Equilibrium constant for general steam reforming: CnHm + n*H2O <-> n*CO +
// (n+m/2)*H2 Kp = (y_CO^n * y_H2^(n+m/2)) / (y_CnHm * y_H2O^n) *
// (P/P0)^(delta_n) where delta_n = n + (n+m/2) - 1 - n = m/2 - 1 + n = (m-2)/2
// + n
static double Kp_reforming(double T, std::size_t i_hc, int n_C, int n_H,
                           const ReformingConfig &cfg) {
  double g_hc = g_over_RT(i_hc, T);
  double g_H2O = g_over_RT(cfg.i_H2O, T);
  double g_CO = g_over_RT(cfg.i_CO, T);
  double g_H2 = g_over_RT(cfg.i_H2, T);

  // CnHm + n*H2O -> n*CO + (n + m/2)*H2
  double n_H2_prod = n_C + 0.5 * n_H; // n + m/2

  // ΔG/RT = n*g_CO + (n+m/2)*g_H2 - g_CnHm - n*g_H2O
  double d_gRT = n_C * g_CO + n_H2_prod * g_H2 - g_hc - n_C * g_H2O;
  return std::exp(-d_gRT);
}

// Solve single hydrocarbon reforming equilibrium at fixed T
// CnHm + n*H2O <-> n*CO + (n+m/2)*H2
static double solve_reforming_isothermal(const std::vector<double> &n0,
                                         double T, double P, std::size_t i_hc,
                                         int n_C, int n_H,
                                         const ReformingConfig &cfg) {
  const double min_species = 1e-12;
  const double m_H2 = n_C + 0.5 * n_H;

  // The extent is bracketed both ways: forward (reforming) until the
  // hydrocarbon or the steam runs out, backward (methanation, the hydrocarbon
  // re-forming from CO + H2) until CO or H2 does. Forward only left a mixture
  // of higher hydrocarbons at 600-800 K 0.2 off in mole fraction: at
  // equilibrium they become methane, which is CH4 reforming run backward.
  const double xi_hi = std::min(n0[i_hc], n0[cfg.i_H2O] / n_C);
  const double xi_lo = -std::min(n0[cfg.i_CO] / n_C, n0[cfg.i_H2] / m_H2);
  if (xi_hi - xi_lo < min_species) {
    return 0.0; // Nothing can move either way
  }

  double Kp = Kp_reforming(T, i_hc, n_C, n_H, cfg);
  double lnKp = std::log(Kp);

  // Reference pressure
  const double P0 = 101325.0;

  // delta_n = products - reactants = (n + (n+m/2)) - (1 + n) = n + m/2 - 1
  double delta_n = n_C + 0.5 * n_H - 1.0;

  const double span = xi_hi - xi_lo;

  // Safeguarded Newton on F(xi) = ln Q - ln Kp over (xi_lo, xi_hi), with the
  // analytic slope (#716)
  //   dF/dxi = nC^2/n_CO + m^2/n_H2 + 1/n_hc + nC^2/n_H2O - dn^2/n_t,  m = nC + nH/2,
  // and bisection where a step leaves the bracket. F rises from -inf (CO or H2
  // used up) to +inf as the hydrocarbon or the steam runs out, so the root is
  // bracketed. The former finite-difference slope used an absolute step of
  // 1e-8 mol: near complete reforming (1e-11 mol left) it stepped past the
  // hydrocarbon and stalled off the root.
  double nt0 = 0.0;
  for (double v : n0)
    nt0 += v;
  auto eval = [&](double xi, double &dF) {
    const double n_hc = n0[i_hc] - xi;
    const double n_H2O = n0[cfg.i_H2O] - n_C * xi;
    const double n_CO = n0[cfg.i_CO] + n_C * xi;
    const double n_H2 = n0[cfg.i_H2] + m_H2 * xi;
    const double nt = nt0 + delta_n * xi;
    const double lnQ = n_C * std::log(n_CO) + m_H2 * std::log(n_H2) - std::log(n_hc) -
                       n_C * std::log(n_H2O) - delta_n * std::log(nt) +
                       delta_n * std::log(P / P0);
    dF = n_C * n_C / n_CO + m_H2 * m_H2 / n_H2 + 1.0 / n_hc + n_C * n_C / n_H2O -
         delta_n * delta_n / nt;
    return lnQ - lnKp;
  };
  double lo = xi_lo;
  double hi = xi_hi;
  double xi = 0.5 * (xi_lo + xi_hi);
  for (int it = 0; it < 200; ++it) {
    double dF = 0.0;
    const double F = eval(xi, dF);
    if (!std::isfinite(F)) {
      break;
    }
    if (std::abs(F) < 1e-12) {
      return xi;
    }
    if (F > 0.0) {
      hi = xi;
    } else {
      lo = xi;
    }
    double xi_new = (dF > 0.0 && std::isfinite(dF)) ? xi - F / dF : 0.5 * (lo + hi);
    if (!(xi_new > lo && xi_new < hi)) {
      xi_new = 0.5 * (lo + hi);
    }
    if (hi - lo <= 1e-15 * span) {
      return xi_new;
    }
    xi = xi_new;
  }
  return std::max(xi_lo, std::min(xi, xi_hi));
}

// Apply reforming extent to composition
static void apply_reforming(std::vector<double> &n, double xi, std::size_t i_hc,
                            int n_C, int n_H, const ReformingConfig &cfg) {
  n[i_hc] -= xi;
  n[cfg.i_H2O] -= n_C * xi;
  n[cfg.i_CO] += n_C * xi;
  n[cfg.i_H2] += (n_C + 0.5 * n_H) * xi;
}

// One round: reform each hydrocarbon in turn, then WGS, each to its own
// equilibrium at the current composition.
static void reforming_wgs_round(std::vector<double> &n_out, double T, double P,
                                const ReformingConfig &cfg) {
  for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
    std::size_t i_hc = cfg.hc_indices[k];
    int n_C = cfg.hc_C[k];
    int n_H = cfg.hc_H[k];

    double xi = solve_reforming_isothermal(n_out, T, P, i_hc, n_C, n_H, cfg);
    apply_reforming(n_out, xi, i_hc, n_C, n_H, cfg);
  }

  WgsConfig wgs_cfg = make_wgs_config();

  double xi_min = std::max(-n_out[wgs_cfg.i_CO2], -n_out[wgs_cfg.i_H2]);
  double xi_max = std::min(n_out[wgs_cfg.i_CO], n_out[wgs_cfg.i_H2O]);

  if (xi_max > xi_min + 1e-15) {
    WgsIsothermalFun F{n_out, T, wgs_cfg};
    double xi_wgs = newton_damped(F, 0.0, xi_min, xi_max, 40);

    n_out[wgs_cfg.i_CO] -= xi_wgs;
    n_out[wgs_cfg.i_H2O] -= xi_wgs;
    n_out[wgs_cfg.i_CO2] += xi_wgs;
    n_out[wgs_cfg.i_H2] += xi_wgs;
  }
}

// Largest change of any amount between two rounds, relative to the total.
static double reforming_round_change(const std::vector<double> &before,
                                     const std::vector<double> &after) {
  double total = 0.0, d = 0.0;
  for (std::size_t j = 0; j < after.size(); ++j) {
    total += after[j];
    d = std::max(d, std::abs(after[j] - before[j]));
  }
  return d / std::max(total, 1e-300);
}

// Solve general reforming + WGS equilibrium at fixed T: rounds until the
// composition stops moving (see REFORMING_PASS_RTOL). Returns whether it did.
static bool solve_reforming_wgs_isothermal(const std::vector<double> &n0,
                                           double T, double P,
                                           const ReformingConfig &cfg,
                                           std::vector<double> &n_out) {
  n_out = n0;
  std::vector<double> before;
  for (int pass = 0; pass < REFORMING_MAX_PASSES; ++pass) {
    before = n_out;
    reforming_wgs_round(n_out, T, P, cfg);
    if (reforming_round_change(before, n_out) <= REFORMING_PASS_RTOL) {
      return true;
    }
  }
  return false;
}

// Adiabatic solver functor for general reforming + WGS
struct ReformingAdiabaticFun {
  const std::vector<double> &n_in;
  double H_mass_target;
  double P;
  ReformingConfig cfg;

  mutable std::vector<double> n;
  mutable std::vector<double> X;

  ReformingAdiabaticFun(const std::vector<double> &n_in_, double H_mass_target_,
                        double P_, const ReformingConfig &cfg_)
      : n_in(n_in_), H_mass_target(H_mass_target_), P(P_), cfg(cfg_),
        n(n_in_.size()), X(n_in_.size()) {}

  double f(double T) const {
    solve_reforming_wgs_isothermal(n_in, T, P, cfg, n);
    molefractions(n, X);
    double H_mass = h(T, X) / mwmix(X);
    return H_mass - H_mass_target;
  }
};

// Public API: General reforming + WGS equilibrium (isothermal)
EquilibriumResult reforming_equilibrium(const State &in) {
  ReformingConfig cfg = make_reforming_config();

  // Check if any hydrocarbon is present
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (in.X[i_hc] > 1e-10) {
      has_hc = true;
      break;
    }
  }

  if (!has_hc) {
    // No hydrocarbons - fall back to WGS only
    return wgs_equilibrium(in);
  }

  std::vector<double> n_out;
  const bool conv = solve_reforming_wgs_isothermal(in.X, in.T, in.P, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(in.T, in.P, X_out);
  return EquilibriumResult{out, in.T, 0.0, conv};
}

// Public API: General reforming + WGS equilibrium (adiabatic)
EquilibriumResult reforming_equilibrium_adiabatic(const State &in) {
  ReformingConfig cfg = make_reforming_config();

  // Check if any hydrocarbon is present
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (in.X[i_hc] > 1e-10) {
      has_hc = true;
      break;
    }
  }

  if (!has_hc) {
    // No hydrocarbons - fall back to WGS only
    return wgs_equilibrium_adiabatic(in);
  }

  // Target enthalpy
  double H_mass_target = h(in.T, in.X) / mwmix(in.X);

  // Solve for adiabatic temperature
  ReformingAdiabaticFun F{in.X, H_mass_target, in.P, cfg};
  double T_ad = newton_damped(F, in.T, 300.0, 4500.0, 50);

  // Get equilibrium composition at T_ad
  std::vector<double> n_out;
  const bool conv = solve_reforming_wgs_isothermal(in.X, T_ad, in.P, cfg, n_out);

  std::vector<double> X_out;
  molefractions(n_out, X_out);

  State out;
  out.set_TPX(T_ad, in.P, X_out);
  return EquilibriumResult{out, in.T, T_ad - in.T, conv};
}

// -------------------------------------------------------------
// Combustion + Equilibrium (convenience function)
// -------------------------------------------------------------

EquilibriumResult combustion_equilibrium(const State &in, bool smooth_phi0, double k0) {
  // Step 1: Complete combustion (adiabatic)
  // We use Phi=0 smoothing (if requested) but NOT Phi=1 smoothing for equilibrium,
  // as equilibrium is naturally smooth at Phi=1.
  State burned = complete_combustion(in, smooth_phi0, false, k0, SMOOTHING_K_PHI1);

  // Step 2: Reforming + WGS equilibrium on combustion products
  // For rich mixtures: unburned HC + H2O -> CO + H2
  // For all mixtures: CO + H2O <-> CO2 + H2
  EquilibriumResult result = reforming_equilibrium_adiabatic(burned);
  return EquilibriumResult{result.state, in.T, result.state.T - in.T,
                           result.converged};
}

// -------------------------------------------------------------
// Combustion + equilibrium with analytic derivatives (#716)
// -------------------------------------------------------------

namespace {

// A forward-mode dual with a runtime number of partials: the solver's inputs
// (here the burned moles n0, T, P) plus one slot for the extent being
// differentiated implicitly.
struct DV {
  double v = 0.0;
  std::vector<double> d;
};
DV dv_const(double v, std::size_t N) { return DV{v, std::vector<double>(N, 0.0)}; }
DV dv_seed(double v, std::size_t N, std::size_t i) {
  DV r = dv_const(v, N);
  r.d[i] = 1.0;
  return r;
}
DV operator+(const DV &a, const DV &b) {
  DV r{a.v + b.v, a.d};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] += b.d[i];
  return r;
}
DV operator-(const DV &a, const DV &b) {
  DV r{a.v - b.v, a.d};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] -= b.d[i];
  return r;
}
DV operator*(double c, const DV &a) {
  DV r{c * a.v, a.d};
  for (double &x : r.d) x *= c;
  return r;
}
DV operator*(const DV &a, const DV &b) {
  DV r{a.v * b.v, std::vector<double>(a.d.size())};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
DV operator/(const DV &a, const DV &b) {
  DV r{a.v / b.v, std::vector<double>(a.d.size())};
  for (std::size_t i = 0; i < r.d.size(); ++i) r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
  return r;
}
DV dv_log(const DV &a) {
  DV r{std::log(a.v), a.d};
  for (double &x : r.d) x /= a.v;
  return r;
}
// f(a) with slope s, lifted.
DV dv_lift(double f, double s, const DV &a) {
  DV r{f, a.d};
  for (double &x : r.d) x *= s;
  return r;
}
DV dv_min(const DV &a, const DV &b) { return (a.v <= b.v) ? a : b; }

// An extent that the value path solved for: its implicit derivative where it
// is interior, the bound's where it sits on one. resid(xi) is the stage's
// residual with xi seeded in slot ix.
template <typename Resid>
DV implicit_extent(double xi, const DV &lo, const DV &hi, std::size_t ix, Resid &&resid) {
  const double span = std::max(hi.v - lo.v, 1e-300);
  if (xi >= hi.v - 1e-12 * span) {
    return hi;
  }
  if (xi <= lo.v + 1e-12 * span) {
    return lo;
  }
  const DV f = resid(dv_seed(xi, lo.d.size(), ix));
  DV out = dv_const(xi, lo.d.size());
  const double f_xi = f.d[ix];
  if (std::abs(f_xi) > 0.0) {
    for (std::size_t q = 0; q < out.d.size(); ++q) {
      if (q != ix) out.d[q] = -f.d[q] / f_xi;
    }
  }
  return out;
}

// d(g/RT)/dT = -h/(R T^2) = -h_RT / T.
DV g_over_RT_dv(std::size_t i, const DV &T) {
  return dv_lift(g_over_RT(i, T.v), -h_RT(i, T.v) / T.v, T);
}

}  // namespace

EquilibriumJacobian combustion_equilibrium_and_jacobian(const State &in, bool smooth_phi0,
                                                        double k0) {
  const std::size_t ns = in.X.size();
  const std::size_t iT = ns, iP = ns + 1, ix = ns + 2, N = ns + 3;
  EquilibriumJacobian out;
  const EquilibriumResult res = combustion_equilibrium(in, smooth_phi0, k0);
  out.state = res.state;
  const double T_ad = res.state.T;
  const double P = in.P;

  // Stage 1: complete combustion (phi = 1 unsmoothed, as the value path).
  const CompleteCombustionJacobian cj =
      complete_combustion_products_and_jacobian(in.X, smooth_phi0, false, k0, SMOOTHING_K_PHI1);
  const std::vector<double> &n0 = cj.X_b;

  // Stage 2: the same reforming + WGS rounds at T_ad, on duals.
  const DV T = dv_seed(T_ad, N, iT);
  const DV Pd = dv_seed(P, N, iP);
  const double P0 = 101325.0;
  std::vector<DV> n(ns);
  for (std::size_t j = 0; j < ns; ++j) n[j] = dv_seed(n0[j], N, j);
  auto values = [&]() {
    std::vector<double> v(ns);
    for (std::size_t j = 0; j < ns; ++j) v[j] = n[j].v;
    return v;
  };
  auto total = [&]() {
    DV t = dv_const(0.0, N);
    for (const DV &x : n) t = t + x;
    return t;
  };
  const ReformingConfig cfg = make_reforming_config();
  // As reforming_equilibrium_adiabatic: reform if the BURNED state has any
  // hydrocarbon, else WGS alone.
  bool has_hc = false;
  for (std::size_t i_hc : cfg.hc_indices) {
    if (n0[i_hc] > 1e-10) has_hc = true;
  }
  // Rounds until the composition stops moving, exactly as the value path
  // (solve_reforming_wgs_isothermal): the derivatives of the fixed point are
  // the limit of the rounds' chain rule.
  for (int pass = 0; pass < REFORMING_MAX_PASSES; ++pass) {
    const std::vector<double> before = values();
    if (has_hc) {
      for (std::size_t k = 0; k < cfg.hc_indices.size(); ++k) {
        const std::size_t i_hc = cfg.hc_indices[k];
        const int nC = cfg.hc_C[k];
        const int nH = cfg.hc_H[k];
        const std::vector<double> nv = values();
        const double xi = solve_reforming_isothermal(nv, T_ad, P, i_hc, nC, nH, cfg);
        DV xi_d = dv_const(0.0, N);
        if (xi != 0.0) {
          const DV hi = dv_min(n[i_hc], (1.0 / nC) * n[cfg.i_H2O]);
          const DV lo = dv_const(0.0, N) -
                        dv_min((1.0 / nC) * n[cfg.i_CO], (1.0 / (nC + 0.5 * nH)) * n[cfg.i_H2]);
          const double dn = nC + 0.5 * nH - 1.0;
          auto resid = [&](const DV &x) {
            const DV nt = total() + dn * x;
            const DV y_hc = (n[i_hc] - x) / nt;
            const DV y_h2o = (n[cfg.i_H2O] - static_cast<double>(nC) * x) / nt;
            const DV y_co = (n[cfg.i_CO] + static_cast<double>(nC) * x) / nt;
            const DV y_h2 = (n[cfg.i_H2] + (nC + 0.5 * nH) * x) / nt;
            const DV lnQ = static_cast<double>(nC) * dv_log(y_co) + (nC + 0.5 * nH) * dv_log(y_h2) -
                           dv_log(y_hc) - static_cast<double>(nC) * dv_log(y_h2o) +
                           dn * dv_log((1.0 / P0) * Pd);
            const DV dg = static_cast<double>(nC) * g_over_RT_dv(cfg.i_CO, T) +
                          (nC + 0.5 * nH) * g_over_RT_dv(cfg.i_H2, T) - g_over_RT_dv(i_hc, T) -
                          static_cast<double>(nC) * g_over_RT_dv(cfg.i_H2O, T);
            return lnQ + dg;  // lnQ - lnKp, lnKp = -dg
          };
          xi_d = implicit_extent(xi, lo, hi, ix, resid);
        }
        n[i_hc] = n[i_hc] - xi_d;
        n[cfg.i_H2O] = n[cfg.i_H2O] - static_cast<double>(nC) * xi_d;
        n[cfg.i_CO] = n[cfg.i_CO] + static_cast<double>(nC) * xi_d;
        n[cfg.i_H2] = n[cfg.i_H2] + (nC + 0.5 * nH) * xi_d;
      }
    }
    {
      const WgsConfig w = make_wgs_config();
      const std::vector<double> nv = values();
      const double lo_v = std::max(-nv[w.i_CO2], -nv[w.i_H2]);
      const double hi_v = std::min(nv[w.i_CO], nv[w.i_H2O]);
      if (hi_v > lo_v + 1e-15) {
        WgsIsothermalFun F{nv, T_ad, w};
        const double xi = newton_damped(F, 0.0, lo_v, hi_v, 40);
        const DV lo = (-nv[w.i_CO2] >= -nv[w.i_H2]) ? dv_const(0.0, N) - n[w.i_CO2]
                                                    : dv_const(0.0, N) - n[w.i_H2];
        const DV hi = dv_min(n[w.i_CO], n[w.i_H2O]);
        auto resid = [&](const DV &x) {
          const DV y_co2 = n[w.i_CO2] + x;
          const DV y_h2 = n[w.i_H2] + x;
          const DV y_co = n[w.i_CO] - x;
          const DV y_h2o = n[w.i_H2O] - x;
          const DV dg = g_over_RT_dv(w.i_CO2, T) + g_over_RT_dv(w.i_H2, T) -
                        g_over_RT_dv(w.i_CO, T) - g_over_RT_dv(w.i_H2O, T);
          // ln Q - ln Kp; WGS keeps the moles, so nt cancels in Q.
          return dv_log(y_co2) + dv_log(y_h2) - dv_log(y_co) - dv_log(y_h2o) + dg;
        };
        const DV xi_d = implicit_extent(xi, lo, hi, ix, resid);
        n[w.i_CO] = n[w.i_CO] - xi_d;
        n[w.i_H2O] = n[w.i_H2O] - xi_d;
        n[w.i_CO2] = n[w.i_CO2] + xi_d;
        n[w.i_H2] = n[w.i_H2] + xi_d;
      }
    }

    if (reforming_round_change(before, values()) <= REFORMING_PASS_RTOL) {
      break;
    }
  }
  // Stage 3: the adiabatic root F(T) = h(T, X_out)/MW_out - H_in = 0.
  const DV nt = total();
  DV h_mol = dv_const(0.0, N);
  DV mw = dv_const(0.0, N);
  for (std::size_t j = 0; j < ns; ++j) {
    const DV hj = dv_lift(h_species(j, T_ad), cp_species(j, T_ad), T);
    h_mol = h_mol + n[j] * hj;
    mw = mw + species_molar_mass(j) * n[j];
  }
  const DV H_mass = h_mol / mw;  // n/nt cancels in h/MW
  (void)nt;
  const double F_T = H_mass.d[iT];
  const double mw_in = mwmix(in.X);
  const double h_in = h(in.T, in.X);
  const double H_in = h_in / mw_in;

  out.dT_dT_in = cp(in.T, in.X) / mw_in / F_T;
  out.dT_dP = -H_mass.d[iP] / F_T;
  out.dT_dX.assign(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    double dF = 0.0;
    for (std::size_t j = 0; j < ns; ++j) dF += H_mass.d[j] * cj.dXb_dX[j][k];
    const double dH = (h_species(k, in.T) - H_in * species_molar_mass(k)) / mw_in;
    out.dT_dX[k] = (dH - dF) / F_T;
  }

  // Product mole fractions: n_out(n0, T, P), T = T_ad(X, T_in, P).
  std::vector<double> nvals(ns);
  double ntv = 0.0;
  for (std::size_t j = 0; j < ns; ++j) {
    nvals[j] = n[j].v;
    ntv += n[j].v;
  }
  auto dX_from_dn = [&](const std::vector<double> &dn) {
    double dN = 0.0;
    for (double x : dn) dN += x;
    std::vector<double> dX(ns);
    for (std::size_t j = 0; j < ns; ++j) dX[j] = (dn[j] - nvals[j] / ntv * dN) / ntv;
    return dX;
  };
  out.dXout_dX.assign(ns, std::vector<double>(ns, 0.0));
  for (std::size_t k = 0; k < ns; ++k) {
    std::vector<double> dn(ns, 0.0);
    for (std::size_t j = 0; j < ns; ++j) {
      double v = n[j].d[iT] * out.dT_dX[k];
      for (std::size_t l = 0; l < ns; ++l) v += n[j].d[l] * cj.dXb_dX[l][k];
      dn[j] = v;
    }
    const auto dX = dX_from_dn(dn);
    for (std::size_t j = 0; j < ns; ++j) out.dXout_dX[j][k] = dX[j];
  }
  std::vector<double> dn_T(ns), dn_P(ns);
  for (std::size_t j = 0; j < ns; ++j) {
    dn_T[j] = n[j].d[iT] * out.dT_dT_in;
    dn_P[j] = n[j].d[iT] * out.dT_dP + n[j].d[iP];
  }
  out.dXout_dT_in = dX_from_dn(dn_T);
  out.dXout_dP = dX_from_dn(dn_P);
  return out;
}

} // namespace combaero
