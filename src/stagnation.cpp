// Isentropic stagnation/static conversion utilities.
// Uses the variable-cp thermo backend (thermo.h) throughout.

#include "../include/stagnation.h"
#include "../include/composition.h"
#include "../include/compressible.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace combaero {

// -------------------------------------------------------------
// Kinetic energy / stagnation enthalpy
// -------------------------------------------------------------

double h0_from_static(double h_static_J_per_kg, double v) {
    return h_static_J_per_kg + 0.5 * v * v;
}

double v_from_h0(double h0, double h_static_J_per_kg) {
    double dh = h0 - h_static_J_per_kg;
    return (dh > 0.0) ? std::sqrt(2.0 * dh) : 0.0;
}

// -------------------------------------------------------------
// Mach number
// -------------------------------------------------------------

double mach_number(double v, double T, const std::vector<double>& X) {
    double a = speed_of_sound(T, X);
    if (a <= 0.0) {
        throw std::invalid_argument("mach_number: speed of sound is non-positive");
    }
    return v / a;
}

// -------------------------------------------------------------
// Stagnation temperature
// -------------------------------------------------------------

// Helper: h_mass [J/kg] at temperature T for composition X
static double h_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return h(T, X) / mw_g * 1000.0;  // J/mol / (g/mol) * 1000 = J/kg
}

// Helper: cp_mass [J/(kg*K)] at temperature T for composition X
static double cp_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return cp(T, X) / mw_g * 1000.0;  // J/(mol*K) / (g/mol) * 1000 = J/(kg*K)
}

double T0_from_static_v(double T, double v, const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (T <= 0.0) {
        throw std::invalid_argument("T0_from_static_v: T must be positive");
    }
    if (v < 0.0) {
        throw std::invalid_argument("T0_from_static_v: v must be non-negative");
    }

    double h_static = h_mass_at(T, X);
    double h0 = h_static + 0.5 * v * v;

    // Newton iteration: find T0 such that h(T0) = h0
    // F(T0) = h(T0) - h0 = 0,  dF/dT0 = cp(T0)
    double T0 = T + 0.5 * v * v / cp_mass_at(T, X);  // ideal-gas initial guess

    constexpr double T_MIN = thermo::T_FLOOR;  // as compressible.cpp solve_T_isentropic

    for (std::size_t it = 0; it < max_iter; ++it) {
        double F = h_mass_at(T0, X) - h0;
        double dF = cp_mass_at(T0, X);

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T0;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T0 += dT;
        T0 = std::max(T0, T_MIN);

        if (std::abs(dT) < tol * (1.0 + std::abs(T0))) {
            return T0;
        }
    }
    return T0;
}

double T0_from_static(double T, double M, const std::vector<double>& X,
                       double tol, std::size_t max_iter) {
    if (M < 0.0) {
        throw std::invalid_argument("T0_from_static: M must be non-negative");
    }
    double a = speed_of_sound(T, X);
    double v = M * a;
    return T0_from_static_v(T, v, X, tol, max_iter);
}

double T_from_stagnation(double T0, double M, const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    if (T0 <= 0.0) {
        throw std::invalid_argument("T_from_stagnation: T0 must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("T_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return T0;

    double h0 = h_mass_at(T0, X);

    // Newton iteration: find T such that h(T) + M^2*a(T)^2/2 = h0
    // F(T) = h(T) + M^2*a(T)^2/2 - h0
    // dF/dT = cp(T) + M^2 * a * da/dT
    // da/dT = (1/(2a)) * d(gamma*R*T)/dT  -- approximate as a/(2T) for ideal gas
    double T = T0 / (1.0 + 0.2 * M * M);  // calorically perfect gas initial guess

    constexpr double T_MIN = thermo::T_FLOOR;  // as compressible.cpp solve_T_isentropic

    for (std::size_t it = 0; it < max_iter; ++it) {
        double a_val = speed_of_sound(T, X);
        double v = M * a_val;
        double F = h_mass_at(T, X) + 0.5 * v * v - h0;
        // dF/dT ~ cp + M^2 * a^2 / (2T)
        double dF = cp_mass_at(T, X) + M * M * a_val * a_val / (2.0 * T);

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T += dT;
        T = std::max(T, T_MIN);

        if (std::abs(dT) < tol * (1.0 + std::abs(T))) {
            return T;
        }
    }
    return T;
}

// -------------------------------------------------------------
// Stagnation pressure (isentropic)
// -------------------------------------------------------------

double P0_from_static(double P, double T, double M,
                       const std::vector<double>& X,
                       double tol, std::size_t max_iter) {
    if (P <= 0.0) {
        throw std::invalid_argument("P0_from_static: P must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("P0_from_static: M must be non-negative");
    }
    if (M == 0.0) return P;

    constexpr double P_REF = 101325.0;
    double s0 = s(T, X, P, P_REF);
    double T0 = T0_from_static(T, M, X, tol, max_iter);

    double R_specific = specific_gas_constant(X);  // J/(kg*K)
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    double s_T0_Pref = s(T0, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass   = s0 / mw_kg;                        // J/(kg*K)
    double ln_ratio = (s_T0_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

double P_from_stagnation(double P0, double T0, double M,
                          const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    if (P0 <= 0.0) {
        throw std::invalid_argument("P_from_stagnation: P0 must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("P_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return P0;

    constexpr double P_REF = 101325.0;
    double s0 = s(T0, X, P0, P_REF);
    double T_static = T_from_stagnation(T0, M, X, tol, max_iter);

    double R_specific = specific_gas_constant(X);  // J/(kg*K)
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    double s_T_Pref = s(T_static, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass  = s0 / mw_kg;                              // J/(kg*K)
    double ln_ratio = (s_T_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

// -------------------------------------------------------------
// Adiabatic wall temperature
// -------------------------------------------------------------

double recovery_factor(double Pr, bool turbulent) {
    if (Pr <= 0.0) {
        throw std::invalid_argument("recovery_factor: Pr must be positive");
    }
    return turbulent ? std::cbrt(Pr) : std::sqrt(Pr);
}

double T_adiabatic_wall(double T_static, double v,
                         double T, double P, const std::vector<double>& X,
                         bool turbulent) {
    double T_safe = std::max(T_static, 10.0);
    double v_abs = std::abs(v);
    if (v_abs < 1e-12) return T_safe;


    double Pr = prandtl(T, P, X);
    double r = recovery_factor(Pr, turbulent);
    double cp_val = cp_mass_at(T, X);

    return T_static + r * 0.5 * v * v / cp_val;
}

double T_adiabatic_wall_mach(double T_static, double M,
                               double T, double P, const std::vector<double>& X,
                               bool turbulent) {
    if (M < 0.0) {
        throw std::invalid_argument("T_adiabatic_wall_mach: M must be non-negative");
    }
    if (M == 0.0) return T_static;

    double a = speed_of_sound(T, X);
    double v = M * a;
    return T_adiabatic_wall(T_static, v, T, P, X, turbulent);
}

std::tuple<double, double> stagnation_from_static(
    double T, double P, double v,
    const std::vector<double>& X,
    double tol, std::size_t max_iter)
{
    if (T <= 0.0)
        throw std::invalid_argument("stagnation_from_static: T must be positive");
    if (P <= 0.0)
        throw std::invalid_argument("stagnation_from_static: P must be positive");
    if (v < 0.0)
        throw std::invalid_argument("stagnation_from_static: v must be non-negative");
    if (v == 0.0) return {T, P};

    // One Newton loop for Tt.
    double Tt = T0_from_static_v(T, v, X, tol, max_iter);

    // Pt from isentropic entropy conservation: s(T, P) = s(Tt, Pt).
    // Bypasses P0_from_static to avoid its internal T0_from_static call,
    // which would repeat the Newton loop and re-evaluate speed_of_sound(T, X).
    constexpr double P_REF = 101325.0;
    const double mw_kg  = mwmix(X) / 1000.0;
    const double R_spec = specific_gas_constant(X);
    const double s_static = s(T,  X, P,     P_REF) / mw_kg;
    const double s_stag   = s(Tt, X, P_REF, P_REF) / mw_kg;
    double Pt = P_REF * std::exp((s_stag - s_static) / R_spec);

    return {Tt, Pt};
}

namespace {
// Integral of f over [a, b] by 8-point Gauss-Legendre, split at the
// NASA-9 range boundaries (1000 K, 6000 K) where cp has a small jump (#549).
// Used for h(T0) - h(T) and s(T0) - s(T) at the SAME composition: the direct
// differences cancel ~4 significant digits near T0 and left ~1e-9 Pa of noise
// in P0, which finite-difference Jacobian checks with small steps read as
// derivative errors (#579 on Linux).
template <class F>
double gl8_split(F f, double a, double b) {
  static const double xg[4] = {0.1834346424956498, 0.5255324099163290,
                               0.7966664774136267, 0.9602898564975363};
  static const double wg[4] = {0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  auto seg = [&](double lo, double hi) {
    const double c = 0.5 * (lo + hi), r = 0.5 * (hi - lo);
    double acc = 0.0;
    for (int i = 0; i < 4; ++i) {
      acc += wg[i] * (f(c - r * xg[i]) + f(c + r * xg[i]));
    }
    return acc * r;
  };
  const double sign = b >= a ? 1.0 : -1.0;
  double lo = std::min(a, b), hi = std::max(a, b);
  double total = 0.0;
  for (double bp : {1000.0, 6000.0}) {
    if (lo < bp && bp < hi) {
      total += seg(lo, bp);
      lo = bp;
    }
  }
  total += seg(lo, hi);
  return sign * total;
}

}  // namespace

ChamberStaticState chamber_static_state(double T0, double P, double G,
                                        const std::vector<double> &X) {
  const double R = specific_gas_constant(X);
  const double k = (G * R / P) * (G * R / P);  // v^2 = k T^2
  auto cp_of = [&](double t) { return cp_mass(t, X); };
  auto cp_over_T = [&](double t) { return cp_mass(t, X) / t; };
  // f(T) = [h(T0) - h(T)] - k T^2 / 2, the enthalpy difference integrated.
  double T = T0;
  for (int it = 0; it < 100; ++it) {
    const double f = gl8_split(cp_of, T, T0) - 0.5 * k * T * T;
    const double fT = -cp_mass(T, X) - k * T;
    const double step = f / fT;
    T -= step;
    if (std::abs(step) <= 1e-13 * T0) {
      break;
    }
  }
  ChamberStaticState out{};
  out.T = T;
  const double cp_T = cp_mass(T, X);
  const double cp_T0 = cp_mass(T0, X);
  const double fT = -cp_T - k * T;
  // df/dx for x = T0, P, G; dT/dx = -(df/dx) / fT.
  out.dT_dT0 = -cp_T0 / fT;
  out.dT_dP = -(k * T * T / P) / fT;
  out.dT_dG = (G > 0.0) ? -(-k * T * T / G) / fT : 0.0;
  // ln(P0 / P) = (s(T0) - s(T)) / R = (1/R) integral of cp / T.
  out.lnP0P = gl8_split(cp_over_T, T, T0) / R;
  const double dln_dTs = -cp_T / (R * T);
  out.dln_dT0 = cp_T0 / (R * T0) + dln_dTs * out.dT_dT0;
  out.dln_dP = dln_dTs * out.dT_dP;
  out.dln_dG = dln_dTs * out.dT_dG;
  const double a = speed_of_sound(T, X);
  out.M = (a > 1e-9) ? G * R * T / (P * a) : 0.0;
  return out;
}

ChamberStaticStateDX chamber_static_state_dX(double T0, double P, double G,
                                            const std::vector<double> &X) {
  const std::size_t n = X.size();
  ChamberStaticStateDX out;
  out.dT_dX.assign(n, 0.0);
  out.dlnP0P_dX.assign(n, 0.0);
  const ChamberStaticState cs = chamber_static_state(T0, P, G, X);
  const double T = cs.T;
  const double R = specific_gas_constant(X);
  const double MW = mwmix(X) / 1000.0;  // kg/mol
  const double Ru = R * MW;
  const double cpT = cp_mass(T, X);
  const double v = G * R * T / P;
  const double fT = -cpT - (G * R / P) * (G * R / P) * T;
  double dH = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    dH += X[j] * (h_species(j, T0) - h_species(j, T));
  }
  dH /= MW;
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    // h_mass difference and R = R_u / MW both move with X_j.
    const double fX = (h_species(j, T0) - h_species(j, T) - dH * W) / MW + v * v * W / MW;
    out.dT_dX[j] = -fX / fT;
    out.dlnP0P_dX[j] =
        (s_species(j, T0) - s_species(j, T)) / Ru - cpT / (R * T) * out.dT_dX[j];
  }
  return out;
}

SonicStaticState sonic_static_state(double T0, double G,
                                    const std::vector<double> &X) {
  const double R = specific_gas_constant(X);
  const double mw_g = mwmix(X);  // g/mol
  // cp, gamma and their T derivatives on the mass basis.
  auto cp_m = [&](double t) { return cp_mass(t, X); };
  auto gamma_and_d = [&](double t, double &dg) {
    const double c = cp_mass(t, X);
    const double dc = J_per_mol_to_J_per_kg(dcp_dT(t, X), mw_g);
    const double cv = c - R;
    dg = -R * dc / (cv * cv);
    return c / cv;
  };
  // F(T) = [h(T0) - h(T)] - a(T)^2 / 2, a^2 = gamma R T; F falls in T.
  double dg0 = 0.0;
  const double g0 = gamma_and_d(T0, dg0);
  double T = 2.0 * T0 / (g0 + 1.0);
  double dFdT = -1.0;
  for (int it = 0; it < 100; ++it) {
    double dg = 0.0;
    const double g = gamma_and_d(T, dg);
    const double F = gl8_split(cp_m, T, T0) - 0.5 * g * R * T;
    dFdT = -cp_mass(T, X) - 0.5 * R * (g + T * dg);
    const double step = F / dFdT;
    T -= step;
    if (std::abs(step) <= 1e-13 * T0) {
      break;
    }
  }
  double dg = 0.0;
  const double g = gamma_and_d(T, dg);
  dFdT = -cp_mass(T, X) - 0.5 * R * (g + T * dg);
  SonicStaticState out{};
  out.T = T;
  out.dT_dT0 = cp_mass(T0, X) / -dFdT;
  const double Gabs = std::abs(G);
  out.P = Gabs * std::sqrt(R * T / g);
  // d ln P* / d T* = 1/(2 T*) - gamma'/(2 gamma).
  const double dlnP_dT = 0.5 / T - 0.5 * dg / g;
  out.dP_dT0 = out.P * dlnP_dT * out.dT_dT0;
  out.dP_dG = (Gabs > 0.0) ? out.P / Gabs * (G >= 0.0 ? 1.0 : -1.0) : std::sqrt(R * T / g);
  return out;
}

DuctStaticState duct_static_state(double T0, double P0, double G,
                                  const std::vector<double> &X) {
  // An unphysical iterate is evaluated at the floors, finite; the partials
  // in a floored argument are then zero.
  const bool T0_floored = !(T0 > thermo::T_FLOOR);
  const bool P0_floored = !(P0 > thermo::P_FLOOR);
  T0 = T0_floored ? thermo::T_FLOOR : T0;
  P0 = P0_floored ? thermo::P_FLOOR : P0;
  const double R = specific_gas_constant(X);
  const double lnP0 = std::log(P0);
  G = std::abs(G);
  // L(u) = lnP0P(T0, e^u, G) + u - ln P0, increasing in u on [ln P*, ln P0].
  auto eval = [&](double u, ChamberStaticState &cs) {
    cs = chamber_static_state(T0, std::exp(u), G, X);
    return cs.lnP0P + u - lnP0;
  };
  DuctStaticState out{};
  ChamberStaticState cs{};
  const double P_star = critical_pressure_ratio(T0, P0, X) * P0;
  double u_lo = std::log(P_star);
  double u_hi = lnP0;
  if (eval(u_lo, cs) > 0.0) {
    // Beyond choke: the sonic state, at the choked flux.
    const double G_star = mass_flux_isentropic(T0, P0, P_star, X);
    cs = chamber_static_state(T0, P_star, G_star, X);
    out.T = cs.T;
    out.P = P_star;
    out.rho = P_star / (R * cs.T);
    out.M = cs.M;
    out.choked = true;
    out.G_eff = G;
    out.dGeff_dG = 1.0;
    const double r = P_star / P0;
    out.dP_dP0 = r;
    out.dT_dT0 = cs.dT_dT0;
    out.dT_dP0 = cs.dT_dP * r;
    out.drho_dT0 = -out.rho * out.dT_dT0 / out.T;
    out.drho_dP0 = out.rho * (out.dP_dP0 / out.P - out.dT_dP0 / out.T);
    if (T0_floored) {
      out.dT_dT0 = out.dP_dT0 = out.drho_dT0 = 0.0;
    }
    if (P0_floored) {
      out.dT_dP0 = out.dP_dP0 = out.drho_dP0 = 0.0;
    }
    return out;
  }
  double u = u_hi;
  double L = eval(u, cs);
  for (int it = 0; it < 100 && L != 0.0; ++it) {
    if (L > 0.0) {
      u_hi = u;
    } else {
      u_lo = u;
    }
    const double dL = 1.0 + std::exp(u) * cs.dln_dP;
    double u_new = (dL > 0.0) ? u - L / dL : 0.5 * (u_lo + u_hi);
    if (!(u_new > u_lo && u_new < u_hi)) {
      u_new = 0.5 * (u_lo + u_hi);
    }
    const bool done = std::abs(u_new - u) <= 1e-15 * std::abs(lnP0) + 1e-15;
    u = u_new;
    L = eval(u, cs);
    if (done) {
      break;
    }
  }
  const double P = std::exp(u);
  const double dL_du = 1.0 + P * cs.dln_dP;
  // du/dx = -(dL/dx) / (dL/du); dL/dT0 = dln_dT0, dL/dP0 = -1/P0, dL/dG = dln_dG.
  const double du_dT0 = -cs.dln_dT0 / dL_du;
  const double du_dP0 = (1.0 / P0) / dL_du;
  const double du_dG = -cs.dln_dG / dL_du;
  out.P = P;
  out.T = cs.T;
  out.rho = P / (R * cs.T);
  out.M = cs.M;
  out.choked = false;
  out.dP_dT0 = P * du_dT0;
  out.dP_dP0 = P * du_dP0;
  out.dP_dG = P * du_dG;
  out.dT_dT0 = cs.dT_dT0 + cs.dT_dP * out.dP_dT0;
  out.dT_dP0 = cs.dT_dP * out.dP_dP0;
  out.dT_dG = cs.dT_dG + cs.dT_dP * out.dP_dG;
  out.drho_dT0 = out.rho * (out.dP_dT0 / P - out.dT_dT0 / out.T);
  out.drho_dP0 = out.rho * (out.dP_dP0 / P - out.dT_dP0 / out.T);
  out.drho_dG = out.rho * (out.dP_dG / P - out.dT_dG / out.T);
  out.G_eff = G;
  out.dGeff_dG = 1.0;
  if (T0_floored) {
    out.dT_dT0 = out.dP_dT0 = out.drho_dT0 = 0.0;
  }
  if (P0_floored) {
    out.dT_dP0 = out.dP_dP0 = out.drho_dP0 = 0.0;
  }
  return out;
}

DuctStaticState duct_htc_static_state(double T0, double P0, double G,
                                      const std::vector<double> &X) {
  G = std::abs(G);
  if (!(T0 > thermo::T_FLOOR) || !(P0 > thermo::P_FLOOR)) {
    return duct_static_state(T0, P0, G, X);
  }
  auto sonic_flux = [&](double T, double P) {
    return mass_flux_isentropic(T, P, critical_pressure_ratio(T, P, X) * P, X);
  };
  const double Gs = sonic_flux(T0, P0);
  const double w = kDuctHtcSonicMargin;
  const double a = 1.0 - 2.0 * w;
  const double g = G / Gs;
  if (g <= a) {
    return duct_static_state(T0, P0, G, X);
  }
  // g_eff = g - (g - a)^2 / (4 w) on [a, 1]: value and slope meet g at a and
  // the cap 1 - w (slope 0) at 1.
  double ge = 1.0 - w;
  double dge = 0.0;
  if (g < 1.0) {
    ge = g - (g - a) * (g - a) / (4.0 * w);
    dge = 1.0 - (g - a) / (2.0 * w);
  }
  DuctStaticState s = duct_static_state(T0, P0, ge * Gs, X);
  // G_eff = Gs(T0, P0) ge(G / Gs): d/dG = dge, d/dGs = ge - g dge. The sonic
  // flux is proportional to P0 at fixed T0 (ideal gas: P*/P0 depends on T0
  // only); its T0 slope is differenced (smooth, step 1e-6).
  const double hT = 1e-6 * T0;
  const double dGs_dT0 = (sonic_flux(T0 + hT, P0) - sonic_flux(T0 - hT, P0)) / (2.0 * hT);
  const double dGs_dP0 = Gs / P0;
  const double dGe_dGs = ge - g * dge;
  s.G_eff = ge * Gs;
  s.dGeff_dG = dge;
  s.dGeff_dT0 = dGe_dGs * dGs_dT0;
  s.dGeff_dP0 = dGe_dGs * dGs_dP0;
  // Chain the state's G partials through G_eff.
  auto chain = [&](double &d_T0, double &d_P0, double &d_G) {
    d_T0 += d_G * s.dGeff_dT0;
    d_P0 += d_G * s.dGeff_dP0;
    d_G *= s.dGeff_dG;
  };
  chain(s.dT_dT0, s.dT_dP0, s.dT_dG);
  chain(s.dP_dT0, s.dP_dP0, s.dP_dG);
  chain(s.drho_dT0, s.drho_dP0, s.drho_dG);
  return s;
}

DuctStaticStateDY duct_htc_static_state_dY(double T0, double P0, double G,
                                           const std::vector<double> &X,
                                           const DuctStaticState &state) {
  const std::size_t n = X.size();
  DuctStaticStateDY out;
  out.dT_dY.assign(n, 0.0);
  out.dP_dY.assign(n, 0.0);
  out.drho_dY.assign(n, 0.0);
  out.dGeff_dY.assign(n, 0.0);
  G = std::abs(G);
  if (!(T0 > thermo::T_FLOOR) || !(P0 > thermo::P_FLOOR)) {
    return out;
  }
  // The caller's solved state. Below the bend it IS the plain state at G
  // (G_eff = G exactly) and nothing is solved again; in the bend, the raw
  // state at G_eff gives its own G partials, before the bend's chain.
  const DuctStaticState &s = state;
  const bool bent = s.G_eff != G;
  const DuctStaticState raw = bent ? duct_static_state(T0, P0, s.G_eff, X) : s;
  auto sonic_flux = [&](const std::vector<double> &Xv) {
    return mass_flux_isentropic(T0, P0, critical_pressure_ratio(T0, P0, Xv) * P0, Xv);
  };
  const double w = kDuctHtcSonicMargin;
  const double a = 1.0 - 2.0 * w;
  const double T = s.T;
  const double P = s.P;
  const double R = specific_gas_constant(X);
  const double MW = mwmix(X) / 1000.0;  // kg/mol
  const double Ru = R * MW;
  const double cpT = cp_mass(T, X);
  const double v = s.G_eff * R * T / P;
  const double kk = (s.G_eff * R / P) * (s.G_eff * R / P);
  // Mass-specific enthalpy drop, molar species terms over MW.
  // Species enthalpy and entropy drops T0 -> T, once.
  std::vector<double> dh_j(n), ds_j(n);
  double dH = 0.0;
  for (std::size_t j = 0; j < n; ++j) {
    dh_j[j] = h_species(j, T0) - h_species(j, T);
    ds_j[j] = s_species(j, T0) - s_species(j, T);
    dH += X[j] * dh_j[j];
  }
  dH /= MW;
  // Jacobian in (T, u = ln P): f_T, f_u; L_T, L_u = 1.
  const double fT = -cpT - kk * T;
  const double fu = kk * T * T;
  const double LT = -cpT / (R * T);
  const double det = fT - fu * LT;
  std::vector<double> dT_dX(n, 0.0), du_dX(n, 0.0), drho_dX(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    // dR/dX_j = -R W/MW enters f through -(v^2 / R) dR.
    const double fX = (dh_j[j] - dH * W) / MW + v * v * W / MW;
    const double LX = ds_j[j] / Ru;
    dT_dX[j] = -(fX - fu * LX) / det;
    du_dX[j] = -(fT * LX - LT * fX) / det;
  }
  // The bent flux near choke: G_eff = ge G* with G*(T0, P0, X).
  std::vector<double> dGe_dX(n, 0.0);
  if (bent) {
    const double Gs = sonic_flux(X);
    const double g = G / Gs;
    double ge = 1.0 - w;
    double dge = 0.0;
    if (g < 1.0) {
      ge = g - (g - a) * (g - a) / (4.0 * w);
      dge = 1.0 - (g - a) / (2.0 * w);
    }
    const double dGe_dGs = ge - g * dge;
    for (std::size_t j = 0; j < n; ++j) {
      // One-sided off an absent species, as it enters.
      const double hX = 1e-6;
      std::vector<double> Xp = X, Xm = X;
      Xp[j] += hX;
      double dGs = 0.0;
      if (X[j] > hX) {
        Xm[j] -= hX;
        dGs = (sonic_flux(Xp) - sonic_flux(Xm)) / (2.0 * hX);
      } else {
        dGs = (sonic_flux(Xp) - Gs) / hX;
      }
      dGe_dX[j] = dGe_dGs * dGs;
      dT_dX[j] += raw.dT_dG * dGe_dX[j];
      du_dX[j] += raw.dP_dG / P * dGe_dX[j];
    }
  }
  std::vector<double> dP_dX(n, 0.0);
  for (std::size_t j = 0; j < n; ++j) {
    const double W = species_molar_mass(j) / 1000.0;
    dP_dX[j] = P * du_dX[j];
    // rho = P MW / (R_u T)
    drho_dX[j] = s.rho * (du_dX[j] - dT_dX[j] / T + W / MW);
  }
  out.dT_dY = mole_gradient_to_mass(X, dT_dX);
  out.dP_dY = mole_gradient_to_mass(X, dP_dX);
  out.drho_dY = mole_gradient_to_mass(X, drho_dX);
  out.dGeff_dY = mole_gradient_to_mass(X, dGe_dX);
  return out;
}

} // namespace combaero
