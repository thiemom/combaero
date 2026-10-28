// 1-D supersonic ejector model (critical / subcritical / unchoked jet-pump
// regimes). See include/ejector.h for paper references and the scope note.

#include "ejector.h"

#include "dual_number.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <utility>

namespace combaero::solver {

double ejector_area_over_astar(double mach, double gamma) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  return (1.0 / mach) *
         std::pow((2.0 / gp1) * (1.0 + 0.5 * gm1 * mach * mach), gp1 / (2.0 * gm1));
}

double ejector_mach_from_area_ratio_supersonic(double area_ratio, double gamma) {
  double lo = 1.0 + 1e-12;
  double hi = 50.0;
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (lo + hi);
    if (ejector_area_over_astar(mid, gamma) < area_ratio) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

double ejector_choked_mass_flow(double Pt, double Tt, double area_throat,
                                double gamma, double r_gas, double eta) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double choke_const = std::sqrt(gamma / r_gas * std::pow(2.0 / gp1, gp1 / gm1));
  return (Pt * area_throat / std::sqrt(Tt)) * choke_const * std::sqrt(eta);
}

std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian(double Pt, double Tt, double area_throat,
                                      double gamma, double r_gas, double eta) {
  double mdot = ejector_choked_mass_flow(Pt, Tt, area_throat, gamma, r_gas, eta);
  // mdot = C * Pt / sqrt(Tt)  =>  d/dPt = mdot/Pt,  d/dTt = -mdot/(2*Tt)
  return {mdot, mdot / Pt, -0.5 * mdot / Tt};
}

namespace {

// Forward C-D nozzle flow from stagnation p_up to exit static p_down <= p_up:
// the subsonic exit flux smooth-min'd with the sonic-throat cap.
double cd_nozzle_forward(double p_up, double Tt, double p_down, double area_throat,
                         double area_exit, double gamma, double r_gas, double eta,
                         double eps_frac) {
  double cap = ejector_choked_mass_flow(p_up, Tt, area_throat, gamma, r_gas, eta);
  double exit_flux = 0.0;
  double r = p_down / p_up;
  if (r < 1.0) {
    double gm1 = gamma - 1.0;
    double gp1 = gamma + 1.0;
    double mach = std::sqrt(2.0 / gm1 * (std::pow(r, -gm1 / gamma) - 1.0));
    // Below the critical ratio the throat is choked: hold the exit flux at
    // its sonic maximum (flux'(M = 1) = 0, so C1). The supersonic-exit flux
    // FALLS with M and dipped under the cap at deep expansion (r < 0.16 at
    // A_e/A_t = 1.5), a non-monotone nozzle flow in P_py (#619).
    mach = std::min(mach, 1.0);
    double flux = mach * std::pow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
    exit_flux = area_exit * p_up / std::sqrt(Tt) * std::sqrt(gamma / r_gas) * flux * std::sqrt(eta);
  }
  double eps = eps_frac * cap;
  double diff = exit_flux - cap;
  return 0.5 * (exit_flux + cap - std::sqrt(diff * diff + eps * eps));
}

} // namespace

double ejector_cd_nozzle_mass_flow(double Pt, double Tt, double p_static_down,
                                   double area_throat, double area_exit,
                                   double gamma, double r_gas, double eta,
                                   double eps_frac) {
  if (p_static_down <= Pt) {
    return cd_nozzle_forward(Pt, Tt, p_static_down, area_throat, area_exit, gamma, r_gas, eta,
                             eps_frac);
  }
  // Exit static above the supply (#619): the nozzle passes REVERSE flow, the
  // same closure from p_static_down back into Pt, shifted by the smooth-min's
  // r = 1 value so the flow is continuous there. The former zero exit flux
  // was flat in P_py, which let zero flow at any P_py above both supplies
  // close the ejector rows (a spurious root at P_py ~ 1e13 Pa).
  double cap = ejector_choked_mass_flow(Pt, Tt, area_throat, gamma, r_gas, eta);
  double eps = eps_frac * cap;
  double offset = std::sqrt(cap * cap + eps * eps) - cap;
  return -cd_nozzle_forward(p_static_down, Tt, Pt, area_throat, area_exit, gamma, r_gas, eta,
                            eps_frac) -
         offset;
}

double ejector_q_lambda(double lam, double gamma) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  return std::pow(1.0 - (gm1 / gp1) * lam * lam, 1.0 / gm1) * std::pow(gp1 / 2.0, 1.0 / gm1) * lam;
}

// The generic forward-mode dual used below (DualN<N>) now lives in
// include/dual_number.h, so the junction closure can share one implementation
// rather than a second copy (issue #271). The fixed 4-seed Dual4 above is
// deliberately left here: it backs the validated critical-mode path.

EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian(
    double Pt_val, double Tt_val, double p_py_val, double area_throat,
    double area_exit, double gamma, double r_gas, double eta, double eps_frac) {
  // Seeds: 0 = p0 (Pt), 1 = t0 (Tt), 2 = p_py.
  using D = DualN<3>;
  D p0 = D::seed(Pt_val, 0);
  D t0 = D::seed(Tt_val, 1);
  D p_py = D::seed(p_py_val, 2);

  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double choke_const = std::sqrt(gamma / r_gas * std::pow(2.0 / gp1, gp1 / gm1));
  auto cap_of = [&](const D& p_up) {
    return (choke_const * area_throat * std::sqrt(eta)) * p_up / dsqrt(t0);
  };
  // eps = eps_frac * cap depends on (p_up, t0); differentiate it too (FD of
  // the Python value function does), else the Jacobian would be inconsistent.
  auto forward = [&](const D& p_up, const D& p_down) {
    D cap = cap_of(p_up);
    D exit_flux = D::constant(0.0);
    if (p_down.v < p_up.v) {
      D r = p_down / p_up;
      D mach_sq = (2.0 / gm1) * (dpow(r, -gm1 / gamma) - 1.0);
      D mach = dsqrt(mach_sq);
      if (mach.v > 1.0) mach = D::constant(1.0); // choked: see ejector_cd_nozzle_mass_flow
      D flux = mach * dpow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
      exit_flux = (area_exit * std::sqrt(gamma / r_gas) * std::sqrt(eta)) * p_up / dsqrt(t0) * flux;
    }
    D eps = eps_frac * cap;
    D diff = exit_flux - cap;
    return 0.5 * (exit_flux + cap - dsqrt(diff * diff + eps * eps));
  };

  D mdot = D::constant(0.0);
  if (p_py_val <= Pt_val) {
    mdot = forward(p0, p_py);
  } else {
    // Reverse flow above the supply (#619), continuous at P_py = p0: see
    // ejector_cd_nozzle_mass_flow.
    D cap = cap_of(p0);
    D eps = eps_frac * cap;
    mdot = -forward(p_py, p0) - (dsqrt(cap * cap + eps * eps) - cap);
  }
  return {mdot.v, mdot.d[0], mdot.d[1], mdot.d[2]};
}

EjectorCDExitStaticJacobian ejector_cd_nozzle_exit_static_and_jacobian(
    double m_dot, double Pt, double Tt, double area_throat, double area_exit,
    double gamma, double r_gas, double eta, double eps_frac) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double r_crit = std::pow(2.0 / gp1, gamma / gm1);
  // cd_nozzle_mass_flow is monotone decreasing in P_py: bisect for the root.
  double lo = Pt * r_crit;
  double hi = Pt * (1.0 - 1e-12);
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (lo + hi);
    double f = ejector_cd_nozzle_mass_flow(Pt, Tt, mid, area_throat, area_exit, gamma, r_gas,
                                           eta, eps_frac);
    if (f > m_dot) {
      lo = mid; // too much flow -> raise P_py
    } else {
      hi = mid;
    }
  }
  double p_py = 0.5 * (lo + hi);
  // Implicit function theorem: F(P_py) = cd_nozzle(P_py) - m_dot = 0, so
  // dP_py/dx = -(dF/dx)/(dF/dP_py), using the C-D nozzle's own Jacobian.
  EjectorCDNozzleJacobian j = ejector_cd_nozzle_mass_flow_and_jacobian(
      Pt, Tt, p_py, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  double inv = 1.0 / j.dmdot_dp_py;
  return {p_py, inv, -j.dmdot_dp0 * inv, -j.dmdot_dt0 * inv};
}

EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    double p_py_val, double omega_val, double gamma, double recovery_efficiency,
    double eta_primary, double eta_secondary) {
  // Seeds: 0=p_g, 1=t_g, 2=p_e, 3=t_e, 4=p_py, 5=omega.
  using D = DualN<6>;
  D p_g = D::seed(p_g_val, 0);
  D t_g = D::seed(t_g_val, 1);
  D p_e = D::seed(p_e_val, 2);
  D t_e = D::seed(t_e_val, 3);
  D p_py = D::seed(p_py_val, 4);
  D gam = D::seed(omega_val, 5);

  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double k = gm1 / gamma;

  // Both streams expand isentropically to the common mixing static P_py:
  // lambda = sqrt((g+1)/(g-1) * (1 - (P_py/P0)^((g-1)/g))). Subsonic (< 1).
  //
  // Out of domain (#481): a solver probe can put P_py above a supply's Pt
  // (no expansion: lambda^2 < 0), at or below zero, or reverse a flow. The
  // arguments then continue through a C1 positive floor (dpos_floor), so the
  // discharge stays finite and smooth instead of NaN -- which otherwise
  // poisoned the whole blended residual over most of the low-primary range.
  // Exact wherever the closure is defined.
  constexpr double kRatioFloor = 1e-12;
  constexpr double kLambdaSqFloor = 1e-8;
  auto lambda_of = [&](const D& p0) {
    const D r = dpos_floor(p_py / p0, kRatioFloor);
    return dsqrt(dpos_floor((gp1 / gm1) * (1.0 - dpow(r, k)), kLambdaSqFloor));
  };
  // The velocities the nozzles deliver (#588): lambda = sqrt(eta) lambda_is,
  // each stream at the stagnation pressure that velocity implies at P_py,
  // p0 = p_supply pi(lambda_is) / pi(lambda) (= P_py / pi(lambda) in the
  // domain; = p_supply exactly at eta = 1, floors included).
  const double kpi = gamma / gm1;
  auto pi_of = [&](const D& lam) {
    return dpow(dpos_floor(1.0 - (gm1 / gp1) * lam * lam, kRatioFloor), kpi);
  };
  const D lambda1_is = lambda_of(p_g);
  const D lambda2_is = lambda_of(p_e);
  D lambda1 = std::sqrt(eta_primary) * lambda1_is;
  D lambda2 = std::sqrt(eta_secondary) * lambda2_is;
  if (eta_primary != 1.0) {
    p_g = p_g * pi_of(lambda1_is) / pi_of(lambda1);
  }
  if (eta_secondary != 1.0) {
    p_e = p_e * pi_of(lambda2_is) / pi_of(lambda2);
  }
  gam = dpos_floor(gam, kRatioFloor);
  D theta21 = t_e / t_g;

  // Kracik & Dvorak mixing (their Eqs. 7-13), identical to the critical path's
  // chain -- only lambda1 is subsonic here instead of supersonic.
  D z1 = lambda1 + 1.0 / lambda1;
  D z2 = lambda2 + 1.0 / lambda2;
  D z3 = (z1 + gam * dsqrt(theta21) * z2) / dsqrt((1.0 + gam) * (1.0 + gam * theta21));
  // Subsonic root (Eq. 12), (z3 - sqrt(z3^2 - 4)) / 2 in the cancellation-free
  // form: a floored lambda makes z3 large, where the difference rounds to 0.
  // z3 >= 2 analytically; the floor only catches rounding at z3 = 2.
  D lambda3 = 2.0 / (z3 + dsqrt(dpos_floor(z3 * z3 - 4.0, kRatioFloor)));

  double k7 = std::pow(gp1 / 2.0, 1.0 / gm1);
  auto q_of = [&](const D& lam) {
    return dpow(dpos_floor(1.0 - (gm1 / gp1) * lam * lam, kRatioFloor), 1.0 / gm1) * k7 * lam;
  };
  D q1 = q_of(lambda1);
  D q2 = q_of(lambda2);
  D q3 = q_of(lambda3);

  D p03 = p_g * dsqrt((1.0 + gam) * (1.0 + gam * theta21)) /
          (1.0 + (p_g / p_e) * gam * dsqrt(theta21) * q1 / q2) * q1 / q3;
  D out = recovery_efficiency * p03;
  return {out.v, out.d[0], out.d[1], out.d[2], out.d[3], out.d[4], out.d[5]};
}

// -----------------------------------------------------------------------------
// Forward-mode dual number carrying a value plus its 4 partial derivatives
// w.r.t. (p_g, t_g, p_e, t_e). Every operator below implements the textbook
// analytic differentiation rule (sum/product/quotient/chain rule) -- this is
// bookkeeping for hand-derived-and-sympy-cross-checked analytic derivatives
// (see scripts/derive_ejector_jacobian.py), not automatic differentiation
// via a library and not finite differences. Writing the physics chain in
// terms of Dual4 lets this mirror _ejector_huang1999.py's Python almost
// line-for-line while the compiler enforces correct chain-rule composition
// at every step, rather than 80+ individually hand-tracked partials.
// -----------------------------------------------------------------------------
namespace {

struct Dual4 {
  double v = 0.0;
  double dpg = 0.0, dtg = 0.0, dpe = 0.0, dte = 0.0;

  static Dual4 constant(double c) { return Dual4{c, 0.0, 0.0, 0.0, 0.0}; }
};

Dual4 operator+(const Dual4& a, const Dual4& b) {
  return {a.v + b.v, a.dpg + b.dpg, a.dtg + b.dtg, a.dpe + b.dpe, a.dte + b.dte};
}
Dual4 operator+(const Dual4& a, double c) { return {a.v + c, a.dpg, a.dtg, a.dpe, a.dte}; }
Dual4 operator+(double c, const Dual4& a) { return a + c; }

Dual4 operator-(const Dual4& a, double c) { return {a.v - c, a.dpg, a.dtg, a.dpe, a.dte}; }
Dual4 operator-(const Dual4& a, const Dual4& b) {
  return {a.v - b.v, a.dpg - b.dpg, a.dtg - b.dtg, a.dpe - b.dpe, a.dte - b.dte};
}
Dual4 operator-(double c, const Dual4& a) { return {c - a.v, -a.dpg, -a.dtg, -a.dpe, -a.dte}; }

Dual4 operator*(const Dual4& a, const Dual4& b) {
  return {a.v * b.v, a.dpg * b.v + a.v * b.dpg, a.dtg * b.v + a.v * b.dtg,
          a.dpe * b.v + a.v * b.dpe, a.dte * b.v + a.v * b.dte};
}
Dual4 operator*(const Dual4& a, double c) { return {a.v * c, a.dpg * c, a.dtg * c, a.dpe * c, a.dte * c}; }
Dual4 operator*(double c, const Dual4& a) { return a * c; }

Dual4 operator/(const Dual4& a, const Dual4& b) {
  double inv = 1.0 / b.v;
  double inv2 = inv * inv;
  return {a.v * inv, (a.dpg * b.v - a.v * b.dpg) * inv2, (a.dtg * b.v - a.v * b.dtg) * inv2,
          (a.dpe * b.v - a.v * b.dpe) * inv2, (a.dte * b.v - a.v * b.dte) * inv2};
}
Dual4 operator/(const Dual4& a, double c) { return a * (1.0 / c); }
Dual4 operator/(double c, const Dual4& a) { return Dual4::constant(c) / a; }

Dual4 dsqrt(const Dual4& a) {
  double v = std::sqrt(a.v);
  double coef = 0.5 / v;
  return {v, a.dpg * coef, a.dtg * coef, a.dpe * coef, a.dte * coef};
}

// x^c for a CONSTANT real exponent c. Every exponent in this chain is a
// function of gamma alone (a fixed parameter here, not a Newton unknown),
// so this covers every power in the chain -- no Dual4^Dual4 is ever needed.
Dual4 dpow(const Dual4& a, double c) {
  double v = std::pow(a.v, c);
  double coef = c * std::pow(a.v, c - 1.0);
  return {v, a.dpg * coef, a.dtg * coef, a.dpe * coef, a.dte * coef};
}

// C1 positive floor (dpos_floor in dual_number.h, for Dual4): exact for
// a >= eps, eps^2 / (2 eps - a) below. Keeps the critical closure finite when
// a probe puts it outside its domain (#481).
Dual4 dpos_floor4(const Dual4& a, double eps) {
  if (a.v >= eps) {
    return a;
  }
  const double den = 2.0 * eps - a.v;
  const double slope = eps * eps / (den * den);
  return {eps * eps / den, a.dpg * slope, a.dtg * slope, a.dpe * slope, a.dte * slope};
}

// Out of the critical closure's domain -- the secondary's choke pressure at or
// above the primary's nozzle exit pressure (M_py^2 <= 0), or a primary core
// that fills the mixing area (A_sy <= 0) -- the arguments continue through
// dpos_floor4: exact wherever the closure is defined. The former hard switch
// M_py = 0 there divided by zero (A_py = inf) and the blended element residual
// went NaN.
constexpr double kMachSqFloor = 1e-8;
constexpr double kAreaFloorFrac = 1e-6;  // of the mixing area ratio

} // namespace

EjectorEntrainmentJacobian ejector_entrainment_ratio_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    const EjectorGeometry& geom, double gamma) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double g_over_gm1 = gamma / gm1;

  Dual4 p_g{p_g_val, 1.0, 0.0, 0.0, 0.0};
  Dual4 t_g{t_g_val, 0.0, 1.0, 0.0, 0.0};
  Dual4 p_e{p_e_val, 0.0, 0.0, 1.0, 0.0};
  Dual4 t_e{t_e_val, 0.0, 0.0, 0.0, 1.0};

  // Primary flow: nozzle-exit Mach (Eq. 2) -- exact zero derivative w.r.t.
  // p_g/t_g/p_e/t_e, see include/ejector.h's header comment.
  double mach_p1_val = ejector_mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, gamma);
  Dual4 mach_p1 = Dual4::constant(mach_p1_val);

  Dual4 p_p1 = p_g / dpow(1.0 + 0.5 * gm1 * mach_p1 * mach_p1, g_over_gm1);

  // Entrained flow chokes at y-y (M_sy = 1); pressure matching P_py = P_sy.
  double k2 = std::pow(gp1 / 2.0, g_over_gm1);
  Dual4 p_sy = p_e / k2;
  Dual4 p_py = p_sy;

  // Invert Eq. 4 for the primary-core Mach at y-y.
  Dual4 core = (1.0 + 0.5 * gm1 * mach_p1 * mach_p1) * dpow(p_py / p_p1, -gm1 / gamma);
  Dual4 mach_py_sq = (core - 1.0) / (0.5 * gm1);
  Dual4 mach_py = dsqrt(dpos_floor4(mach_py_sq, kMachSqFloor));

  // Areas at y-y (Eqs. 5, 8).
  Dual4 area_py =
      ejector::phi_p * (1.0 / mach_py) *
      dpow(2.0 / gp1 * (1.0 + 0.5 * gm1 * mach_py * mach_py), gp1 / (2.0 * gm1));
  Dual4 area_sy = dpos_floor4(geom.area_ratio_mix - area_py, kAreaFloorFrac * geom.area_ratio_mix);

  Dual4 omega = (p_e / p_g) * area_sy * dsqrt(t_g / t_e) * std::sqrt(ejector::eta_s / ejector::eta_p);

  EjectorEntrainmentResult value{
      omega.v,       mach_p1.v, mach_py.v, area_py.v, area_sy.v, p_py.v,
  };
  return {value, omega.dpg, omega.dtg, omega.dpe, omega.dte};
}

EjectorEntrainmentResult ejector_entrainment_ratio(double p_g, double t_g, double p_e,
                                                    double t_e, const EjectorGeometry& geom,
                                                    double gamma) {
  return ejector_entrainment_ratio_and_jacobian(p_g, t_g, p_e, t_e, geom, gamma).value;
}

namespace {

// Kracik & Dvorak's mixing closure from the critical y-y state, mixed at the
// critical entrainment ratio (P_c*) or at zero entrainment (the dead-head
// P_b0, their gam = 0): ONE chain, so the two anchors of the subcritical droop
// cannot drift apart (1:1 with _yy_lambda_state + _mixed_flow_stagnation).
EjectorCriticalPressureJacobian mixed_stagnation_from_yy(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    const EjectorGeometry& geom, double gamma, double r_gas, double recovery_efficiency,
    bool dead_head) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double g_over_gm1 = gamma / gm1;

  Dual4 p_g{p_g_val, 1.0, 0.0, 0.0, 0.0};
  Dual4 t_g{t_g_val, 0.0, 1.0, 0.0, 0.0};
  Dual4 p_e{p_e_val, 0.0, 0.0, 1.0, 0.0};
  Dual4 t_e{t_e_val, 0.0, 0.0, 0.0, 1.0};

  double mach_p1_val = ejector_mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, gamma);
  Dual4 mach_p1 = Dual4::constant(mach_p1_val);

  Dual4 p_p1 = p_g / dpow(1.0 + 0.5 * gm1 * mach_p1 * mach_p1, g_over_gm1);
  double k2 = std::pow(gp1 / 2.0, g_over_gm1);
  Dual4 p_sy = p_e / k2;
  Dual4 p_py = p_sy;

  Dual4 core = (1.0 + 0.5 * gm1 * mach_p1 * mach_p1) * dpow(p_py / p_p1, -gm1 / gamma);
  Dual4 mach_py_sq = (core - 1.0) / (0.5 * gm1);
  Dual4 mach_py = dsqrt(dpos_floor4(mach_py_sq, kMachSqFloor));

  Dual4 area_py =
      ejector::phi_p * (1.0 / mach_py) *
      dpow(2.0 / gp1 * (1.0 + 0.5 * gm1 * mach_py * mach_py), gp1 / (2.0 * gm1));
  Dual4 area_sy = dpos_floor4(geom.area_ratio_mix - area_py, kAreaFloorFrac * geom.area_ratio_mix);

  Dual4 omega = (p_e / p_g) * area_sy * dsqrt(t_g / t_e) * std::sqrt(ejector::eta_s / ejector::eta_p);

  // Static temperatures and velocities at y-y (Huang Eqs. 9-10, 13-14).
  Dual4 t_py = t_g / (1.0 + 0.5 * gm1 * mach_py * mach_py);
  Dual4 t_sy = t_e / (gp1 / 2.0);
  Dual4 v_py = mach_py * dsqrt(gamma * r_gas * t_py);
  Dual4 v_sy = dsqrt(gamma * r_gas * t_sy); // M_sy = 1

  Dual4 c1_star = dsqrt(2.0 * gamma / gp1 * r_gas * t_g);
  Dual4 c2_star = dsqrt(2.0 * gamma / gp1 * r_gas * t_e);
  Dual4 lambda1 = v_py / c1_star;
  Dual4 lambda2 = v_sy / c2_star; // == 1 identically (M_sy = 1); kept computed for value parity

  Dual4 theta21 = t_e / t_g;
  Dual4 gam = dead_head ? Dual4::constant(0.0) : omega;

  Dual4 z1 = lambda1 + 1.0 / lambda1;
  Dual4 z2 = lambda2 + 1.0 / lambda2;
  Dual4 z3 = (z1 + gam * dsqrt(theta21) * z2) / dsqrt((1.0 + gam) * (1.0 + gam * theta21));
  // Subsonic root (Eq. 12) in the cancellation-free form 2/(z3 + sqrt(z3^2-4)).
  Dual4 lambda3 = 2.0 / (z3 + dsqrt(dpos_floor4(z3 * z3 - 4.0, 1e-12)));

  double gm1_local = gm1;
  double gp1_local = gp1;
  double k7 = std::pow(gp1_local / 2.0, 1.0 / gm1_local);
  auto q_of = [&](const Dual4& lam) {
    Dual4 term = dpos_floor4(1.0 - (gm1_local / gp1_local) * lam * lam, 1e-12);
    return dpow(term, 1.0 / gm1_local) * k7 * lam;
  };
  Dual4 q1 = q_of(lambda1);
  Dual4 q2 = q_of(lambda2); // == 1 identically (lambda2 == 1)
  Dual4 q3 = q_of(lambda3);

  // Mixed flow's stagnation pressure and temperature after mixing (Eqs. 7, 11).
  Dual4 p03 = p_g * dsqrt((1.0 + gam) * (1.0 + gam * theta21)) /
              (1.0 + (p_g / p_e) * gam * dsqrt(theta21) * q1 / q2) * q1 / q3;
  Dual4 t03 = t_g * (1.0 + gam * theta21) / (1.0 + gam);

  double mach3_sq = 2.0 * lambda3.v * lambda3.v / (gp1 - gm1 * lambda3.v * lambda3.v);
  double mach3 = std::sqrt(mach3_sq);

  EjectorCriticalPressureResult value{
      recovery_efficiency * p03.v, p03.v, t03.v, lambda3.v, mach3,
  };
  return {value, recovery_efficiency * p03.dpg, recovery_efficiency * p03.dtg,
          recovery_efficiency * p03.dpe, recovery_efficiency * p03.dte};
}

} // namespace

EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return mixed_stagnation_from_yy(p_g, t_g, p_e, t_e, geom, gamma, r_gas, recovery_efficiency,
                                  false);
}

EjectorCriticalPressureJacobian ejector_dead_head_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return mixed_stagnation_from_yy(p_g, t_g, p_e, t_e, geom, gamma, r_gas, recovery_efficiency,
                                  true);
}

namespace {

// The critical y-y state (Huang Eqs. 1-8), as Dual4 in (p_g, t_g, p_e, t_e):
// the same chain ejector_entrainment_ratio_and_jacobian and
// mixed_stagnation_from_yy evaluate.
struct YYState4 {
  Dual4 p_g, t_g, p_e, t_e;
  Dual4 p_py, mach_py, omega;
};

YYState4 critical_yy_state4(double p_g_val, double t_g_val, double p_e_val, double t_e_val,
                            const EjectorGeometry& geom, double gamma) {
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  const double g_over_gm1 = gamma / gm1;
  YYState4 s;
  s.p_g = Dual4{p_g_val, 1.0, 0.0, 0.0, 0.0};
  s.t_g = Dual4{t_g_val, 0.0, 1.0, 0.0, 0.0};
  s.p_e = Dual4{p_e_val, 0.0, 0.0, 1.0, 0.0};
  s.t_e = Dual4{t_e_val, 0.0, 0.0, 0.0, 1.0};
  const double mach_p1 = ejector_mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, gamma);
  const double exit_factor = 1.0 + 0.5 * gm1 * mach_p1 * mach_p1;
  Dual4 p_p1 = s.p_g / std::pow(exit_factor, g_over_gm1);
  s.p_py = s.p_e / std::pow(gp1 / 2.0, g_over_gm1);
  Dual4 core = exit_factor * dpow(s.p_py / p_p1, -gm1 / gamma);
  s.mach_py = dsqrt(dpos_floor4((core - 1.0) / (0.5 * gm1), kMachSqFloor));
  Dual4 area_py = ejector::phi_p * (1.0 / s.mach_py) *
                  dpow(2.0 / gp1 * (1.0 + 0.5 * gm1 * s.mach_py * s.mach_py), gp1 / (2.0 * gm1));
  Dual4 area_sy =
      dpos_floor4(geom.area_ratio_mix - area_py, kAreaFloorFrac * geom.area_ratio_mix);
  s.omega = (s.p_e / s.p_g) * area_sy * dsqrt(s.t_g / s.t_e) *
            std::sqrt(ejector::eta_s / ejector::eta_p);
  return s;
}

} // namespace

double ejector_huang_mixing_coefficient(double area_ratio_mix) {
  if (area_ratio_mix > ejector::phi_m_upper_area_ratio) {
    return ejector::phi_m_large;
  }
  return area_ratio_mix >= ejector::phi_m_lower_area_ratio ? ejector::phi_m_mid
                                                           : ejector::phi_m_small;
}

EjectorCriticalPressureJacobian ejector_huang_critical_back_pressure_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  YYState4 s = critical_yy_state4(p_g_val, t_g_val, p_e_val, t_e_val, geom, gamma);

  // Static temperatures and velocities at y-y (Huang Eqs. 9-10, 13-14).
  Dual4 t_py = s.t_g / (1.0 + 0.5 * gm1 * s.mach_py * s.mach_py);
  Dual4 t_sy = s.t_e / (gp1 / 2.0);
  Dual4 v_py = s.mach_py * dsqrt(gamma * r_gas * t_py);
  Dual4 v_sy = dsqrt(gamma * r_gas * t_sy);  // M_sy = 1

  // Constant-pressure mixing with phi_m (Eq. 11) and energy (Eq. 12).
  const double phi_m = ejector_huang_mixing_coefficient(geom.area_ratio_mix);
  Dual4 v_m = phi_m * (v_py + s.omega * v_sy) / (1.0 + s.omega);
  const double cp = gamma * r_gas / gm1;
  Dual4 t03 = (s.t_g + s.omega * s.t_e) / (1.0 + s.omega);
  Dual4 t_m = dpos_floor4(t03 - 0.5 * v_m * v_m / cp, 1e-6 * t03.v);
  Dual4 mach_m_sq = v_m * v_m / (gamma * r_gas * t_m);

  // Normal shock in the mixed flow (Eqs. 16-17), C2 through M_m = 1.
  Dual4 p3 = s.p_py;
  Dual4 mach3_sq = mach_m_sq;
  if (mach_m_sq.v > 1.0) {
    p3 = s.p_py * (1.0 + 2.0 * gamma / gp1 * (mach_m_sq - 1.0));
    mach3_sq = (1.0 + 0.5 * gm1 * mach_m_sq) / (gamma * mach_m_sq - 0.5 * gm1);
  }
  // Isentropic diffuser (Eq. 18).
  Dual4 p_c = p3 * dpow(1.0 + 0.5 * gm1 * mach3_sq, gamma / gm1);

  const double lambda3 = std::sqrt(gp1 * mach3_sq.v / (2.0 + gm1 * mach3_sq.v));
  EjectorCriticalPressureResult value{
      recovery_efficiency * p_c.v, p_c.v, t03.v, lambda3, std::sqrt(mach3_sq.v),
  };
  return {value, recovery_efficiency * p_c.dpg, recovery_efficiency * p_c.dtg,
          recovery_efficiency * p_c.dpe, recovery_efficiency * p_c.dte};
}

EjectorCriticalPressureJacobian ejector_sudden_expansion_dead_head_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  (void)t_e_val;  // the secondary is dead-headed
  const double gm1 = gamma - 1.0;
  const double gp1 = gamma + 1.0;
  const double k = 0.5 * gm1;
  const double a1 = geom.area_ratio_nozzle;
  const double a3 = geom.area_ratio_mix;
  Dual4 p_g{p_g_val, 1.0, 0.0, 0.0, 0.0};
  Dual4 t_g{t_g_val, 0.0, 1.0, 0.0, 0.0};
  Dual4 p_e{p_e_val, 0.0, 0.0, 1.0, 0.0};

  const double mach_p1 = ejector_mach_from_area_ratio_supersonic(a1, gamma);
  Dual4 p_p1 = p_g / std::pow(1.0 + k * mach_p1 * mach_p1, gamma / gm1);
  // Choked primary mass flow per unit throat area (isentropic).
  Dual4 m = p_g * dsqrt(gamma / (r_gas * t_g)) * std::pow(2.0 / gp1, gp1 / (2.0 * gm1));
  Dual4 impulse = p_p1 * (a1 * (1.0 + gamma * mach_p1 * mach_p1)) + p_e * (a3 - a1);
  Dual4 phi = impulse / (m * dsqrt(r_gas * t_g / gamma));
  // Subsonic root of (k Phi^2 - g^2) X^2 + (Phi^2 - 2g) X - 1 = 0, X = M^2.
  Dual4 b = phi * phi - 2.0 * gamma;
  Dual4 d = dpos_floor4(b * b + 4.0 * (k * phi * phi - gamma * gamma),
                        ejector::dead_head_discriminant_floor);
  Dual4 x = 2.0 / (b + dsqrt(d));
  Dual4 t3 = t_g / (1.0 + k * x);
  Dual4 p3 = m * dsqrt(r_gas * t3 / gamma) / (a3 * dsqrt(x));
  Dual4 p_b0 = p3 * dpow(1.0 + k * x, gamma / gm1);

  const double lambda3 = std::sqrt(gp1 * x.v / (2.0 + gm1 * x.v));
  EjectorCriticalPressureResult value{
      recovery_efficiency * p_b0.v, p_b0.v, t_g.v, lambda3, std::sqrt(x.v),
  };
  return {value, recovery_efficiency * p_b0.dpg, recovery_efficiency * p_b0.dtg,
          recovery_efficiency * p_b0.dpe, recovery_efficiency * p_b0.dte};
}

EjectorCriticalPressureResult ejector_critical_back_pressure(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return ejector_critical_back_pressure_and_jacobian(p_g, t_g, p_e, t_e, geom, gamma, r_gas,
                                                      recovery_efficiency)
      .value;
}

// -----------------------------------------------------------------------------
// Whole-element (f, J): the 4-row operating-regime residual system.
//
// A transliteration of EjectorElement.residuals() (ejector_element.py): the
// scalar closures above are chained into a forward-mode DualN<9> exactly as the
// Python `_chain` wraps each closure's (value, partials), then blended by the
// two smootherstep weights. Kept in 1:1 correspondence so the two can be
// cross-validated (tests/test_ejector_jacobians.cpp FD-checks this directly).
// -----------------------------------------------------------------------------
namespace {

// _chain: value + sum(partial * grad(input dual)) -- 1:1 with Python _chain.
template <int N>
DualN<N> ejector_chain(double value,
                       std::initializer_list<std::pair<double, DualN<N>>> pairs) {
  DualN<N> r = DualN<N>::constant(value);
  for (const auto& pr : pairs) {
    for (int i = 0; i < N; ++i) r.d[i] += pr.first * pr.second.d[i];
  }
  return r;
}

// C1 clamped smootherstep saturating to EXACTLY 0/1 outside [lo, hi] -- 1:1
// with Python _smootherstep (so s_choke = 1 exactly at the choke threshold and
// critical mode reduces exactly).
template <int N>
DualN<N> ejector_smootherstep(const DualN<N>& t, double lo, double hi) {
  DualN<N> x = (t - lo) / (hi - lo);
  if (x.v <= 0.0) return DualN<N>::constant(0.0);
  if (x.v >= 1.0) return DualN<N>::constant(1.0);
  double val = x.v * x.v * x.v * (x.v * (x.v * 6.0 - 15.0) + 10.0);
  double dval = 30.0 * x.v * x.v * (x.v - 1.0) * (x.v - 1.0);
  return ejector_chain<N>(val, {{dval, x}});
}

// The subcritical droop fades in over (P_b0 - P_c*)/P_c* in [lo, hi].
constexpr double kDroopWindowLo = 1e-3;
constexpr double kDroopWindowHi = 1e-2;

} // namespace

EjectorElementResidualJacobian ejector_element_residuals_and_jacobian(
    double mp_v, double ms_v, double mo_v, double p_g_v, double t_g_v,
    double p_e_v, double t_e_v, double p_out_v, double p_py_v,
    const EjectorGeometry& geom, double area_throat, double area_nozzle_exit,
    double area_secondary, double gamma, double r_gas, double eta_primary,
    double eta_secondary, double recovery_efficiency, double eps_frac,
    double s_choke_lo, double s_choke_hi, double s_sub_lo, double s_sub_hi,
    EjectorAnchorClosure anchor_closure) {
  using D = DualN<9>;
  const bool huang = anchor_closure == EjectorAnchorClosure::huang1999;

  // Seed the nine Newton unknowns (see header for the index convention).
  const D mp = D::seed(mp_v, 0);
  const D ms = D::seed(ms_v, 1);
  const D mdot_out = D::seed(mo_v, 2);
  const D p_g = D::seed(p_g_v, 3);
  const D t_g = D::seed(t_g_v, 4);
  const D p_e = D::seed(p_e_v, 5);
  const D t_e = D::seed(t_e_v, 6);
  const D p_out = D::seed(p_out_v, 7);
  const D p_py = D::seed(p_py_v, 8); // owned unknown == mixing-plane static

  // s_choke = smootherstep(mp / cap), cap = choked_mass_flow(Pt_p).
  auto [cap_v, dcap_dp0, dcap_dt0] =
      ejector_choked_mass_flow_and_jacobian(p_g_v, t_g_v, area_throat, gamma, r_gas, eta_primary);
  D cap = ejector_chain<9>(cap_v, {{dcap_dp0, p_g}, {dcap_dt0, t_g}});
  D s = ejector_smootherstep<9>(mp / cap, s_choke_lo, s_choke_hi);

  // R0: primary C-D nozzle (choked/unchoked via the closure's own smooth-min).
  auto cd = ejector_cd_nozzle_mass_flow_and_jacobian(
      p_g_v, t_g_v, p_py_v, area_throat, area_nozzle_exit, gamma, r_gas, eta_primary, eps_frac);
  D cd_d = ejector_chain<9>(
      cd.mdot, {{cd.dmdot_dp0, p_g}, {cd.dmdot_dt0, t_g}, {cd.dmdot_dp_py, p_py}});
  D r0 = mp - cd_d;

  // Blend the critical (choked-primary) and jet-pump (unchoked) closures by
  // s = s_choke. The critical closures return NaN for an unchoked primary, and
  // the jet-pump closure is symmetric at the choked end; smootherstep is flat
  // at s in {0, 1} so the inactive branch contributes zero, but 0 * NaN = NaN
  // would still poison the value/Jacobian -- so evaluate each branch ONLY where
  // its weight is nonzero, leaving a hard-zero dual otherwise (C1-consistent
  // with the flat endpoint).
  // ACTUAL ratio for the discharge (R1 pins ms itself). A probe can drive
  // either flow to zero or reverse it: both are floored C1 at 1e-9 of the
  // choked flow, so omega stays finite and positive (#481).
  const double m_floor = 1e-9 * cap.v;
  D omega_eff = dpos_floor(ms, m_floor) / dpos_floor(mp, m_floor);
  const bool crit_active = s.v > 0.0;
  const bool jet_active = s.v < 1.0;

  D omega_crit = D::constant(0.0);
  D p03_crit = D::constant(0.0);
  D s_sub = D::constant(1.0);
  // Subcritical droop (Phase 1A, #619): between P_c* and the dead-head P_b0
  // the choked ejector's entrainment falls linearly to 0 (omega_choked) and the
  // outlet floats; the outlet pin switches in at P_b0 (s_b0), not at P_c*.
  D omega_choked = D::constant(0.0);
  D p03_dead = D::constant(0.0);
  D s_b0 = D::constant(1.0);
  if (crit_active) {
    auto entr = ejector_entrainment_ratio_and_jacobian(p_g_v, t_g_v, p_e_v, t_e_v, geom, gamma);
    omega_crit = ejector_chain<9>(entr.value.omega,
                                  {{entr.domega_dp_g, p_g},
                                   {entr.domega_dt_g, t_g},
                                   {entr.domega_dp_e, p_e},
                                   {entr.domega_dt_e, t_e}});
    auto pc = huang ? ejector_huang_critical_back_pressure_and_jacobian(
                          p_g_v, t_g_v, p_e_v, t_e_v, geom, gamma, r_gas, 1.0)
                    : ejector_critical_back_pressure_and_jacobian(p_g_v, t_g_v, p_e_v, t_e_v,
                                                                  geom, gamma, r_gas, 1.0);
    p03_crit = ejector_chain<9>(pc.value.p_mixed_stagnation_pa,
                                {{pc.dpc_dp_g, p_g},
                                 {pc.dpc_dt_g, t_g},
                                 {pc.dpc_dp_e, p_e},
                                 {pc.dpc_dt_e, t_e}});
    s_sub = ejector_smootherstep<9>(p_out / (recovery_efficiency * p03_crit), s_sub_lo, s_sub_hi);

    // The droop needs a window P_b0 > P_c*; it collapses in the M_py < 1
    // corner (an unchoked primary, the jet pump's). d fades the droop out
    // there, back to the plateau-and-pin-at-P_c* form, continuously: P_b0
    // depends on Newton unknowns, so a hard switch would be a jump in F.
    auto pb = huang ? ejector_sudden_expansion_dead_head_and_jacobian(
                          p_g_v, t_g_v, p_e_v, t_e_v, geom, gamma, r_gas, 1.0)
                    : ejector_dead_head_back_pressure_and_jacobian(p_g_v, t_g_v, p_e_v, t_e_v,
                                                                   geom, gamma, r_gas, 1.0);
    D p03_b0 = ejector_chain<9>(pb.value.p_mixed_stagnation_pa,
                                {{pb.dpc_dp_g, p_g},
                                 {pb.dpc_dt_g, t_g},
                                 {pb.dpc_dp_e, p_e},
                                 {pb.dpc_dt_e, t_e}});
    D window = (p03_b0 - p03_crit) / p03_crit;
    D d = ejector_smootherstep<9>(window, kDroopWindowLo, kDroopWindowHi);
    D s_b0_live = D::constant(0.0);
    D omega_sub = omega_crit;
    if (d.v > 0.0) {
      // Tier-1 chord (P_c*, omega_crit) -> (P_b0, 0) in the outlet Pt;
      // window >= kDroopWindowLo keeps the division bounded.
      omega_sub = omega_crit * (recovery_efficiency * p03_b0 - p_out) /
                  (recovery_efficiency * (p03_b0 - p03_crit));
      s_b0_live = ejector_smootherstep<9>(p_out / (recovery_efficiency * p03_b0), s_sub_lo,
                                          s_sub_hi);
    }
    // Plateau below P_c* exactly (s_sub = 0), the chord above it: critical
    // results are unchanged, unlike a smooth-min, which shaves omega_crit.
    omega_choked = omega_crit + d * s_sub * (omega_sub - omega_crit);
    s_b0 = d * s_b0_live + (1.0 - d) * s_sub;
    p03_dead = d * p03_b0 + (1.0 - d) * p03_crit;
  }

  D ms_sec = D::constant(0.0);
  D p03_jet = D::constant(0.0);
  if (jet_active) {
    auto sec = ejector_cd_nozzle_mass_flow_and_jacobian(
        p_e_v, t_e_v, p_py_v, area_secondary, area_secondary, gamma, r_gas, eta_secondary,
        eps_frac);
    ms_sec = ejector_chain<9>(
        sec.mdot, {{sec.dmdot_dp0, p_e}, {sec.dmdot_dt0, t_e}, {sec.dmdot_dp_py, p_py}});
    auto jp = ejector_jetpump_discharge_and_jacobian(
        p_g_v, t_g_v, p_e_v, t_e_v, p_py_v, omega_eff.v, gamma, 1.0, eta_primary,
        eta_secondary);
    p03_jet = ejector_chain<9>(jp.p03,
                               {{jp.dp03_dp_g, p_g},
                                {jp.dp03_dt_g, t_g},
                                {jp.dp03_dp_e, p_e},
                                {jp.dp03_dt_e, t_e},
                                {jp.dp03_dp_py, p_py},
                                {jp.dp03_domega, omega_eff}});
  }

  // R1: blended entrainment  ms = s*omega_crit*mp + (1-s)*ms_sec.
  D ms_model = s * omega_choked * mp + (1.0 - s) * ms_sec;
  D r1 = ms - ms_model;

  // R2: mass conservation.
  D r2 = mdot_out - mp - ms;

  // R3: regime-dependent closure of the owned unknown (see ejector_element.py
  // for the w_pin table). discharge = recovery * blend(P_c*, P03_jet) by
  // s_choke; w_pin = 1 - s*(1 - s_sub) is the weight on the outlet-pin form.
  // Choked: P_py = recovery*P_c* (a diagnostic) while the outlet floats,
  // below P_b0; the outlet pinned at recovery*P_b0 above it, which excludes
  // the spurious double-choked root (Pt_primary below the back pressure is
  // only possible above P_b0, since P_b0 < Pt_primary).
  D discharge = recovery_efficiency *
                (s * ((1.0 - s_b0) * p03_crit + s_b0 * p03_dead) + (1.0 - s) * p03_jet);
  D w_pin = 1.0 - s * (1.0 - s_b0);
  D lhs = w_pin * p_out + (1.0 - w_pin) * p_py;
  D r3 = lhs - discharge;

  EjectorElementResidualJacobian out;
  const D rows[4] = {r0, r1, r2, r3};
  for (int i = 0; i < 4; ++i) {
    out.residuals[i] = rows[i].v;
    for (int j = 0; j < 9; ++j) out.jacobian[i][j] = rows[i].d[j];
  }
  return out;
}

EjectorJetPumpIncompressibleResult ejector_jetpump_incompressible_residuals_and_jacobian(
    double mp, double ms, double mo, double pt_p, double v_p, double pt_s,
    double v_s, double pt_o, double p3, double area_nozzle_exit,
    double area_secondary, double area_mixing, double k_nozzle,
    double k_suction, double k_throat, double k_diffuser, double recovery) {
  using D = DualN<9>;
  const D Mp = D::seed(mp, 0), Ms = D::seed(ms, 1), Mo = D::seed(mo, 2);
  const D Ptp = D::seed(pt_p, 3), Vp = D::seed(v_p, 4);
  const D Pts = D::seed(pt_s, 5), Vs = D::seed(v_s, 6);
  const D Pto = D::seed(pt_o, 7), P3 = D::seed(p3, 8);
  const double Ae = area_nozzle_exit, As = area_secondary, A3 = area_mixing;

  // m|m| as m * m * sgn(m): exact derivative 2|m|, C1 through zero.
  auto msq = [](const D &m) { return m * m * (m.v >= 0.0 ? 1.0 : -1.0); };
  const D Mt = Mp + Ms;
  const double sgn_t = (Mt.v >= 0.0) ? 1.0 : -1.0;
  // mt|mt| v_m = |mt| (mp v_p + ms v_s)
  const D flux_t = Mt * sgn_t * (Mp * Vp + Ms * Vs);

  const D R0 = Ptp - P3 - (1.0 + k_nozzle) * 0.5 * Vp * msq(Mp) / (Ae * Ae);
  const D R1 = Pts - P3 - (1.0 + k_suction) * 0.5 * Vs * msq(Ms) / (As * As);
  const D R2 = Mo - Mp - Ms;
  const D q4 = 0.5 * flux_t / (A3 * A3);
  const D p4 = P3 + (Vp * msq(Mp) / Ae + Vs * msq(Ms) / As - flux_t / A3) / A3 - k_throat * q4;
  const D R3 = Pto - recovery * (p4 + (1.0 - k_diffuser) * q4);

  EjectorJetPumpIncompressibleResult out;
  const D rows[4] = {R0, R1, R2, R3};
  for (int r = 0; r < 4; ++r) {
    out.residuals[r] = rows[r].v;
    for (int k = 0; k < 9; ++k) {
      out.jacobian[r][k] = rows[r].d[k];
    }
  }
  out.p4 = p4.v;
  out.q4 = q4.v;
  return out;
}

} // namespace combaero::solver
