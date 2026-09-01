// 1-D critical-mode supersonic ejector model. See include/ejector.h for
// paper references and the Jacobian-scope note.

#include "ejector.h"
#include <array>
#include <cmath>

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

double ejector_choked_mass_flow(double p0, double t0, double area_throat,
                                double gamma, double r_gas, double eta) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double choke_const = std::sqrt(gamma / r_gas * std::pow(2.0 / gp1, gp1 / gm1));
  return (p0 * area_throat / std::sqrt(t0)) * choke_const * std::sqrt(eta);
}

std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian(double p0, double t0, double area_throat,
                                      double gamma, double r_gas, double eta) {
  double mdot = ejector_choked_mass_flow(p0, t0, area_throat, gamma, r_gas, eta);
  // mdot = C * p0 / sqrt(t0)  =>  d/dp0 = mdot/p0,  d/dt0 = -mdot/(2*t0)
  return {mdot, mdot / p0, -0.5 * mdot / t0};
}

double ejector_cd_nozzle_mass_flow(double p0, double t0, double p_static_down,
                                   double area_throat, double area_exit,
                                   double gamma, double r_gas, double eta,
                                   double eps_frac) {
  double cap = ejector_choked_mass_flow(p0, t0, area_throat, gamma, r_gas, eta);
  double exit_flux = 0.0;
  double r = p_static_down / p0;
  if (r < 1.0) {
    double gm1 = gamma - 1.0;
    double gp1 = gamma + 1.0;
    double mach = std::sqrt(2.0 / gm1 * (std::pow(r, -gm1 / gamma) - 1.0));
    double flux = mach * std::pow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
    exit_flux = area_exit * p0 / std::sqrt(t0) * std::sqrt(gamma / r_gas) * flux * std::sqrt(eta);
  }
  double eps = eps_frac * cap;
  double diff = exit_flux - cap;
  return 0.5 * (exit_flux + cap - std::sqrt(diff * diff + eps * eps));
}

double ejector_q_lambda(double lam, double gamma) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  return std::pow(1.0 - (gm1 / gp1) * lam * lam, 1.0 / gm1) * std::pow(gp1 / 2.0, 1.0 / gm1) * lam;
}

// -----------------------------------------------------------------------------
// Generic forward-mode dual (N seeds) for the operating-regime closures, whose
// Newton unknowns include P_py and the outlet pressure in addition to the four
// thermodynamic inputs -- so the fixed 4-seed Dual4 above (kept untouched for
// the validated critical-mode path) does not fit. Same analytic chain-rule
// bookkeeping; DualN<N>::seed(v, i) sets the i-th partial to 1.
// -----------------------------------------------------------------------------
namespace {

template <int N> struct DualN {
  double v = 0.0;
  std::array<double, N> d{};
  static DualN constant(double c) {
    DualN r;
    r.v = c;
    return r;
  }
  static DualN seed(double c, int i) {
    DualN r;
    r.v = c;
    r.d[i] = 1.0;
    return r;
  }
};

template <int N> DualN<N> operator+(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v + b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] + b.d[i];
  return r;
}
template <int N> DualN<N> operator+(const DualN<N>& a, double c) {
  DualN<N> r = a;
  r.v += c;
  return r;
}
template <int N> DualN<N> operator+(double c, const DualN<N>& a) { return a + c; }
template <int N> DualN<N> operator-(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v - b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] - b.d[i];
  return r;
}
template <int N> DualN<N> operator-(const DualN<N>& a, double c) {
  DualN<N> r = a;
  r.v -= c;
  return r;
}
template <int N> DualN<N> operator-(double c, const DualN<N>& a) {
  DualN<N> r;
  r.v = c - a.v;
  for (int i = 0; i < N; ++i) r.d[i] = -a.d[i];
  return r;
}
template <int N> DualN<N> operator*(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  r.v = a.v * b.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
template <int N> DualN<N> operator*(const DualN<N>& a, double c) {
  DualN<N> r;
  r.v = a.v * c;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * c;
  return r;
}
template <int N> DualN<N> operator*(double c, const DualN<N>& a) { return a * c; }
template <int N> DualN<N> operator/(const DualN<N>& a, const DualN<N>& b) {
  DualN<N> r;
  double inv = 1.0 / b.v;
  double inv2 = inv * inv;
  r.v = a.v * inv;
  for (int i = 0; i < N; ++i) r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) * inv2;
  return r;
}
template <int N> DualN<N> operator/(const DualN<N>& a, double c) { return a * (1.0 / c); }
template <int N> DualN<N> operator/(double c, const DualN<N>& a) {
  DualN<N> r;
  double inv = 1.0 / a.v;
  r.v = c * inv;
  double coef = -c * inv * inv;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}
template <int N> DualN<N> dsqrt(const DualN<N>& a) {
  DualN<N> r;
  r.v = std::sqrt(a.v);
  double coef = 0.5 / r.v;
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}
template <int N> DualN<N> dpow(const DualN<N>& a, double c) {
  DualN<N> r;
  r.v = std::pow(a.v, c);
  double coef = c * std::pow(a.v, c - 1.0);
  for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * coef;
  return r;
}

} // namespace

EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian(
    double p0_val, double t0_val, double p_py_val, double area_throat,
    double area_exit, double gamma, double r_gas, double eta, double eps_frac) {
  // Seeds: 0 = p0, 1 = t0, 2 = p_py.
  using D = DualN<3>;
  D p0 = D::seed(p0_val, 0);
  D t0 = D::seed(t0_val, 1);
  D p_py = D::seed(p_py_val, 2);

  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double choke_const = std::sqrt(gamma / r_gas * std::pow(2.0 / gp1, gp1 / gm1));
  D cap = (choke_const * area_throat * std::sqrt(eta)) * p0 / dsqrt(t0);

  D exit_flux = D::constant(0.0);
  if (p_py_val < p0_val) {
    D r = p_py / p0;
    D mach_sq = (2.0 / gm1) * (dpow(r, -gm1 / gamma) - 1.0);
    D mach = dsqrt(mach_sq);
    D flux = mach * dpow(1.0 + 0.5 * gm1 * mach * mach, -gp1 / (2.0 * gm1));
    exit_flux = (area_exit * std::sqrt(gamma / r_gas) * std::sqrt(eta)) * p0 / dsqrt(t0) * flux;
  }

  // eps = eps_frac * cap depends on (p0, t0); differentiate it too (FD of the
  // Python value function does), else the Jacobian would be inconsistent.
  D eps = eps_frac * cap;
  D diff = exit_flux - cap;
  D mdot = 0.5 * (exit_flux + cap - dsqrt(diff * diff + eps * eps));
  return {mdot.v, mdot.d[0], mdot.d[1], mdot.d[2]};
}

EjectorCDExitStaticJacobian ejector_cd_nozzle_exit_static_and_jacobian(
    double m_dot, double p0, double t0, double area_throat, double area_exit,
    double gamma, double r_gas, double eta, double eps_frac) {
  double gm1 = gamma - 1.0;
  double gp1 = gamma + 1.0;
  double r_crit = std::pow(2.0 / gp1, gamma / gm1);
  // cd_nozzle_mass_flow is monotone decreasing in P_py: bisect for the root.
  double lo = p0 * r_crit;
  double hi = p0 * (1.0 - 1e-12);
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (lo + hi);
    double f = ejector_cd_nozzle_mass_flow(p0, t0, mid, area_throat, area_exit, gamma, r_gas,
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
      p0, t0, p_py, area_throat, area_exit, gamma, r_gas, eta, eps_frac);
  double inv = 1.0 / j.dmdot_dp_py;
  return {p_py, inv, -j.dmdot_dp0 * inv, -j.dmdot_dt0 * inv};
}

EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    double p_py_val, double omega_val, double gamma, double recovery_efficiency) {
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
  D lambda1 = dsqrt((gp1 / gm1) * (1.0 - dpow(p_py / p_g, k)));
  D lambda2 = dsqrt((gp1 / gm1) * (1.0 - dpow(p_py / p_e, k)));
  D theta21 = t_e / t_g;

  // Kracik & Dvorak mixing (their Eqs. 7-13), identical to the critical path's
  // chain -- only lambda1 is subsonic here instead of supersonic.
  D z1 = lambda1 + 1.0 / lambda1;
  D z2 = lambda2 + 1.0 / lambda2;
  D z3 = (z1 + gam * dsqrt(theta21) * z2) / dsqrt((1.0 + gam) * (1.0 + gam * theta21));
  D lambda3 = (z3 - dsqrt(z3 * z3 - 4.0)) * 0.5; // subsonic root (Eq. 12)

  double k7 = std::pow(gp1 / 2.0, 1.0 / gm1);
  auto q_of = [&](const D& lam) {
    return dpow(1.0 - (gm1 / gp1) * lam * lam, 1.0 / gm1) * k7 * lam;
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

Dual4 operator-(const Dual4& a, const Dual4& b) {
  return {a.v - b.v, a.dpg - b.dpg, a.dtg - b.dtg, a.dpe - b.dpe, a.dte - b.dte};
}
Dual4 operator-(const Dual4& a, double c) { return {a.v - c, a.dpg, a.dtg, a.dpe, a.dte}; }
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
  Dual4 mach_py = mach_py_sq.v > 0.0 ? dsqrt(mach_py_sq) : Dual4::constant(0.0);

  // Areas at y-y (Eqs. 5, 8).
  Dual4 area_py =
      ejector::phi_p * (1.0 / mach_py) *
      dpow(2.0 / gp1 * (1.0 + 0.5 * gm1 * mach_py * mach_py), gp1 / (2.0 * gm1));
  Dual4 area_sy = geom.area_ratio_mix - area_py;

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

EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian(
    double p_g_val, double t_g_val, double p_e_val, double t_e_val,
    const EjectorGeometry& geom, double gamma, double r_gas, double recovery_efficiency) {
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
  Dual4 mach_py = mach_py_sq.v > 0.0 ? dsqrt(mach_py_sq) : Dual4::constant(0.0);

  Dual4 area_py =
      ejector::phi_p * (1.0 / mach_py) *
      dpow(2.0 / gp1 * (1.0 + 0.5 * gm1 * mach_py * mach_py), gp1 / (2.0 * gm1));
  Dual4 area_sy = geom.area_ratio_mix - area_py;

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
  Dual4 gam = omega;

  Dual4 z1 = lambda1 + 1.0 / lambda1;
  Dual4 z2 = lambda2 + 1.0 / lambda2;
  Dual4 z3 = (z1 + gam * dsqrt(theta21) * z2) / dsqrt((1.0 + gam) * (1.0 + gam * theta21));
  Dual4 lambda3 = (z3 - dsqrt(z3 * z3 - 4.0)) * 0.5; // subsonic root (Eq. 12)

  double gm1_local = gm1;
  double gp1_local = gp1;
  double k7 = std::pow(gp1_local / 2.0, 1.0 / gm1_local);
  auto q_of = [&](const Dual4& lam) {
    Dual4 term = 1.0 - (gm1_local / gp1_local) * lam * lam;
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

EjectorCriticalPressureResult ejector_critical_back_pressure(
    double p_g, double t_g, double p_e, double t_e, const EjectorGeometry& geom,
    double gamma, double r_gas, double recovery_efficiency) {
  return ejector_critical_back_pressure_and_jacobian(p_g, t_g, p_e, t_e, geom, gamma, r_gas,
                                                      recovery_efficiency)
      .value;
}

} // namespace combaero::solver
