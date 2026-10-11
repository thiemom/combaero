#include "side_stream_momentum.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "compressible.h"
#include "solver_interface.h"
#include "stagnation.h"
#include "thermo.h"

namespace combaero {

namespace {
// Floor on the impulse quadratic's discriminant, relative to Pi^2, so a
// choked face stays finite and differentiable (flagged, not hidden).
constexpr double MERGE_DISCRIMINANT_FLOOR = 1e-8;
// Central-difference steps for the smooth temperature sensitivities.
constexpr double T_STEP_REL = 1e-6;
constexpr double P_STEP_REL = 1e-7;
}  // namespace

StationHalfDrop station_half_drop(double m_a, double m_b, double P, double T,
                                  const std::vector<double>& X, double area,
                                  double kappa) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument("station_half_drop: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double k = 0.5 / (rho * area * area);
  const double f = m_b * std::abs(m_b) - m_a * std::abs(m_a) + kappa * m_a * (m_a - m_b);
  StationHalfDrop out;
  out.dP = f * k;
  out.d_dm_a = (-2.0 * std::abs(m_a) + kappa * (2.0 * m_a - m_b)) * k;
  out.d_dm_b = (2.0 * std::abs(m_b) - kappa * m_a) * k;
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  return out;
}

StationHalfDrop channel_entry_drop(double m, double P, double T,
                                   const std::vector<double>& X, double area,
                                   double K_in) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument("channel_entry_drop: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double k = 0.5 * (1.0 + K_in) / (rho * area * area);
  StationHalfDrop out;
  out.dP = m * std::abs(m) * k;
  out.d_dm_a = 2.0 * std::abs(m) * k;  // the entering flow is "m_a"
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  return out;
}

namespace {
struct FaceCore {
  double P_f, Pt_f, M_f;
  bool choked;
  // analytic partials of P_f wrt Pi and c, and of Pt_f wrt P_f and m_main
  double dPf_dPi, dPf_dc, dPt_dPf, dPt_dm;
};

FaceCore face_core(double m_main, double T_main, const std::vector<double>& X_main,
                   double m_out, double P, double T, const std::vector<double>& X,
                   double J, double A, double* dPi_dP_out, double* dPi_dmout_out,
                   double* c_out) {
  const double rho = density(T, P, X);
  const double Pi = P + m_out * std::abs(m_out) / (rho * A * A) - J / A;
  const double RT = specific_gas_constant(X_main) * T_main;
  const double c = m_main * m_main * RT / (A * A);
  const double D_raw = Pi * Pi - 4.0 * c;
  const double D_min = MERGE_DISCRIMINANT_FLOOR * Pi * Pi;
  const bool choked = D_raw < D_min;
  const double sD = std::sqrt(std::max(D_raw, D_min));
  FaceCore f{};
  f.choked = choked;
  f.P_f = 0.5 * (Pi + sD);
  // Held at the floor, the root no longer moves with c; with Pi it moves as
  // the floored discriminant does.
  f.dPf_dPi = choked ? 0.5 * (1.0 + std::sqrt(MERGE_DISCRIMINANT_FLOOR))
                     : 0.5 * (1.0 + Pi / sD);
  f.dPf_dc = choked ? 0.0 : -1.0 / sD;
  const double a = speed_of_sound(T_main, X_main);
  const double u = std::abs(m_main) * RT / (f.P_f * A);
  f.M_f = a > 0.0 ? u / a : 0.0;
  auto [P0, dP0_dM] = solver::P0_from_static_and_jacobian_M(f.P_f, T_main, f.M_f, X_main);
  f.Pt_f = P0;
  // P0 = P exp(...) at fixed (M, T): dP0/dP_f = P0/P_f, and M ~ 1/P_f.
  f.dPt_dPf = P0 / f.P_f + dP0_dM * (-f.M_f / f.P_f);
  const double sgn = m_main >= 0.0 ? 1.0 : -1.0;
  f.dPt_dm = (a > 0.0 && f.P_f > 0.0) ? dP0_dM * sgn * RT / (f.P_f * A * a) : 0.0;
  // Ideal gas: drho/dP = rho/P at fixed T.
  if (dPi_dP_out) *dPi_dP_out = 1.0 - m_out * std::abs(m_out) / (rho * A * A) / P;
  if (dPi_dmout_out) *dPi_dmout_out = 2.0 * std::abs(m_out) / (rho * A * A);
  if (c_out) *c_out = c;
  return f;
}
}  // namespace

MergeFaceState chamber_merge_face_state(double m_main, double T_main,
                                        const std::vector<double>& X_main,
                                        double m_out, double P, double T,
                                        const std::vector<double>& X,
                                        double side_momentum, double area) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0) || !(T_main > 0.0)) {
    throw std::invalid_argument(
        "chamber_merge_face_state: area, P, T and T_main must be positive");
  }
  double dPi_dP = 0.0, dPi_dmout = 0.0, c = 0.0;
  const FaceCore f = face_core(m_main, T_main, X_main, m_out, P, T, X,
                               side_momentum, area, &dPi_dP, &dPi_dmout, &c);
  MergeFaceState out;
  out.P_face = f.P_f;
  out.Pt_face = f.Pt_f;
  out.M_face = f.M_f;
  out.choked = f.choked;

  // P_f(Pi, c): Pi carries P, m_out and J; c carries m_main.
  const double dc_dmm = 2.0 * c / (m_main != 0.0 ? m_main : 1.0);
  out.dPf_dP = f.dPf_dPi * dPi_dP;
  out.dPf_dm_out = f.dPf_dPi * dPi_dmout;
  out.dPf_dJ = f.dPf_dPi * (-1.0 / area);
  out.dPf_dm_main = m_main != 0.0 ? f.dPf_dc * dc_dmm : 0.0;
  out.dPtf_dP = f.dPt_dPf * out.dPf_dP;
  out.dPtf_dm_out = f.dPt_dPf * out.dPf_dm_out;
  out.dPtf_dJ = f.dPt_dPf * out.dPf_dJ;
  out.dPtf_dm_main = f.dPt_dPf * out.dPf_dm_main + f.dPt_dm;

  // Temperatures: central differences of the whole (smooth) evaluation.
  auto eval = [&](double Tm, double Tc) {
    return face_core(m_main, Tm, X_main, m_out, P, Tc, X, side_momentum, area,
                     nullptr, nullptr, nullptr);
  };
  const double hm = std::max(1e-4, T_main * T_STEP_REL);
  const FaceCore mp = eval(T_main + hm, T), mm = eval(T_main - hm, T);
  out.dPf_dT_main = (mp.P_f - mm.P_f) / (2.0 * hm);
  out.dPtf_dT_main = (mp.Pt_f - mm.Pt_f) / (2.0 * hm);
  const double hc = std::max(1e-4, T * T_STEP_REL);
  const FaceCore cp = eval(T_main, T + hc), cm = eval(T_main, T - hc);
  out.dPf_dT = (cp.P_f - cm.P_f) / (2.0 * hc);
  out.dPtf_dT = (cp.Pt_f - cm.Pt_f) / (2.0 * hc);
  return out;
}

namespace {
struct JetW {
  double w;
  bool choked;
};

JetW jet_w(double Pt, double Tt, double P, const std::vector<double>& X) {
  if (!(Pt > 0.0) || !(Tt > 0.0) || !(P > 0.0) || P >= Pt) {
    return {0.0, false};
  }
  const double P_star = critical_pressure_ratio(Tt, Pt, X) * Pt;
  if (P >= P_star) {
    const double M = mach_from_pressure_ratio(Tt, Pt, P, X);
    const double Ts = T_from_stagnation(Tt, M, X);
    return {M * speed_of_sound(Ts, X), false};
  }
  const double T_star = T_from_stagnation(Tt, 1.0, X);
  const double u_star = speed_of_sound(T_star, X);
  const double rho_star = density(T_star, P_star, X);
  return {u_star + (P_star - P) / (rho_star * u_star), true};
}
}  // namespace

JetImpulse jet_impulse(double m, double Pt, double Tt, double P,
                       const std::vector<double>& X) {
  const JetW base = jet_w(Pt, Tt, P, X);
  JetImpulse out;
  out.J = m * base.w;
  out.dJ_dm = base.w;
  out.choked = base.choked;
  if (base.w == 0.0) {
    return out;
  }
  // w(Pt, Tt, P) is smooth (C1 across choking); central differences, kept
  // inside C++ like the chamber's own stagnation sensitivities.
  const double hP = std::max(1.0, Pt * P_STEP_REL);
  const double hT = std::max(1e-4, Tt * T_STEP_REL);
  out.dJ_dPt = m * (jet_w(Pt + hP, Tt, P, X).w - jet_w(Pt - hP, Tt, P, X).w) / (2.0 * hP);
  out.dJ_dTt = m * (jet_w(Pt, Tt + hT, P, X).w - jet_w(Pt, Tt - hT, P, X).w) / (2.0 * hT);
  out.dJ_dP = m * (jet_w(Pt, Tt, P + hP, X).w - jet_w(Pt, Tt, P - hP, X).w) / (2.0 * hP);
  return out;
}

CrossflowRatio crossflow_velocity_ratio(double m_a, double m_b, double P,
                                        double T, const std::vector<double>& X,
                                        double P_down, double area) {
  if (!(area > 0.0) || !(P > 0.0) || !(T > 0.0)) {
    throw std::invalid_argument(
        "crossflow_velocity_ratio: area, P and T must be positive");
  }
  auto [rho, drho_dT, drho_dP] = solver::density_and_jacobians(T, P, X);
  const double m_mean = 0.5 * (m_a + m_b);
  const double sgn = m_mean >= 0.0 ? 1.0 : -1.0;
  const double U1 = std::abs(m_mean) / (rho * area);
  // Unit-flow jet impulse: J = w, the isentropic velocity from (P, T) to P_down.
  const JetImpulse jet = jet_impulse(1.0, P, T, P_down, X);
  const double Vi = std::sqrt(jet.J * jet.J + CROSSFLOW_VI_FLOOR * CROSSFLOW_VI_FLOOR);
  const double dVi_dw = Vi > 0.0 ? jet.J / Vi : 0.0;

  CrossflowRatio out;
  out.U1 = U1;
  out.Vi = Vi;
  out.U1_over_Vi = U1 / Vi;
  const double dU1_dm = 0.5 * sgn / (rho * area);
  out.d_dm_a = dU1_dm / Vi;
  out.d_dm_b = dU1_dm / Vi;
  const double dU1_drho = -U1 / rho;
  const double r_over_Vi = out.U1_over_Vi / Vi;
  out.d_dP = dU1_drho * drho_dP / Vi - r_over_Vi * dVi_dw * jet.dJ_dPt;
  out.d_dT = dU1_drho * drho_dT / Vi - r_over_Vi * dVi_dw * jet.dJ_dTt;
  out.d_dP_down = -r_over_Vi * dVi_dw * jet.dJ_dP;
  return out;
}

}  // namespace combaero
