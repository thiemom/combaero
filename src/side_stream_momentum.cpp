#include "side_stream_momentum.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "composition.h"
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
  // Only an entering flow accelerates out of the reservoir; a reversed one
  // discharges into it at static continuity (see the header).
  const double m_in = std::max(m, 0.0);
  StationHalfDrop out;
  out.dP = m_in * m_in * k;
  out.d_dm_a = 2.0 * m_in * k;  // the entering flow is "m_a"
  const double d_drho = -out.dP / rho;
  out.d_dP = d_drho * drho_dP;
  out.d_dT = d_drho * drho_dT;
  return out;
}

MergeFaceState chamber_merge_face_state(double m_main, double T0_main,
                                        const std::vector<double>& X_main,
                                        double m_out, double P, double T0,
                                        const std::vector<double>& X,
                                        double side_momentum, double area,
                                        bool incompressible) {
  if (!(area > 0.0) || !(P > 0.0) || !(T0 > 0.0) || !(T0_main > 0.0)) {
    throw std::invalid_argument(
        "chamber_merge_face_state: area, P, T and T_main must be positive");
  }

  if (incompressible) {
    // Low-Mach limit (#588): rho = P / (R T0) for both streams, at the
    // chamber pressure. Linear in P_f; every sensitivity is closed form.
    MergeFaceState out;
    const double R = specific_gas_constant(X);
    const double Rm = specific_gas_constant(X_main);
    const double A2 = area * area;
    const double a = m_main * m_main / A2;  // main impulse factor
    const double b = m_out * m_out / A2;    // outlet impulse factor
    const double vm = Rm * T0_main / P;     // 1 / rho_main
    const double vo = R * T0 / P;           // 1 / rho
    out.P_face = P + b * vo - side_momentum / area - a * vm;
    out.Pt_face = out.P_face + 0.5 * a * vm;
    out.M_face = std::abs(m_main) * vm / area / speed_of_sound(T0_main, X_main);
    out.choked = false;
    out.dPf_dm_main = -2.0 * m_main * vm / A2;
    out.dPf_dm_out = 2.0 * m_out * vo / A2;
    out.dPf_dJ = -1.0 / area;
    out.dPf_dP = 1.0 - b * vo / P + a * vm / P;
    out.dPf_dT = b * R / P;
    out.dPf_dT_main = -a * Rm / P;
    out.dPtf_dm_main = -m_main * vm / A2;
    out.dPtf_dm_out = out.dPf_dm_out;
    out.dPtf_dJ = out.dPf_dJ;
    out.dPtf_dP = 1.0 - b * vo / P + 0.5 * a * vm / P;
    out.dPtf_dT = out.dPf_dT;
    out.dPtf_dT_main = -0.5 * a * Rm / P;
    // Composition: R = R_u / MW in both specific volumes, dR/dX_j = -R W_j / MW.
    {
      const std::size_t n = X.size();
      const double MW = mwmix(X) / 1000.0;
      const double MWm = mwmix(X_main) / 1000.0;
      std::vector<double> gPf(n), gPfm(n), gPtm(n);
      for (std::size_t j = 0; j < n; ++j) {
        const double W = species_molar_mass(j) / 1000.0;
        gPf[j] = -b * vo * W / MW;
        gPfm[j] = a * vm * W / MWm;
        gPtm[j] = 0.5 * a * vm * W / MWm;
      }
      out.dPf_dY = mole_gradient_to_mass(X, gPf);
      out.dPtf_dY = out.dPf_dY;
      out.dPf_dY_main = mole_gradient_to_mass(X_main, gPfm);
      out.dPtf_dY_main = mole_gradient_to_mass(X_main, gPtm);
    }
    return out;
  }
  const double A = area;
  const double J = side_momentum;

  // Chamber (outlet) side: static T from the chamber's stagnation T0 and the
  // outlet mass flux (#578), the same solve as the chamber's own closure.
  const double R = specific_gas_constant(X);
  const double s_out = m_out >= 0.0 ? 1.0 : -1.0;
  const ChamberStaticState cs = chamber_static_state(T0, P, std::abs(m_out) / A, X);
  // Outlet momentum flux per area, m|m| / (rho A^2) = m|m| R Ts / (P A^2).
  const double q = m_out * std::abs(m_out) * R * cs.T / (P * A * A);
  const double Pi = P + q - J / A;
  const double dPi_dP = 1.0 + q * (cs.dT_dP / cs.T - 1.0 / P);
  const double dPi_dmout =
      2.0 * std::abs(m_out) * R * cs.T / (P * A * A) + q * (cs.dT_dG * s_out / A) / cs.T;
  const double dPi_dT0 = q * cs.dT_dT0 / cs.T;
  const double dPi_dJ = -1.0 / A;

  // Main face: static (P_f, T_f) from the impulse and the main stream's
  // stagnation T0_main (#578; T_main was read as static before):
  //   F1 = P_f + c T_f / P_f - Pi = 0,   c = m_main^2 R_m / A^2
  //   F2 = T_f - Ts(T0_main, P_f, |m_main| / A) = 0
  // Fixed point from T_f = T0_main: the quadratic in P_f on its subsonic
  // root, then the static T at that pressure. Contracts on the subsonic
  // branch (the face Mach is well below 1 wherever the quadratic has a root).
  const double Rm = specific_gas_constant(X_main);
  const double c = m_main * m_main * Rm / (A * A);
  const double Gm = std::abs(m_main) / A;
  double T_f = T0_main;
  double P_f = Pi;
  bool choked = false;
  ChamberStaticState fs{};
  for (int it = 0; it < 100; ++it) {
    const double D_raw = Pi * Pi - 4.0 * c * T_f;
    const double D_min = MERGE_DISCRIMINANT_FLOOR * Pi * Pi;
    choked = D_raw < D_min;
    P_f = 0.5 * (Pi + std::sqrt(std::max(D_raw, D_min)));
    fs = chamber_static_state(T0_main, P_f, Gm, X_main);
    const double dT = fs.T - T_f;
    T_f = fs.T;
    if (std::abs(dT) <= 1e-12 * T0_main) {
      break;
    }
  }
  MergeFaceState out;
  out.P_face = P_f;
  out.choked = choked;
  const double a_f = speed_of_sound(T_f, X_main);
  const double u_f = Gm * Rm * T_f / P_f;
  out.M_face = a_f > 0.0 ? u_f / a_f : 0.0;
  // Pt_f = P_f exp((s(T0_main) - s(T_f)) / R_m), consistent with the main
  // stream's own stagnation state; the entropy difference is the integrated
  // one from the static-state solve (no cancellation, #579).
  const double Pt_f = P_f * std::exp(fs.lnP0P);
  out.Pt_face = Pt_f;

  // Sensitivities of (P_f, T_f) by the implicit function theorem on
  // (F1, F2); with the discriminant floored, P_f follows Pi alone.
  // M = [[F1_Pf, F1_Tf], [F2_Pf, 1]].
  double dPf_dPi, dPf_dc, dPf_dT0m, dPf_dGm;
  if (choked) {
    dPf_dPi = 0.5 * (1.0 + std::sqrt(MERGE_DISCRIMINANT_FLOOR));
    dPf_dc = dPf_dT0m = dPf_dGm = 0.0;
  } else {
    const double F1_Pf = 1.0 - c * T_f / (P_f * P_f);
    const double F1_Tf = c / P_f;
    const double F2_Pf = -fs.dT_dP;
    const double det = F1_Pf - F1_Tf * F2_Pf;
    // d(P_f) = -(F1_x - F1_Tf * F2_x) / det for an input x with partials
    // (F1_x, F2_x); F2 has unit T_f coefficient.
    auto dPf = [&](double F1_x, double F2_x) { return -(F1_x - F1_Tf * F2_x) / det; };
    dPf_dPi = dPf(-1.0, 0.0);
    dPf_dc = dPf(T_f / P_f, 0.0);
    dPf_dT0m = dPf(0.0, -fs.dT_dT0);
    dPf_dGm = dPf(0.0, -fs.dT_dG);
  }
  // T_f = Ts(T0_main, P_f, Gm): its total derivative per input follows P_f.
  const double dTf_dPf = fs.dT_dP;
  const double cp_f = cp_mass(T_f, X_main);
  const double cp_0 = cp_mass(T0_main, X_main);
  // dPt_f = (Pt_f / P_f) dP_f + Pt_f dL, dL = -cp_f/(R T_f) dT_f + cp_0/(R T0) dT0.
  auto dPt = [&](double dPf_dx, double dTf_dx, double dT0m_dx) {
    return (Pt_f / P_f) * dPf_dx +
           Pt_f * (-cp_f / (Rm * T_f) * dTf_dx + cp_0 / (Rm * T0_main) * dT0m_dx);
  };

  const double sgn_m = m_main >= 0.0 ? 1.0 : -1.0;
  const double dc_dmm = 2.0 * m_main * Rm / (A * A);
  const double dGm_dmm = sgn_m / A;

  // Inputs that reach the face only through Pi.
  auto through_Pi = [&](double dPi_dx, double& dPf_out, double& dPt_out) {
    dPf_out = dPf_dPi * dPi_dx;
    dPt_out = dPt(dPf_out, dTf_dPf * dPf_out, 0.0);
  };
  through_Pi(dPi_dP, out.dPf_dP, out.dPtf_dP);
  through_Pi(dPi_dmout, out.dPf_dm_out, out.dPtf_dm_out);
  through_Pi(dPi_dJ, out.dPf_dJ, out.dPtf_dJ);
  through_Pi(dPi_dT0, out.dPf_dT, out.dPtf_dT);

  // m_main: through c and through the face mass flux Gm.
  out.dPf_dm_main = dPf_dc * dc_dmm + dPf_dGm * dGm_dmm;
  {
    const double dTf = dTf_dPf * out.dPf_dm_main + fs.dT_dG * dGm_dmm;
    out.dPtf_dm_main = dPt(out.dPf_dm_main, dTf, 0.0);
  }
  // T0_main: through T_f (F2) and through s(T0_main) in Pt_f.
  out.dPf_dT_main = dPf_dT0m;
  {
    const double dTf = dTf_dPf * dPf_dT0m + fs.dT_dT0;
    out.dPtf_dT_main = dPt(dPf_dT0m, dTf, 1.0);
  }

  // Composition (#743). Outlet gas X: through Pi (R and the outlet static T
  // in q). Main X: through c = m_main^2 R_m / A^2 (F1), the face static T
  // (F2) and the entropy difference in Pt_f.
  {
    const std::size_t n = X.size();
    const double MW = mwmix(X) / 1000.0;
    const double MWm = mwmix(X_main) / 1000.0;
    const ChamberStaticStateDX cx = chamber_static_state_dX(T0, P, std::abs(m_out) / A, X);
    const ChamberStaticStateDX fx = chamber_static_state_dX(T0_main, P_f, Gm, X_main);
    std::vector<double> gPf(n), gPt(n), gPfm(n), gPtm(n);
    for (std::size_t j = 0; j < n; ++j) {
      const double W = species_molar_mass(j) / 1000.0;
      through_Pi(q * (-W / MW + cx.dT_dX[j] / cs.T), gPf[j], gPt[j]);
      const double dc = -c * W / MWm;
      double dPf_x = 0.0;
      if (!choked) {
        const double F1_Pf = 1.0 - c * T_f / (P_f * P_f);
        const double F1_Tf = c / P_f;
        const double det = F1_Pf - F1_Tf * (-fs.dT_dP);
        // F1_x = (T_f / P_f) dc, F2_x = -dTs/dX at fixed P_f.
        dPf_x = -((T_f / P_f) * dc - F1_Tf * (-fx.dT_dX[j])) / det;
      }
      const double dTf = dTf_dPf * dPf_x + fx.dT_dX[j];
      gPfm[j] = dPf_x;
      // ln(P0/P) at the face moves with X directly and through T_f.
      const double dL = fx.dlnP0P_dX[j] - fx.dT_dX[j] * (-cp_f / (Rm * T_f)) +
                        (-cp_f / (Rm * T_f)) * dTf;
      gPtm[j] = (Pt_f / P_f) * dPf_x + Pt_f * dL;
    }
    out.dPf_dY = mole_gradient_to_mass(X, gPf);
    out.dPtf_dY = mole_gradient_to_mass(X, gPt);
    out.dPf_dY_main = mole_gradient_to_mass(X_main, gPfm);
    out.dPtf_dY_main = mole_gradient_to_mass(X_main, gPtm);
  }
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
