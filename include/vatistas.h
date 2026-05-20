#pragma once
#include <tuple>

// -------------------------------------------------------------------------
// Vatistas n-vortex model
// Reference: Vatistas, Kozel, Mih (1991), Exp. Fluids 11, 73-76.
//            DOI: 10.1007/BF00198434
//
// Models concentrated vortices with a smooth tangential-velocity profile:
//
//   V_theta(r) = [Gamma / (2*pi*r_c)] * V0_bar(r/r_c; n)
//   V0_bar(r_bar; n) = r_bar / (1 + r_bar^{2n})^{1/n}         [Eq. 2]
//
// The profile peaks at r = r_c for all n >= 1 and recovers the Rankine
// vortex in the limit n -> inf.  The shape parameter n = 2 gives the best
// fit to most experimental swirl data (Eq. 5 / Fig. 2a of the paper).
//
// All functions enforce n >= 1 (hard lower bound: singularities in the
// radial and axial velocities occur for n < 1, paper p.74).
// -------------------------------------------------------------------------

namespace vatistas {
constexpr double n_min = 1.0;     // lower bound on shape parameter
constexpr double n_default = 2.0; // best-fit integer (paper Fig. 2a)
} // namespace vatistas

namespace combaero::solver {

// ---- Normalised shape functions (dimensionless; r_bar = r / r_c) ---------

// V0_bar(r_bar; n) = r_bar / (1 + r_bar^{2n})^{1/n}
// Peak at r_bar = 1 with V0_bar_max = 2^{-1/n}.
double vatistas_v0_bar(double r_bar, double n);

// d(V0_bar)/d(r_bar) = (1 - r_bar^{2n}) / (1 + r_bar^{2n})^{1 + 1/n}
// Analytical; positive for r_bar < 1, zero at r_bar = 1, negative for r_bar > 1.
double vatistas_dv0_bar_drbar(double r_bar, double n);

// Normalised radial (inward) velocity: Vr_bar = V_r * r_c / nu
// Eq. (3):  Vr_bar = -2*(1+n)*r_bar^{2n-1} / (1 + r_bar^{2n})
// Negative by definition (inward flow toward vortex axis).
double vatistas_vr_bar(double r_bar, double n);

// d(Vr_bar)/d(r_bar)  -- analytical
double vatistas_dvr_bar_drbar(double r_bar, double n);

// Raw pressure integral: I(r_bar; n) = integral_0^{r_bar} V0_bar^2/r' dr'
// Dimensional delta_P = V_scale * I  where V_scale = (Gamma/(2*pi*r_c))^2 * rho.
//
// n = 1: r_bar^2 / (2*(1 + r_bar^2))     [closed form]
// n = 2: arctan(r_bar^2) / 2             [closed form, from Eq. 11 antiderivative]
// other: composite Simpson's rule (64 steps per domain segment)
double vatistas_pressure_integral(double r_bar, double n);

// d(I)/d(r_bar) = V0_bar(r_bar)^2 / r_bar  [analytical for all n, by FTC]
double vatistas_d_pressure_integral_drbar(double r_bar, double n);

// ---- Dimensional functions -----------------------------------------------

// Tangential velocity [m/s].  Throws if n < 1, Gamma <= 0, r_c <= 0.
double vatistas_v_theta(double r, double Gamma, double r_c,
                        double n = vatistas::n_default);

// (V_theta [m/s], dV_theta/dr [1/s], dV_theta/dGamma [(m/s)/(m^2/s)],
//  dV_theta/dr_c [(m/s)/m])
std::tuple<double, double, double, double>
vatistas_v_theta_and_jacobians(double r, double Gamma, double r_c,
                               double n = vatistas::n_default);

// Static pressure rise P(r) - P(0) [Pa].
double vatistas_delta_p(double r, double rho, double Gamma, double r_c,
                        double n = vatistas::n_default);

// (delta_P [Pa], d_delta_P/dr [Pa/m], d_delta_P/dGamma [Pa/(m^2/s)],
//  d_delta_P/dr_c [Pa/m])
std::tuple<double, double, double, double>
vatistas_delta_p_and_jacobians(double r, double rho, double Gamma, double r_c,
                               double n = vatistas::n_default);

} // namespace combaero::solver
