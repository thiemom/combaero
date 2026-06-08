#pragma once
#include "math_constants.h"
#include <cmath>

// Momentum-CV junction element constants and inline kernel helpers.
//
// Sanctioned successor to the K-closure tee (docs/junction/momentum cv
// implementation guide.pdf, supersedes the PLM in junction_model_v3.tex).
//
// Architecture (PDF Section 1):
//   Junction = pure conservation. Owns one scalar P_jct. Emits per-port
//   impulse-function residuals R_mom,i = (P_i + rho_i*u_i^2) - P_jct = 0 plus a
//   global mass residual sum_i mdot_i = 0. Sign-free (u_i^2 is direction-
//   invariant). Geometry-only (port angles + areas).
//
//   Loss = separate per-port elements (see BorderCarnotLossElement). The
//   junction itself carries no empirical content; turning / contraction /
//   expansion losses live on the loss elements bolted onto specific ports.
//
// Sign convention:
//   m_dot_i > 0 means flow OUT of the junction through port i (PDF Section 2.2).
//   The sum-mass residual sum_i mdot_i = 0 closes mass conservation.
//   u_i^2 in the impulse residual is sign-free, so the residual itself does
//   not care which way flow goes through any port.

namespace combaero {

// -----------------------------------------------------------------------------
// Constants
// -----------------------------------------------------------------------------

// Hager (1984) effective-angle correction for sharp-edged lateral branches:
// the local outflow angle gamma(v) varies along the branch; integrating along
// v in [0,1] gives an average outflow angle epsilon = delta_geom/4. The
// effective turning angle entering the momentum balance is
// theta_eff = delta_geom - epsilon = (3/4) * delta_geom. PDF Section 3.2.
inline constexpr double HAGER_FRACTION = 0.75;

// Per-port choking threshold. At M_i > MACH_CHOKE_THRESHOLD the port's
// impulse residual is swapped for R_choke,i = mdot_i - mdot_critical,i = 0
// (PDF Section 4, mirrors v3 MACH_CLAMP outer continuation).
inline constexpr double MACH_CHOKE_THRESHOLD = 0.95;

// Border-Carnot dynamic-head prefactor in L_i = 4 * (1 - cos(theta_eff))^2.
// PDF Section 3.1: the squared form (vs the naive 2*(1-cos(theta)) linear
// form) reproduces Hager xi_l and Bassett K_inc exactly at M -> 0.
inline constexpr double BC_LOSS_PREFACTOR = 4.0;

// -----------------------------------------------------------------------------
// Inline kernel helpers
// -----------------------------------------------------------------------------

// Border-Carnot turning-loss coefficient L = prefactor * (1 - cos(theta_eff))^2,
// with the Hager 3/4 correction applied to the geometric branch angle.
// Returns L (dimensionless), the coefficient in dP_loss = L * 0.5*rho*u^2.
// At delta_geom = 0 (straight-through port): L = 0 (no loss).
// At delta_geom = pi/2 (90-deg lateral): theta_eff = 3*pi/8 (67.5 deg);
//   L = 4 * (1 - cos(67.5 deg))^2 ~ 1.39.
inline double border_carnot_L(double delta_geom) {
    const double theta_eff = HAGER_FRACTION * delta_geom;
    const double one_minus_cos = 1.0 - std::cos(theta_eff);
    return BC_LOSS_PREFACTOR * one_minus_cos * one_minus_cos;
}

// Derivative of border_carnot_L w.r.t. delta_geom. delta_geom is geometric
// (not solved), so this is only needed for sensitivity diagnostics, not the
// solver Jacobian. Kept inline for symmetry with the tee_junction.h style.
inline double dborder_carnot_L_ddelta(double delta_geom) {
    const double theta_eff = HAGER_FRACTION * delta_geom;
    const double s = std::sin(theta_eff);
    const double one_minus_cos = 1.0 - std::cos(theta_eff);
    // d/d(delta) [4 * (1 - cos(0.75*delta))^2]
    //   = 4 * 2 * (1 - cos(0.75*delta)) * sin(0.75*delta) * 0.75
    return BC_LOSS_PREFACTOR * 2.0 * one_minus_cos * s * HAGER_FRACTION;
}

}  // namespace combaero
