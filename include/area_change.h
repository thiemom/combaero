#ifndef AREA_CHANGE_H
#define AREA_CHANGE_H

// -----------------------------------------------------------------------------
// Mach number thresholds (shared by all area-change elements)
// -----------------------------------------------------------------------------
//
// Two named levels separate physics validity from numerical safety:
//
//   MACH_ACCURACY_LIMIT  -- upper bound of the calibrated regime for the
//     compressibility correction  k1 = 1 + 0.35 * M^2.  Results above this
//     value are physically extrapolated; the correction continues to grow but
//     is no longer anchored to empirical data.  Use this threshold in
//     post-processing diagnostics and convergence reports.
//
//   MACH_CLAMP  -- hard numerical ceiling applied inside every element.
//     k1 is evaluated at min(|Mach|, MACH_CLAMP), so its value is bounded
//     regardless of what the network solver presents during iteration.  This
//     prevents unbounded k1 growth from derailing Newton steps when a
//     non-converged iterate overshoots into unphysical supersonic territory.
//     The clamp is intentionally set above M=1 so that the element does not
//     create a spurious fixed point at the sonic line.
//
//   AreaChangeResult::mach_clamped is set true whenever |Mach| > MACH_CLAMP.
//   It is always false for physically subsonic inputs and can be aggregated
//   across elements to flag branches that require a compressible-flow model.
// -----------------------------------------------------------------------------

namespace combaero {

constexpr double MACH_ACCURACY_LIMIT = 0.9;  // k1 correlation valid below here
constexpr double MACH_CLAMP          = 1.5;  // hard cap for Newton robustness

// Result struct returned by all area-change elements.
struct AreaChangeResult {
    double dP;            // Total pressure drop [Pa], signed
    double dS_dm;         // d(dP)/d(m_dot)  [Pa*s/kg]
    double dS_drho;       // d(dP)/d(rho)    [Pa*m^3/kg]
    double dS_dmu;        // d(dP)/d(mu)     [Pa*s/(Pa*s)] = [-]
    bool   mach_clamped;  // true if |Mach| was reduced to MACH_CLAMP
};

// -----------------------------------------------------------------------------
// Sharp-Edge Area Change Element (Expansion / Contraction)
// -----------------------------------------------------------------------------
//
// Computes total pressure drop and analytical Jacobian partials for a sudden
// sharp-edge area change (expansion or contraction).  Designed for use with
// scipy.optimize.root (hybr/lm).  Smooth behavior at zero flow is required
// for Newton convergence.
//
// Physics:
//   Two distinct Idelchik/Weisbach correlations selected by geometry,
//   smoothly blended by flow direction via sigmoid:
//
//   Expansion (Borda-Carnot), referenced to F_small:
//     zeta_exp = alpha - 2*ar + ar^2
//     ar       = F_small / F_large,  always in (0, 1]
//     alpha    = 1 + 2.73 / (Re + 10)^0.4
//
//   Contraction (Weisbach sharp-edge), referenced to F_small:
//     zeta_con = 0.5 * (1 - ar)
//
//   Direction blend (sigmoid):
//     w        = sigma(m_dot / m_scale) = 1 / (1 + exp(-m_dot / m_scale))
//     zeta_eff = w * zeta_fwd + (1 - w) * zeta_rev
//
//   Pressure drop:
//     dP = k1 * zeta_eff * m_dot * |m_dot| / (2 * rho * F_small^2)
//     k1 = 1 + 0.35 * M^2,  M = min(|Mach|, MACH_CLAMP)
//     Physically calibrated for M < MACH_ACCURACY_LIMIT.
//
// What this model does NOT cover:
//   - Gradual area change (diffuser/nozzle): needs half-angle input
//   - Rounded or beveled edges: separate Cc correction
//   - Two-phase flow
//   - Re < ~100 laminar contraction: Weisbach zeta_con not validated
//
// References:
//   - Idelchik, I.E. (1960) Handbook of Hydraulic Resistance (various editions)
//   - Weisbach, J. (1855) Die Neue Markscheidekunst
// -----------------------------------------------------------------------------

// Parameters:
//   m_dot   : mass flow rate [kg/s], signed (positive = F0 -> F1)
//   rho     : density [kg/m^3]
//   mu      : dynamic viscosity [Pa*s]
//   F0      : area at node-0 side [m^2]
//   F1      : area at node-1 side [m^2]
//   Mach    : Mach number (default 0); clamped to MACH_CLAMP internally
//   m_scale : sigmoid transition width [kg/s], tune to ~1% of nominal flow
//   D_h     : hydraulic diameter [m] for Re calculation (default 0 = circular,
//             i.e. D_h = sqrt(4*F_small/pi)).  Set explicitly for non-circular
//             cross-sections where D_h and A are independent.
AreaChangeResult sharp_area_change(double m_dot, double rho, double mu,
                                   double F0, double F1,
                                   double Mach = 0.0,
                                   double m_scale = 1e-4,
                                   double D_h = 0.0);

// -----------------------------------------------------------------------------
// Conical (Smooth) Area Change Element (Diffuser / Nozzle)
// -----------------------------------------------------------------------------
//
// Computes total pressure drop and analytical Jacobian partials for a gradual
// conical area change.  Companion to sharp_area_change; shares the same return
// struct and sign convention.
//
// Geometry:
//   Half-angle theta is derived from F0, F1, and the axial length L:
//     r_small = sqrt(F_small / pi)
//     r_large = sqrt(F_large / pi)
//     theta   = atan2(r_large - r_small, L)   [rad, always >= 0]
//   Assumes circular cross-sections.  For non-circular ducts pass the
//   hydraulic-radius-equivalent area: F_eq = pi * (D_h / 2)^2.
//
// Loss coefficients (both purely geometric — no Re dependence):
//
//   Expansion (diffuser), referenced to F_small:
//     zeta_exp = eta(theta) * (1 - ar)^2
//     eta(theta) = 1 - exp(-k_eta * theta^n_eta)
//     k_eta = 3.2,  n_eta = 0.8
//     Limits: theta=0 -> eta=0 (ideal pipe); theta=pi/2 -> eta->1 (Borda-Carnot)
//
//   Contraction (nozzle), referenced to F_small:
//     zeta_con = C_con(theta) * (1 - ar)
//     C_con = 0.04 + 0.46 * sin(theta)^1.5
//     Limits: theta=0  -> C_con=0.04 (smooth, minimal vena contracta)
//             theta=pi/2 -> C_con=0.50 (recovers Weisbach sharp-edge)
//
//   Direction blend and pressure drop: identical to sharp_area_change.
//
// Consistency with sharp_area_change:
//   At theta -> pi/2 (length -> 0):
//     zeta_exp -> (1 - ar)^2          (Borda-Carnot, alpha~=1 at high Re)
//     zeta_con -> 0.5 * (1 - ar)      (Weisbach sharp-edge)
//
// What this model does NOT cover:
//   - Distributed friction (f*L/D): add Darcy-Weisbach pipe in series for L/D > 5
//   - Non-circular cross-sections: hydraulic diameter equivalence assumed
//   - Bell-mouth / curved profiles: Idelchik diagram 5-21, separate element
//   - Re < ~500: laminar diffuser separation physics differ substantially
//   - Two-phase flow
//
// References:
//   - Idelchik, I.E. (1960) Handbook of Hydraulic Resistance, diagrams 5-2, 5-23
// -----------------------------------------------------------------------------

// Parameters:
//   m_dot   : mass flow rate [kg/s], signed (positive = F0 -> F1)
//   rho     : density [kg/m^3]
//   mu      : dynamic viscosity [Pa*s]  (accepted for API uniformity; unused)
//   F0      : area at node-0 side [m^2]
//   F1      : area at node-1 side [m^2]
//   length  : axial length of the conical section [m]
//   Mach    : Mach number (default 0); clamped to MACH_CLAMP internally
//   m_scale : sigmoid transition width [kg/s], tune to ~1% of nominal flow
//
// Note: loss coefficients are purely geometric (no Re dependence).  For long
//   sections (L/D_small > 5) add a Darcy-Weisbach pipe element in series.
//   dS_dmu is always 0.
AreaChangeResult conical_area_change(double m_dot, double rho, double mu,
                                     double F0, double F1,
                                     double length,
                                     double Mach = 0.0,
                                     double m_scale = 1e-4);

} // namespace combaero

#endif // AREA_CHANGE_H
