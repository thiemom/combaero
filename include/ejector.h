#pragma once

// -----------------------------------------------------------------------------
// 1-D critical-mode supersonic ejector model.
//
// C++ port of python/combaero/network/_ejector_huang1999.py -- see that
// module's docstring for the full paper-by-paper research narrative and
// validation/ejector/README.md for the accuracy comparison against
// alternative closures. Kept in exact correspondence with the Python
// reference so the two can be cross-validated row-for-row against
// validation/ejector/data/huang1999_reference_data.h (see
// tests/test_ejector_physics.cpp, tests/test_ejector_jacobians.cpp).
//
// References
// ----------
// [Huang 1999] B.J. Huang, J.M. Chang, C.P. Wang, V.A. Petrenko,
//     "A 1-D analysis of ejector performance",
//     International Journal of Refrigeration 22 (1999) 354-364.
//     Local copy: docs/ejector/Huang_1d_analysis_ejector.pdf
//     Source of ejector_entrainment_ratio (Eqs. 1-8) and
//     ejector_choked_mass_flow (Eqs. 1, 7).
// [Kracik & Dvorak 2016] J. Kracik, V. Dvorak,
//     "Development of an Analytical Method for Predicting Flow in a
//     Supersonic Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
//     DOI: 10.1051/epjconf/201611402059.
//     Local copy: docs/ejector/Development_of_an_Analytical_Method_for.pdf
//     Source of the mixing closure ejector_critical_back_pressure
//     implements (their Eqs. 7-13).
//
// Scope. The critical-mode scalar closures (ejector_entrainment_ratio,
// ejector_critical_back_pressure and their Jacobians) model double-choked
// operation only; the additional operating-regime closures below
// (ejector_cd_nozzle_mass_flow, ejector_jetpump_discharge, ...) extend to the
// subcritical-droop and unchoked-primary jet-pump regimes, and
// ejector_element_residuals_and_jacobian composes ALL of them into the
// network element's full 4-row residual system spanning every regime. Design
// + provenance: validation/ejector/OPERATING_REGIMES_DESIGN.md.
//
// Jacobian scope. The critical-mode scalar closures expose derivatives w.r.t.
// the 4 thermodynamic inputs (p_g, t_g, p_e, t_e) only -- the Newton-solved
// unknowns for the choked-plateau path; gamma, r_gas, geometry, and
// recovery_efficiency are frozen element parameters (deliberate scope
// boundary, not an oversight -- validation/ejector/data/
// huang1999_reference_data.h separately carries central-difference targets
// for those too, for a possible future design-sensitivity extension). The
// operating-regime closures additionally carry the P_py / omega / m_dot
// partials the extended element needs, and the whole-element function returns
// the full Jacobian w.r.t. all nine element unknowns.
//
// Key simplification that makes every critical-mode derivative below exact
// analytic chain rule (no implicit/root-finding differentiation needed there):
// the one internal root-find (ejector_mach_from_area_ratio_supersonic, a
// bisection for the primary nozzle-exit Mach number M_p1) depends ONLY on
// area_ratio_nozzle and gamma -- neither a Newton-solved unknown here -- so
// M_p1 has a PROVABLY exact zero derivative w.r.t. p_g, t_g, p_e, t_e.  It
// is computed once (bisection, same as Python) and then treated as a
// precomputed constant throughout the rest of the (otherwise fully
// closed-form) chain.
// -----------------------------------------------------------------------------

#include <array>
#include <tuple>

namespace ejector {
constexpr double eta_p = 0.95; // primary-nozzle isentropic efficiency (Huang 1999, Section 4)
constexpr double eta_s = 0.85; // entrained-flow isentropic efficiency
constexpr double phi_p = 0.88; // primary-core area-loss coefficient at section y-y (Eq. 5)
} // namespace ejector

namespace combaero::solver {

struct EjectorGeometry {
  double area_ratio_nozzle; // A_p1 / A_t
  double area_ratio_mix;    // A_3 / A_t
};

// -----------------------------------------------------------------------------
// Value-only functions (Python 1:1 parity; see tests/test_ejector_physics.cpp)
// -----------------------------------------------------------------------------

// Isentropic area ratio A/A* for a given Mach number.
double ejector_area_over_astar(double mach, double gamma);

// Supersonic root of A/A* = area_ratio by bisection (inverse of Eq. 2).
// A/A* is monotonically increasing for M > 1, so bisection is
// unconditionally convergent.
double ejector_mach_from_area_ratio_supersonic(double area_ratio, double gamma);

// Choked (sonic-throat) mass flow rate (Huang 1999 Eqs. 1 and 7).
//
// Method: Analytical Chain Rule (see ejector_choked_mass_flow_and_jacobian).
double ejector_choked_mass_flow(double p0, double t0, double area_throat,
                                double gamma, double r_gas, double eta);

// Converging-diverging nozzle mass flow spanning choked and unchoked, the R0
// primary residual of the operating-regime extension (see
// validation/ejector/OPERATING_REGIMES_DESIGN.md sec 3.1). Subsonic exit flux
// (area_exit, exit static p_static_down) smooth-min'd with the sonic-throat cap
// (area_throat); the two are equal at the choke threshold so the min is
// continuous, and the C1 sqrt-smoothing (eps = eps_frac * cap) rounds the
// corner. 1:1 with cd_nozzle_mass_flow in _ejector_huang1999.py.
double ejector_cd_nozzle_mass_flow(double p0, double t0, double p_static_down,
                                   double area_throat, double area_exit,
                                   double gamma, double r_gas, double eta,
                                   double eps_frac);

struct EjectorEntrainmentResult {
  double omega;                     // entrainment ratio ms / mp
  double mach_nozzle_exit;          // M_p1
  double mach_hypothetical_throat;  // M_py
  double area_ratio_primary_core;   // A_py / A_t
  double area_ratio_entrained;      // A_sy / A_t
  double p_mixing_pa;               // P_py = P_sy
};

// Critical-mode entrainment ratio omega (Huang 1999 Eqs. 1-8).
//
// Method: Analytical Chain Rule, composed around one precomputed bisection
// root (M_p1) whose derivative w.r.t. p_g/t_g/p_e/t_e is exactly zero --
// see the header comment above.
EjectorEntrainmentResult ejector_entrainment_ratio(double p_g, double t_g,
                                                    double p_e, double t_e,
                                                    const EjectorGeometry& geom,
                                                    double gamma);

// Aerodynamic mass-flow function q(lambda) (Kracik & Dvorak Eq. 9).
double ejector_q_lambda(double lam, double gamma);

struct EjectorCriticalPressureResult {
  double p_c_pa;                  // critical back pressure P_c*
  double p_mixed_stagnation_pa;   // p03, lossless (recovery_efficiency=1) mixed stagnation pressure
  double temp_mixed_stagnation_k; // T03
  double lambda_mixed;            // lambda3, subsonic root
  double mach_mixed;              // Mach number equivalent of lambda3
};

// Critical back pressure P_c* via Kracik & Dvorak's mixing closure
// (their Eqs. 7-13), continuing from the entrainment-ratio state above.
//
// Method: Analytical Chain Rule (same M_p1 treatment as
// ejector_entrainment_ratio).
EjectorCriticalPressureResult ejector_critical_back_pressure(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma, double r_gas,
    double recovery_efficiency);

// -----------------------------------------------------------------------------
// Jacobian-carrying variants: (value, d/dp_g, d/dt_g, d/dp_e, d/dt_e).
// -----------------------------------------------------------------------------

// (m_dot, d_mdot_dp0, d_mdot_dt0)
std::tuple<double, double, double>
ejector_choked_mass_flow_and_jacobian(double p0, double t0, double area_throat,
                                      double gamma, double r_gas, double eta);

struct EjectorCDNozzleJacobian {
  double mdot;        // == ejector_cd_nozzle_mass_flow
  double dmdot_dp0;   // d/d(primary stagnation pressure)
  double dmdot_dt0;   // d/d(primary stagnation temperature)
  double dmdot_dp_py; // d/d(exit/mixing static pressure P_py)
};

// C-D nozzle mass flow with its analytic Jacobian w.r.t. the three inputs that
// are Newton-solved unknowns for EjectorElement (p0 = primary Pt, t0 = primary
// Tt, p_static_down = P_py). gamma/r_gas/geometry/eta are frozen parameters.
// The secondary entrained flux reuses this with area_throat = area_exit = A_s.
EjectorCDNozzleJacobian ejector_cd_nozzle_mass_flow_and_jacobian(
    double p0, double t0, double p_static_down, double area_throat,
    double area_exit, double gamma, double r_gas, double eta, double eps_frac);

struct EjectorCDExitStaticJacobian {
  double p_py;         // exit static that passes m_dot (== cd_nozzle_exit_static)
  double dp_py_dm_dot; // d/d(mass flow)
  double dp_py_dp0;    // d/d(primary stagnation pressure)
  double dp_py_dt0;    // d/d(primary stagnation temperature)
};

// Inverse of the C-D nozzle: the exit static P_py that passes m_dot, with its
// implicit-function-theorem Jacobian. The derived mixing pressure of the
// operating-regime element's unchoked branch (OPERATING_REGIMES_DESIGN.md sec
// 6b); 1:1 with cd_nozzle_exit_static in _ejector_huang1999.py. The partials
// are ratios of the C-D nozzle's own Jacobian at the solved P_py:
//   dP_py/dm_dot = 1/(dmdot/dP_py),  dP_py/dp0 = -(dmdot/dp0)/(dmdot/dP_py), ...
EjectorCDExitStaticJacobian ejector_cd_nozzle_exit_static_and_jacobian(
    double m_dot, double p0, double t0, double area_throat, double area_exit,
    double gamma, double r_gas, double eta, double eps_frac);

struct EjectorJetPumpDischargeJacobian {
  double p03;      // recovery_efficiency * mixed stagnation pressure
  double dp03_dp_g;
  double dp03_dt_g;
  double dp03_dp_e;
  double dp03_dt_e;
  double dp03_dp_py; // d/d(mixing static pressure P_py)
  double dp03_domega; // d/d(entrainment ratio -- the R1 unknown feeding the mix)
};

// Jet-pump mixed-flow discharge stagnation: the Kracik & Dvorak mixing closure
// (reused from the critical path, ejector_critical_back_pressure) evaluated at
// the SUBSONIC-primary jet-pump state -- both streams expand to the common
// mixing static P_py (lambda1/lambda2 from that expansion), mixed at ratio
// omega. The R3 (outlet-pin) building block for the operating-regime element;
// 1:1 with _mixed_flow_stagnation in _ejector_huang1999.py fed jet-pump
// lambdas. r_gas cancels (lambda/q are dimensionless), so it is not an input.
EjectorJetPumpDischargeJacobian ejector_jetpump_discharge_and_jacobian(
    double p_g, double t_g, double p_e, double t_e, double p_py, double omega,
    double gamma, double recovery_efficiency);

struct EjectorEntrainmentJacobian {
  EjectorEntrainmentResult value;
  double domega_dp_g;
  double domega_dt_g;
  double domega_dp_e;
  double domega_dt_e;
};

EjectorEntrainmentJacobian ejector_entrainment_ratio_and_jacobian(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma);

struct EjectorCriticalPressureJacobian {
  EjectorCriticalPressureResult value;
  double dpc_dp_g;
  double dpc_dt_g;
  double dpc_dp_e;
  double dpc_dt_e;
};

EjectorCriticalPressureJacobian ejector_critical_back_pressure_and_jacobian(
    double p_g, double t_g, double p_e, double t_e,
    const EjectorGeometry& geom, double gamma, double r_gas,
    double recovery_efficiency);

// -----------------------------------------------------------------------------
// Whole-element (f, J): the 4-row operating-regime residual system.
//
// Composes the scalar closures above (choked/C-D nozzle mass flow, critical
// entrainment + back pressure, jet-pump discharge) into the four coupled rows
// the network EjectorElement solves, blending the critical/subcritical/unchoked
// regimes by two smootherstep weights, and returns all four residuals with
// their analytic Jacobian w.r.t. the nine Newton unknowns -- forward-mode dual
// (DualN<9>) through the whole blend. This is the C++ home of the assembly
// (matching the whole-element (f, J) practice of MultiPortChamberElement /
// TeeJunctionElement); the Python EjectorElement is a thin relabeling shim.
// 1:1 with EjectorElement.residuals() in ejector_element.py. Full derivation
// and provenance: validation/ejector/OPERATING_REGIMES_DESIGN.md sec 6c/8.3.
//
// Unknown (seed) order -- both `jac` columns and the physical-flow convention:
//   0 mp        primary mass flow (physical, = -port_mdot_primary)
//   1 ms        secondary (entrained) mass flow (physical)
//   2 mdot_out  outlet mass flow (physical)
//   3 p_g       primary stagnation pressure   Pt_primary
//   4 t_g       primary stagnation temperature Tt_primary
//   5 p_e       secondary stagnation pressure  Pt_secondary
//   6 t_e       secondary stagnation temperature Tt_secondary
//   7 p_out     outlet-node stagnation pressure Pt_outlet
//   8 p_py      the owned unknown, repurposed as the mixing-plane static P_py
struct EjectorElementResidualJacobian {
  std::array<double, 4> residuals;               // R0..R3
  std::array<std::array<double, 9>, 4> jacobian; // jacobian[row][seed]
};

EjectorElementResidualJacobian ejector_element_residuals_and_jacobian(
    double mp, double ms, double mdot_out, double p_g, double t_g, double p_e,
    double t_e, double p_out, double p_py, const EjectorGeometry& geom,
    double area_throat, double area_nozzle_exit, double area_secondary,
    double gamma, double r_gas, double eta_primary, double eta_secondary,
    double recovery_efficiency, double eps_frac, double s_choke_lo,
    double s_choke_hi, double s_sub_lo, double s_sub_hi);

} // namespace combaero::solver
