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
// Scope: critical (double-choked) operation only. The subcritical,
// unchoked-primary, and back-flow regimes are not modeled here; their design
// (choked/unchoked primary R0, subcritical omega droop to the dead-head
// pressure P_b0, and a subsonic jet-pump mode) lives in
// validation/ejector/OPERATING_REGIMES_DESIGN.md.
//
// Jacobian scope (this header): analytic derivatives w.r.t. the 4
// thermodynamic inputs (p_g, t_g, p_e, t_e) only -- these are the only
// inputs that are ever Newton-solved unknowns for the network element
// consuming this module (EjectorElement). gamma, r_gas, geometry, and
// recovery_efficiency are element parameters, not solved unknowns, so
// their derivatives are intentionally not exposed here (deliberate scope
// boundary, not an oversight -- validation/ejector/data/
// huang1999_reference_data.h separately carries central-difference targets
// for those too, for a possible future design-sensitivity extension).
//
// Key simplification that makes every derivative below exact analytic
// chain rule (no implicit/root-finding differentiation needed anywhere):
// the one internal root-find (ejector_mach_from_area_ratio_supersonic, a
// bisection for the primary nozzle-exit Mach number M_p1) depends ONLY on
// area_ratio_nozzle and gamma -- neither a Newton-solved unknown here -- so
// M_p1 has a PROVABLY exact zero derivative w.r.t. p_g, t_g, p_e, t_e.  It
// is computed once (bisection, same as Python) and then treated as a
// precomputed constant throughout the rest of the (otherwise fully
// closed-form) chain.
// -----------------------------------------------------------------------------

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

} // namespace combaero::solver
