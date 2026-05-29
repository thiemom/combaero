#pragma once
#include "math_constants.h"
#include "correlation_status.h"
#include <cmath>
#include <complex>

// Tee junction loss coefficients (Bassett 2001, Proc. IMechE Part C 215(8):861-881).
//
// Notation:
//   q    = m_dot_branch / m_dot_com  (mass flow ratio, design range [0, 1])
//   psi  = F_C / F_B                 (area ratio: main duct / branch, > 0)
//   theta = branch angle [rad]       (validated range: (0, pi/2])
//   Ports: A = main inlet/outlet, B = lateral branch, C = common (total) port
//
// Accuracy limits (Bassett Section 4):
//   theta < pi/6:  joining-flow 1D/3D discrepancy can reach 100%+.
//   Ma < 0.2:      incompressible assumption required.
//   q in [0, 1]:   outside this range equations extrapolate beyond original derivation.
//
// Design for solver robustness:
//   - Raw K1-K12 functions: pure Bassett math, no clamping. Use for testing.
//   - Blend functions (merging_tee_K_*, branching_tee_K_*): apply soft bounds on
//     psi and the K11/K12 denominator to prevent NaN/Inf outside valid inputs.
//   - tee_check_inputs(): non-throwing validity check, returns CorrelationStatus.
//   - All functions are smooth (C-infinity) for all real inputs.

namespace combaero {

// -----------------------------------------------------------------------------
// Input validity
// -----------------------------------------------------------------------------

constexpr double TEE_Q_LO    = 0.0;         // design lower bound on flow ratio
constexpr double TEE_Q_HI    = 1.0;         // design upper bound on flow ratio
constexpr double TEE_PSI_MIN = 0.05;        // minimum area ratio (F_C/F_B)
constexpr double TEE_THETA_MAX = M_PI / 2.0; // max validated branch angle [rad]

struct TeeInputStatus {
    bool q_in_range;   // q in [TEE_Q_LO, TEE_Q_HI]
    bool psi_valid;    // psi >= TEE_PSI_MIN
    bool theta_valid;  // theta in (0, TEE_THETA_MAX]
    bool valid() const { return q_in_range && psi_valid && theta_valid; }
    CorrelationStatus status() const {
        return valid() ? CorrelationStatus::Valid : CorrelationStatus::Extrapolated;
    }
};

inline TeeInputStatus tee_check_inputs(double q, double psi, double theta) {
    return {
        q >= TEE_Q_LO && q <= TEE_Q_HI,
        psi >= TEE_PSI_MIN,
        theta > 0.0 && theta <= TEE_THETA_MAX
    };
}

// -----------------------------------------------------------------------------
// Smooth lower bound: soft_max(x, lo)
// Exact for x >> lo; smoothly approaches lo as x -> 0 or below.
// Derivative: 0.5 * (1 + (x-lo) / sqrt((x-lo)^2 + eps^2)) in [0, 1].
// eps controls transition width (1% of lo by default).
// -----------------------------------------------------------------------------
inline double soft_lower(double x, double lo) {
    // transition width 1% of lo avoids rounding issues for lo > 1e-4
    const double eps = lo * 0.01;
    const double dx = x - lo;
    return 0.5 * (x + lo + std::sqrt(dx * dx + eps * eps));
}

// -----------------------------------------------------------------------------
// Raw loss coefficients (Bassett 2001, Table 2) — pure math, no guards.
// These are provided for testing and diagnostics. Call with inputs validated
// by tee_check_inputs() for physically meaningful results.
// -----------------------------------------------------------------------------

// Separating flows (BranchingTee, flow types 1-3)

// K2/K5: straight-through loss, type 2/5. Independent of theta and psi (Eq. 15).
inline double K5(double q) { return q * q - 1.5 * q + 0.5; }
inline double dK5_dq(double q) { return 2.0 * q - 1.5; }

// K6: main-to-branch, type 3 (standard: branch at theta from main).
inline double K6(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * c;
}
inline double dK6_dq(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    return 2.0 * q * psi * psi - 2.0 * psi * c;
}

// K1: main-to-branch, type 1 (reversed orientation vs K6).
inline double K1(double q, double psi, double theta) {
    const double c = std::cos(0.75 * (M_PI - theta));
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * c;
}
inline double dK1_dq(double q, double psi, double theta) {
    const double c = std::cos(0.75 * (M_PI - theta));
    return 2.0 * q * psi * psi - 2.0 * psi * c;
}

// K4: branch-to-downstream, type 2. Has 1/psi: protect psi before calling.
inline double K4(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    return 1.0 + q * q / (psi * psi) - (2.0 * q / psi) * c;
}
inline double dK4_dq(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    return 2.0 * q / (psi * psi) - 2.0 * c / psi;
}

// K3: branch-to-downstream, type 2 reversed. Has 1/psi: protect psi before calling.
inline double K3(double q, double psi, double theta) {
    const double c = std::cos(0.75 * (M_PI - theta));
    return 1.0 + q * q / (psi * psi) - (2.0 * q / psi) * c;
}
inline double dK3_dq(double q, double psi, double theta) {
    const double c = std::cos(0.75 * (M_PI - theta));
    return 2.0 * q / (psi * psi) - 2.0 * c / psi;
}

// Joining flows (MergingTee, flow types 4-6)
// Angle corrections from Bassett Table 3 are embedded per function.

// K12: branch-to-outlet, type 6. D = psi + 0.5*cos(3/4*theta) in denominator.
// For theta <= pi/2 and psi > 0: D > 0. Protect if theta or psi is extrapolated.
inline double K12(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = psi + 0.5 * c;
    const double N = 1.0 - (1.0 - q) * (1.0 - q) - q * q * psi * c;
    return (2.0 * psi / D) * N + q * q * psi * psi - 1.0;
}
inline double dK12_dq(double q, double psi, double theta) {
    // dN/dq = 2*(1-q) - 2*q*psi*c  (chain rule through Bassett Eq. 12)
    const double c = std::cos(0.75 * theta);
    const double D = psi + 0.5 * c;
    const double dN = 2.0 * (1.0 - q) - 2.0 * q * psi * c;
    return (2.0 * psi / D) * dN + 2.0 * q * psi * psi;
}

// K11: straight-to-outlet, type 6. Same denominator D as K12.
inline double K11(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = psi + 0.5 * c;
    const double N = 1.0 - q * q - (1.0 - q) * (1.0 - q) * psi * c;
    return (2.0 * psi / D) * N + q * q - 1.0;
}
inline double dK11_dq(double q, double psi, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = psi + 0.5 * c;
    const double dN = -2.0 * q + 2.0 * (1.0 - q) * psi * c;
    return (2.0 * psi / D) * dN + 2.0 * q;
}

// K8: branch path, type 4. Has 1/psi: protect psi before calling.
// theta_corr = pi - (3/4)*(pi - theta).
inline double K8(double q, double psi, double theta) {
    const double tc = M_PI - 0.75 * (M_PI - theta);
    const double c = std::cos(tc);
    return 1.0 - q * q + 2.0 * (1.0 - q) * (1.0 - q) * c / psi;
}
inline double dK8_dq(double q, double psi, double theta) {
    const double tc = M_PI - 0.75 * (M_PI - theta);
    const double c = std::cos(tc);
    return -2.0 * q - 4.0 * (1.0 - q) * c / psi;
}

// K7: straight path, type 4. theta_corr = pi - (3/4)*(pi - theta).
inline double K7(double q, double psi, double theta) {
    const double tc = M_PI - 0.75 * (M_PI - theta);
    const double c = std::cos(tc);
    return 4.0 * q - 1.0 + q * q * (psi * psi - 2.0 + 2.0 * psi * c);
}
inline double dK7_dq(double q, double psi, double theta) {
    const double tc = M_PI - 0.75 * (M_PI - theta);
    const double c = std::cos(tc);
    return 4.0 + 2.0 * q * (psi * psi - 2.0 + 2.0 * psi * c);
}

// K9: type 5, first path. Has 1/psi: protect psi before calling. (Eq. 35)
// alpha = (pi-theta)/4, beta = theta/4.
inline double K9(double q, double psi, double theta) {
    const double alpha = (M_PI - theta) / 4.0;
    const double beta = theta / 4.0;
    const double ca = std::cos(theta + alpha);
    const double cb = std::cos(theta - beta);
    return 2.0 * q * q * ca / psi - 2.0 * (1.0 - q) * (1.0 - q) * cb / psi
           + q * q / (psi * psi) + 1.0;
}
inline double dK9_dq(double q, double psi, double theta) {
    const double alpha = (M_PI - theta) / 4.0;
    const double beta = theta / 4.0;
    const double ca = std::cos(theta + alpha);
    const double cb = std::cos(theta - beta);
    return 4.0 * q * ca / psi + 4.0 * (1.0 - q) * cb / psi + 2.0 * q / (psi * psi);
}

// K10: type 5, second path. Has 1/psi: protect psi before calling. (Eq. 37)
inline double K10(double q, double psi, double theta) {
    const double alpha = (M_PI - theta) / 4.0;
    const double beta = theta / 4.0;
    const double ca = std::cos(theta + alpha);
    const double cb = std::cos(theta - beta);
    return 2.0 * (1.0 - q) * (1.0 - q) * ca / psi - 2.0 * q * q * cb / psi
           + q * q / (psi * psi) + 1.0;
}
inline double dK10_dq(double q, double psi, double theta) {
    const double alpha = (M_PI - theta) / 4.0;
    const double beta = theta / 4.0;
    const double ca = std::cos(theta + alpha);
    const double cb = std::cos(theta - beta);
    return -4.0 * (1.0 - q) * ca / psi - 4.0 * q * cb / psi + 2.0 * q / (psi * psi);
}

// -----------------------------------------------------------------------------
// Smooth topology blend
// w(r) = 0.5*(1+tanh(k*r)): 1 at r>>0 (design direction), 0 at r<<0 (reversed).
// k=30 gives a transition width of ~0.1 in q around r=0.
// -----------------------------------------------------------------------------
inline double blend_weight(double r, double k = 30.0) {
    return 0.5 * (1.0 + std::tanh(k * r));
}
inline double blend_weight_deriv(double r, double k = 30.0) {
    // d/dr[0.5*(1+tanh(k*r))] = 0.5*k*sech^2(k*r) = 2*k*w*(1-w)
    const double w = blend_weight(r, k);
    return 2.0 * k * w * (1.0 - w);
}

// -----------------------------------------------------------------------------
// Internal safe helpers for K11/K12 with soft-bounded psi and denominator.
// Used by the blend functions. Not intended for direct user calls.
// -----------------------------------------------------------------------------
namespace detail {

// TEE_PSI_MIN * 0.1 as denominator floor: handles extreme theta extrapolation.
constexpr double D_FLOOR = TEE_PSI_MIN * 0.1;

inline double K11_s(double q, double ps, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = soft_lower(ps + 0.5 * c, D_FLOOR);
    const double N = 1.0 - q * q - (1.0 - q) * (1.0 - q) * ps * c;
    return (2.0 * ps / D) * N + q * q - 1.0;
}
inline double dK11_s_dq(double q, double ps, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = soft_lower(ps + 0.5 * c, D_FLOOR);
    const double dN = -2.0 * q + 2.0 * (1.0 - q) * ps * c;
    return (2.0 * ps / D) * dN + 2.0 * q;
}
inline double K12_s(double q, double ps, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = soft_lower(ps + 0.5 * c, D_FLOOR);
    const double N = 1.0 - (1.0 - q) * (1.0 - q) - q * q * ps * c;
    return (2.0 * ps / D) * N + q * q * ps * ps - 1.0;
}
inline double dK12_s_dq(double q, double ps, double theta) {
    const double c = std::cos(0.75 * theta);
    const double D = soft_lower(ps + 0.5 * c, D_FLOOR);
    const double dN = 2.0 * (1.0 - q) - 2.0 * q * ps * c;
    return (2.0 * ps / D) * dN + 2.0 * q * ps * ps;
}

} // namespace detail

// -----------------------------------------------------------------------------
// MergingTee effective K (type 6 primary; K5/K6 for reversed topology).
// Applies soft lower bound on psi and K11/K12 denominator for robustness.
// -----------------------------------------------------------------------------

inline double merging_tee_K_straight(double q, double psi, double theta,
                                      double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    return w * detail::K11_s(q, ps, theta) + (1.0 - w) * K5(q);
}
inline double merging_tee_dK_straight_dq(double q, double psi, double theta,
                                          double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    const double dw = blend_weight_deriv(q, blend_k);
    const double K11v = detail::K11_s(q, ps, theta);
    const double K5v = K5(q);
    return dw * (K11v - K5v)
           + w * detail::dK11_s_dq(q, ps, theta) + (1.0 - w) * dK5_dq(q);
}

inline double merging_tee_K_branch(double q, double psi, double theta,
                                    double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    return w * detail::K12_s(q, ps, theta) + (1.0 - w) * K6(q, ps, theta);
}
inline double merging_tee_dK_branch_dq(double q, double psi, double theta,
                                        double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    const double dw = blend_weight_deriv(q, blend_k);
    const double K12v = detail::K12_s(q, ps, theta);
    const double K6v = K6(q, ps, theta);
    return dw * (K12v - K6v)
           + w * detail::dK12_s_dq(q, ps, theta) + (1.0 - w) * dK6_dq(q, ps, theta);
}

// -----------------------------------------------------------------------------
// BranchingTee effective K (type 3 primary; K11/K12 for reversed topology).
// Applies same soft bounds as MergingTee for consistent behavior.
// -----------------------------------------------------------------------------

inline double branching_tee_K_straight(double q, double psi, double theta,
                                        double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    return w * K5(q) + (1.0 - w) * detail::K11_s(q, ps, theta);
}
inline double branching_tee_dK_straight_dq(double q, double psi, double theta,
                                            double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    const double dw = blend_weight_deriv(q, blend_k);
    const double K5v = K5(q);
    const double K11v = detail::K11_s(q, ps, theta);
    return dw * (K5v - K11v)
           + w * dK5_dq(q) + (1.0 - w) * detail::dK11_s_dq(q, ps, theta);
}

inline double branching_tee_K_branch(double q, double psi, double theta,
                                      double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    return w * K6(q, ps, theta) + (1.0 - w) * detail::K12_s(q, ps, theta);
}
inline double branching_tee_dK_branch_dq(double q, double psi, double theta,
                                          double blend_k = 30.0) {
    const double ps = soft_lower(psi, TEE_PSI_MIN);
    const double w = blend_weight(q, blend_k);
    const double dw = blend_weight_deriv(q, blend_k);
    const double K6v = K6(q, ps, theta);
    const double K12v = detail::K12_s(q, ps, theta);
    return dw * (K6v - K12v)
           + w * dK6_dq(q, ps, theta) + (1.0 - w) * detail::dK12_s_dq(q, ps, theta);
}

// -----------------------------------------------------------------------------
// Diagnostic utilities
// -----------------------------------------------------------------------------

// Returns true if q is in the design topology window [-epsilon, 1+epsilon].
inline bool tee_topology_valid(double q, double epsilon = 0.05) {
    return q >= -epsilon && q <= 1.0 + epsilon;
}

// Signed flow ratio q = m_dot_branch / m_dot_com, regularised identically to
// the solver: q -> 0 as m_dot_com -> 0 (stalled flow has no meaningful ratio).
inline double tee_flow_ratio(double m_dot_branch, double m_dot_com) {
    constexpr double eps_m = 1e-10; // kg/s, matches solver regularisation
    return m_dot_branch * m_dot_com / (m_dot_com * m_dot_com + eps_m * eps_m);
}

// =============================================================================
// Unified0D compressible junction model
// Mynard & Valen-Sendstad 2015, extended to O(M^2) (see docs/junction/).
// =============================================================================

// Per-branch face state, templated on scalar T (double or complex<double>).
// mdot sign: >0 = supplier (inflow to junction), <0 = collector (outflow).
// gamma_eff and R_gas_eff are effective perfect-gas constants for this branch.
template <typename T>
struct BranchState {
    T p;              // static pressure [Pa]
    T Tk;             // static temperature [K]
    T mdot;           // signed mass flow [kg/s]
    double A;         // cross-section area [m^2]
    double theta;     // centreline angle [rad] (measured from reference direction)
    double gamma_eff; // ratio of specific heats [-]
    double R_gas_eff; // specific gas constant [J/(kg*K)]

    T rho()  const { return p / (R_gas_eff * Tk); }
    T u()    const { return mdot * R_gas_eff * Tk / (p * A); }
    T M2()   const { auto v = u(); return v * v / (gamma_eff * R_gas_eff * Tk); }
    T Phi()  const { return T{1.0} + T{0.5 * (gamma_eff - 1.0)} * M2(); }
    T h0()   const {
        const double cp_eff = gamma_eff * R_gas_eff / (gamma_eff - 1.0);
        auto v = u();
        return T{cp_eff} * Tk + T{0.5} * v * v;
    }
    T p0() const {
        return p * std::pow(Phi(), T{gamma_eff / (gamma_eff - 1.0)});
    }
};

// Closed-form K_dat,j to O(M^2), Variant-2 (momentum-consistent) expansion;
// see Section 5 (eq. Kclosed) of docs/junction/junction_model_v2.tex.
//   x_j   = lambda_j * psi_j  (flow-area ratio product, >= 0)
//   phi_j  = effective inflow angle [rad]  (Hager 3/4 factor already applied)
//   M_dat2 = M_dat^2 of the pseudodatum
//
// kappa = (s/2) * [-1 - 2*x/sqrt(K_inc) + s/(2*K_inc)] with s = mu0^2 - 1 and
// mu0 = x + sqrt(K_inc). The bracket is <= 0 on the entire physical domain, so
// kappa <= 0: compressibility lowers K at fixed (lambda, psi, phi). The earlier
// frozen-xi form kappa = s/2 dropped the bracket and had the wrong sign, which
// drove the Newton residual against its Jacobian (v1 convergence failures).
template <typename T>
T K_dat_j_closed(T x_j, T phi_j, T M_dat2) {
    const T cosphi  = std::cos(phi_j);
    const T K_inc   = T{1.0} + x_j * x_j - T{2.0} * x_j * cosphi;
    const T sqrtK   = std::sqrt(K_inc);
    const T mu0     = x_j + sqrtK;
    const T s       = mu0 * mu0 - T{1.0};
    const T bracket = T{-1.0} - T{2.0} * x_j / sqrtK + s / (T{2.0} * K_inc);
    const T kappa   = T{0.5} * s * bracket;
    return K_inc * (T{1.0} + kappa * M_dat2);
}

// Machine-precision derivative of K_dat_j_closed (or any compatible kernel)
// w.r.t. one of its three arguments, evaluated at a real point.
//   which: 0 = d/dx, 1 = d/dphi, 2 = d/dM2
//   h:     complex-step size (default 1e-30 gives machine-precision results)
template <typename Fn>
double cstep_deriv(Fn f, int which, double x, double phi, double M2,
                   double h = 1e-30)
{
    using C = std::complex<double>;
    C xC{x}, pC{phi}, M2C{M2};
    switch (which) {
        case 0: xC  = C(x,   h); break;
        case 1: pC  = C(phi, h); break;
        case 2: M2C = C(M2,  h); break;
        default: break;
    }
    return std::imag(f(xC, pC, M2C)) / h;
}

} // namespace combaero
