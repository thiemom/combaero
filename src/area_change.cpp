#include "area_change.h"

#include "math_constants.h"  // MSVC compatibility for M_PI

#include <algorithm>
#include <cmath>

namespace combaero {

namespace {

// Division-safe absolute value with smooth gradient through zero.
// m_abs = sqrt(m^2 + eps)  -- never exactly zero, always differentiable.
inline double smooth_abs(double m, double eps = 1e-12) {
    return std::sqrt(m * m + eps);
}

// d(smooth_abs)/dm = m / sqrt(m^2 + eps)
inline double smooth_abs_deriv(double m, double m_abs) {
    return m / m_abs;
}

// Logistic sigmoid:  sigma(x) = 1 / (1 + exp(-x))
// Returns (w, dw/dx) for x = m_dot / m_scale.
inline std::pair<double, double> sigmoid(double x) {
    // Clamp to avoid overflow in exp for large |x|.
    if (x > 30.0) return {1.0, 0.0};
    if (x < -30.0) return {0.0, 0.0};
    double e = std::exp(-x);
    double w = 1.0 / (1.0 + e);
    double dw = w * (1.0 - w);  // sigmoid derivative
    return {w, dw};
}

// Compressibility correction factor k1 = 1 + 0.35 * M_eff^2.
//
// M_eff is min(|Mach|, MACH_CLAMP).  The correction is calibrated for
// M < MACH_ACCURACY_LIMIT; above that it extrapolates but remains bounded,
// preventing unbounded k1 growth from destabilising the Newton solver when
// non-converged iterates produce unphysical Mach numbers.
//
// clamped is set true when the raw |Mach| exceeds MACH_CLAMP so the caller
// can propagate the diagnostic through AreaChangeResult::mach_clamped.
struct MachCorrection {
    double k1;
    bool   clamped;
};

inline MachCorrection mach_k1(double Mach) {
    double M_abs    = std::abs(Mach);
    bool   clamped  = M_abs > MACH_CLAMP;
    double M_eff    = clamped ? MACH_CLAMP : M_abs;
    return {1.0 + 0.35 * M_eff * M_eff, clamped};
}

} // anonymous namespace

AreaChangeResult sharp_area_change(double m_dot, double rho, double mu,
                                   double F0, double F1,
                                   double Mach, double m_scale, double D_h) {
    constexpr double eps = 1e-12;

    // --- Geometry (static, independent of flow direction) ---
    double F_small = std::min(F0, F1);
    double F_large = std::max(F0, F1);
    double ar = F_small / std::max(F_large, eps);  // area ratio in (0, 1]

    // Hydraulic diameter of the smaller section (reference for Re).
    // D_h > 0 overrides the circular assumption for non-circular ducts.
    double D_small = (D_h > 0.0) ? D_h : std::sqrt(4.0 * F_small / M_PI);

    // --- Reynolds number based on F_small (higher-velocity side) ---
    double m_abs = smooth_abs(m_dot);
    double v_small = m_abs / (rho * F_small + eps);
    double Re = rho * v_small * D_small / (mu + eps);

    // --- Loss coefficients (geometry-dependent, direction-independent) ---
    // Expansion (Borda-Carnot): zeta_exp = alpha - 2*ar + ar^2
    double alpha = 1.0 + 2.73 / std::pow(Re + 10.0, 0.4);
    double zeta_exp = alpha - 2.0 * ar + ar * ar;

    // Contraction (Weisbach sharp-edge): zeta_con = 0.5 * (1 - ar)
    double zeta_con = 0.5 * (1.0 - ar);

    // --- Assign forward/reverse zeta based on geometry ---
    // Positive m_dot means flow from F0 to F1.
    // If F0 < F1: forward flow sees expansion, reverse sees contraction.
    // If F0 > F1: forward flow sees contraction, reverse sees expansion.
    double zeta_fwd, zeta_rev;
    if (F0 <= F1) {
        // F0 -> F1 is an expansion
        zeta_fwd = zeta_exp;
        zeta_rev = zeta_con;
    } else {
        // F0 -> F1 is a contraction
        zeta_fwd = zeta_con;
        zeta_rev = zeta_exp;
    }

    // --- Sigmoid blending for smooth bidirectional transition ---
    auto [w, dw_dx] = sigmoid(m_dot / m_scale);
    double zeta_eff = w * zeta_fwd + (1.0 - w) * zeta_rev;

    // --- Compressibility correction ---
    auto [k1, mach_clamped] = mach_k1(Mach);

    // --- Pressure drop: dP = k1 * zeta_eff * m_dot * |m_dot| / (2 * rho * F_small^2) ---
    double K = k1 / (2.0 * rho * F_small * F_small + eps);
    double dP = K * zeta_eff * m_dot * m_abs;

    // --- Jacobians ---
    // dm_abs/dm = m / m_abs  (smooth_abs derivative)
    double dm_abs_dm = smooth_abs_deriv(m_dot, m_abs);

    // dw/dm = dw/dx * dx/dm = dw/dx / m_scale
    double dw_dm = dw_dx / m_scale;

    // d(zeta_exp)/dm via alpha(Re(m_dot)):
    //   alpha = 1 + 2.73 * (Re + 10)^(-0.4)
    //   d(alpha)/dRe = -2.73 * 0.4 * (Re + 10)^(-1.4)
    //   Re = m_abs * D_small / (mu * F_small)
    //   dRe/dm = D_small / (mu * F_small + eps) * dm_abs/dm
    double dalpha_dRe = -2.73 * 0.4 / std::pow(Re + 10.0, 1.4);
    double dRe_dm = D_small / (mu * F_small + eps) * dm_abs_dm;
    double dzeta_exp_dm = dalpha_dRe * dRe_dm;

    // d(zeta_eff)/dm = sigmoid blend + Re-dependent alpha
    //   = (zeta_fwd - zeta_rev) * dw/dm
    //     + w * dzeta_fwd/dm + (1-w) * dzeta_rev/dm
    // Only zeta_exp depends on Re; zeta_con does not.
    double dzeta_fwd_dm = (F0 <= F1) ? dzeta_exp_dm : 0.0;
    double dzeta_rev_dm = (F0 <= F1) ? 0.0 : dzeta_exp_dm;
    double dzeta_dm = (zeta_fwd - zeta_rev) * dw_dm
                      + w * dzeta_fwd_dm + (1.0 - w) * dzeta_rev_dm;

    // dP = K * zeta_eff * m_dot * m_abs
    // dS/dm = K * (zeta_eff * (m_abs + m_dot * dm_abs_dm) + dzeta_dm * m_dot * m_abs)
    double dS_dm = K * (zeta_eff * (m_abs + m_dot * dm_abs_dm)
                        + dzeta_dm * m_dot * m_abs);

    // dS/drho = -dP / rho  (since K ~ 1/rho, dP/drho = -dP/rho)
    double dS_drho = (std::abs(rho) > eps) ? -dP / rho : 0.0;

    // dS/dmu: mu affects dP through Re -> alpha -> zeta_exp
    //   dRe/dmu = -Re / mu
    //   dzeta_exp/dmu = dalpha_dRe * dRe/dmu
    double dRe_dmu = (std::abs(mu) > eps) ? -Re / mu : 0.0;
    double dzeta_exp_dmu = dalpha_dRe * dRe_dmu;
    double dzeta_fwd_dmu = (F0 <= F1) ? dzeta_exp_dmu : 0.0;
    double dzeta_rev_dmu = (F0 <= F1) ? 0.0 : dzeta_exp_dmu;
    double dzeta_eff_dmu = w * dzeta_fwd_dmu + (1.0 - w) * dzeta_rev_dmu;
    double dS_dmu = K * dzeta_eff_dmu * m_dot * m_abs;

    return {dP, dS_dm, dS_drho, dS_dmu, mach_clamped};
}

AreaChangeResult conical_area_change(double m_dot, double rho, double mu,
                                     double F0, double F1,
                                     double length,
                                     double Mach, double m_scale) {
    constexpr double eps    = 1e-12;
    constexpr double k_eta  = 3.2;   // angle-efficiency curve fit constants
    constexpr double n_eta  = 0.8;   // (Idelchik Fig. 5-2, circular diffusers)

    (void)mu;  // accepted for API uniformity; unused in conical loss model

    // --- Geometry (static, flow-independent) ---
    double F_small = std::min(F0, F1);
    double F_large = std::max(F0, F1);
    double ar = F_small / std::max(F_large, eps);  // area ratio in (0, 1]

    double r_small = std::sqrt(F_small / M_PI);
    double r_large = std::sqrt(F_large / M_PI);
    // Half-angle: atan2 keeps theta in [0, pi/2] for any positive length.
    double theta = std::atan2(r_large - r_small, std::max(length, eps));

    // --- Loss coefficients (purely geometric, no Re dependence) ---

    // Expansion (diffuser): eta rises from 0 at theta=0 to ~1 at theta=pi/2,
    // recovering Borda-Carnot in the sharp limit.
    double eta     = 1.0 - std::exp(-k_eta * std::pow(theta, n_eta));
    double zeta_exp = eta * (1.0 - ar) * (1.0 - ar);

    // Contraction (nozzle): C_con interpolates from 0.04 (gradual) to 0.50
    // (sharp, recovers Weisbach) with the half-angle.
    double sin_theta = std::sin(theta);
    double C_con     = 0.04 + 0.46 * std::pow(sin_theta, 1.5);
    double zeta_con  = C_con * (1.0 - ar);

    // --- Assign forward/reverse zeta based on geometry ---
    double zeta_fwd, zeta_rev;
    if (F0 <= F1) {
        // F0 -> F1 is an expansion (diffuser direction)
        zeta_fwd = zeta_exp;
        zeta_rev = zeta_con;
    } else {
        // F0 -> F1 is a contraction (nozzle direction)
        zeta_fwd = zeta_con;
        zeta_rev = zeta_exp;
    }

    // --- Sigmoid blending for smooth bidirectional transition ---
    auto [w, dw_dx] = sigmoid(m_dot / m_scale);
    double zeta_eff = w * zeta_fwd + (1.0 - w) * zeta_rev;

    // --- Compressibility correction ---
    auto [k1, mach_clamped] = mach_k1(Mach);

    // --- Smooth |m_dot| ---
    double m_abs = smooth_abs(m_dot);

    // --- Pressure drop: dP = k1 * zeta_eff * m_dot * |m_dot| / (2*rho*F_small^2) ---
    double K  = k1 / (2.0 * rho * F_small * F_small + eps);
    double dP = K * zeta_eff * m_dot * m_abs;

    // --- Jacobians ---
    // zeta_exp and zeta_con are purely geometric: dzeta_exp/dm = dzeta_con/dm = 0.
    // Therefore d(zeta_eff)/dm reduces to the sigmoid blend term only.
    double dw_dm      = dw_dx / m_scale;
    double dzeta_eff_dm = (zeta_fwd - zeta_rev) * dw_dm;  // Re chain absent

    // d(m_dot * m_abs)/dm = m_abs + m_dot * dm_abs/dm
    double dm_abs_dm  = smooth_abs_deriv(m_dot, m_abs);

    // dS/dm = K * (zeta_eff * (m_abs + m_dot*dm_abs_dm) + dzeta_eff_dm * m_dot*m_abs)
    double dS_dm = K * (zeta_eff * (m_abs + m_dot * dm_abs_dm)
                        + dzeta_eff_dm * m_dot * m_abs);

    // dS/drho = -dP / rho  (K ~ 1/rho)
    double dS_drho = (std::abs(rho) > eps) ? -dP / rho : 0.0;

    // dS/dmu = 0: mu does not appear in the conical loss model.
    double dS_dmu = 0.0;

    return {dP, dS_dm, dS_drho, dS_dmu, mach_clamped};

}

} // namespace combaero
