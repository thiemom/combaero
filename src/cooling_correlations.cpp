#include "cooling_correlations.h"
#include "dual_number.h"
#include "heat_transfer.h"
#include "math_constants.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace combaero::cooling {

double thermal_performance_factor(double Nu_ratio, double f_ratio) {
    // Webb & Eckert (1972) thermal-hydraulic performance factor.
    // eta = (Nu/Nu0) / (f/f0)^(1/3)
    // Compares heat transfer gain against pumping power penalty at equal velocity.
    // eta > 1: net improvement over smooth; eta < 1: pressure-drop generator.
    //
    // Cross-check against Singh & Ekkad tabulated data (e/D=0.0625, P/e=10, 90 deg):
    //   Re=30k:  2.45 / 6.50^(1/3) = 2.45 / 1.866 = 1.31  (table: 1.31)
    //   Re=100k: 1.95 / 6.30^(1/3) = 1.95 / 1.849 = 1.05  (table: 1.05)
    //   Re=400k: 1.55 / 6.20^(1/3) = 1.55 / 1.838 = 0.84  (table: 0.84)
    //
    // Reference: Webb, R.L. & Eckert, E.R.G. (1972)
    //   Int. J. Heat Mass Transfer 15(8), 1647-1658
    return Nu_ratio / std::pow(f_ratio, 1.0 / 3.0);
}

double adiabatic_wall_temperature(double T_hot, double T_coolant, double eta) {
    if (eta < 0.0 || eta > 1.0) {
        throw std::runtime_error(
            "Cooling effectiveness eta must be in range [0, 1], got " + std::to_string(eta)
        );
    }
    return T_hot - eta * (T_hot - T_coolant);
}

double cooled_wall_heat_flux(double T_hot, double T_coolant,
                             double h_hot, double h_coolant,
                             double eta,
                             double t_wall, double k_wall) {
    if (h_hot <= 0.0 || h_coolant <= 0.0) {
        throw std::runtime_error("Heat transfer coefficients must be positive");
    }
    if (t_wall < 0.0 || k_wall <= 0.0) {
        throw std::runtime_error("Wall thickness must be non-negative and conductivity positive");
    }

    const double T_aw = adiabatic_wall_temperature(T_hot, T_coolant, eta);
    const double U = overall_htc_wall(h_hot, h_coolant, t_wall, k_wall);
    return U * (T_aw - T_coolant);
}

// -------------------------------------------------------------
// Baldauf et al. (2002) film-cooling effectiveness
// -------------------------------------------------------------

namespace {

using combaero::solver::DualN;

// A power whose EXPONENT may itself carry derivatives. dpow() in
// dual_number.h takes a constant exponent, but Eqs. (37) and (40) raise to
// 1/eta_s, b_1 c_1 and xi_s, all of which depend on M and P. a^b =
// exp(b ln a) keeps that exact.
inline double pow_tt(double a, double b) { return std::pow(a, b); }

template <int N>
inline DualN<N> pow_tt(const DualN<N>& a, const DualN<N>& b) {
    return combaero::solver::dexp(b * combaero::solver::dlog(a));
}

inline double dpow_t(double a, double e) { return std::pow(a, e); }
template <int N>
inline DualN<N> dpow_t(const DualN<N>& a, double e) {
    return combaero::solver::dpow(a, e);
}

inline double as_const(double v, double) { return v; }
template <int N>
inline DualN<N> as_const(double v, const DualN<N>&) { return DualN<N>::constant(v); }

// The whole correlation, written once. Only M and P are templated: geometry
// and turbulence are fixed per element, so they stay double and contribute
// no partials. One chain of equations means the derivative cannot drift out
// of sync with the value.
template <typename T>
T baldauf_eta_impl(double x_over_D, const T& M, const T& P, double alpha_deg,
                   double sD, double Tu) {
    namespace B = combaero::cooling::baldauf2002;

    // Degrees in, radians inside every trig function -- established from
    // seven of the paper's own Table 4 coefficients, not assumed.
    const double a = alpha_deg * M_PI / 180.0;

    const T U = M / P;                                                   // Eq. (16)
    const T one = as_const(1.0, M);

    // --- geometry-only coefficients (all plain doubles) -------------------
    const double xi_c = 0.6 + 0.4 * (2.0 - std::cos(a))
                              / (1.0 + std::pow((sD - 1.0) / 3.3, 6.0));  // Eq. (18)
    const double a_pk = 0.2;                                              // Eq. (11)
    const double b_pk = std::exp(1.92 - 7.5 * std::pow(sD, -1.5));        // Eq. (12)
    const double c_pk = 0.7 + 336.0 * std::exp(-1.85 * sD);               // Eq. (13)
    const double mu_0 = 0.125 + 0.063 * std::pow(sD, 1.8);                // Eq. (14)
    const double eta_c0 = 0.465 / (1.0 + 0.048 * sD * sD);                // Eq. (15)
    const double g = 0.75 * (1.0 - std::exp(-0.8 * (sD - 1.0)));          // Eq. (24)
    const double k = 2.0 * (1.0 - std::exp(0.57 * (1.0 - sD)))
                   + 0.91 * std::pow(std::cos(a), 0.65);                  // Eq. (25)
    const double xi_hat =
        1.17 * (1.0 - (sD - 1.0) / (1.0 + 0.2 * (sD - 1.0) * (sD - 1.0)))
             * (std::cos(2.3 * a) + 2.45);                                // Eq. (26)
    const double eta_hat = 0.022 * (sD + 1.0) * (0.9 - std::sin(2.0 * a))
                         - (0.08 + 0.46 / (1.0 + (sD - 3.2) * (sD - 3.2)));  // Eq. (27)

    // Eq. (31), IMPLEMENTED AS PRINTED. It disagrees with the paper's own
    // Table 4 by 36% (0.83612467 against 0.61626073) while all 17 other
    // checkable coefficients reproduce to 2e-6. See the header note.
    const double b_0 = 0.8 - 0.014 * sD * sD
                     + (1.5 - 2.0 / std::sqrt(sD))
                       * std::sin(0.86 * a * (1.0 + 0.754 / (1.0 + 0.87 * sD * sD)));
    const double c_1 = 7.5 + sD;                                          // Eq. (33)

    const double b_T = 0.7 * (1.0 + (1.22 / (1.0 + 7.0 * std::pow(sD - 1.0, -7.0))
                                     + 0.87 + std::cos(2.5 * a))
                              * std::exp(2.6 * Tu - 0.0012 / (Tu * Tu) - 1.76));  // (34)
    const double eta_0T = B::eta_0T_c0
        * std::pow(B::eta_0T_c1 / B::eta_0T_c0, b_T / 0.7);               // Eq. (35)

    // --- flow-dependent coefficients --------------------------------------
    const T mu = U * dpow_t(P, 0.8)
               * (one - as_const((0.03 + 0.11 * (5.0 - sD)) * std::cos(a), M));  // (9)

    // Shared transition weight of Eqs. (22) and (23).
    const T tr = one / (one + dpow_t(U * dpow_t(P, g) / k, -5.0));
    const T xi_s = one + xi_hat * tr;                                     // Eq. (22)
    const T eta_s = one + eta_hat * tr;                                   // Eq. (23)
    const T a_1 = as_const(0.04 + 0.23 * sD, M)
                + as_const((0.95 - 0.19 * sD) * std::cos(1.5 * a), M);    // Eq. (30)
    const T xi_1 = 65.0 / pow_tt(M / 2.5, a_1);                           // Eq. (29)
    const T b_1 = b_0 / (one + dpow_t(M, -3.0));                          // Eq. (32)

    // Eq. (40) inverted in closed form -- no iteration:
    //   x/D = xi'^(1/xi_s) (pi/4) U^E / ((s/D) xi_c)
    const double E = std::pow(sD / 3.0, -0.75);
    const T base = (x_over_D * sD * xi_c) / ((M_PI / 4.0) * dpow_t(U, E));
    const T xi_p = pow_tt(base, xi_s);

    // Eq. (36): the turbulence-dependent base curve.
    const T r0 = xi_p / B::xi_0;
    const T eta_sp = eta_0T * dpow_t(r0, B::a_star)
        / dpow_t(one + dpow_t(r0, (B::a_star + b_T) * B::c_star), 1.0 / B::c_star);

    // Eq. (37): undo the adjacent-jet downstream branch.
    const T r1 = xi_p / xi_1;
    const T eta_st = 0.1 * pow_tt(eta_sp / 0.1, one / eta_s)
                   * pow_tt(one + pow_tt(r1, b_1 * c_1), one / (one * c_1));

    // Eq. (38): blow up to this case's peak.
    const T rm = mu / mu_0;
    const T eta_c = eta_c0 * eta_st * dpow_t(rm, a_pk)
        / dpow_t(one + dpow_t(rm, (a_pk + b_pk) * c_pk), 1.0 / c_pk);

    // Eq. (39): backscale out the geometry normalisation.
    return eta_c * dpow_t(P, 0.9 / sD)
         / as_const(std::pow(std::sin(a), 0.06 * sD), M);
}

void check_inputs(double x_over_D, double M, double P, double alpha_deg,
                  double sD, double Tu) {
    if (x_over_D <= 0.0) {
        throw std::invalid_argument(
            "film_effectiveness_baldauf_2002: x/D must be positive; the "
            "correlation is defined downstream of the ejection point");
    }
    if (M <= 0.0 || P <= 0.0 || sD <= 0.0) {
        throw std::invalid_argument(
            "film_effectiveness_baldauf_2002: M, P and s/D must be positive");
    }
    if (Tu <= 0.0) {
        // Eq. (34) carries exp(-0.0012/Tu^2), which underflows to a
        // different branch rather than erroring, so catch it here.
        throw std::invalid_argument(
            "film_effectiveness_baldauf_2002: Tu must be positive; Eq. (34) "
            "contains -0.0012/Tu^2 and is singular at zero turbulence");
    }
    if (alpha_deg <= 0.0 || alpha_deg > 90.0) {
        throw std::invalid_argument(
            "film_effectiveness_baldauf_2002: alpha_deg is the ejection angle "
            "to the SURFACE and must lie in (0, 90]");
    }
}

}  // anonymous namespace

double film_effectiveness_baldauf_2002(double x_over_D, double M, double P,
                                       double alpha_deg, double s_over_D,
                                       double Tu) {
    check_inputs(x_over_D, M, P, alpha_deg, s_over_D, Tu);
    return baldauf_eta_impl<double>(x_over_D, M, P, alpha_deg, s_over_D, Tu);
}

std::tuple<double, double, double> film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu) {
    check_inputs(x_over_D, M, P, alpha_deg, s_over_D, Tu);
    using D = DualN<2>;
    const D dM = D::seed(M, 0);
    const D dP = D::seed(P, 1);
    const D r = baldauf_eta_impl<D>(x_over_D, dM, dP, alpha_deg, s_over_D, Tu);
    return {r.v, r.d[0], r.d[1]};
}

}  // namespace combaero::cooling
