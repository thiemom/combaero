#include "cooling_correlations.h"
#include "dual_number.h"
#include "heat_transfer.h"
#include "correlation_status.h"
#include "math_constants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

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
                   double sD, double Tu, double b_0_override){
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
    //
    // A caller may substitute its own b_0 to test the Table 4 reading. The
    // library offers no second FORMULA because the paper contains none:
    // Table 4 states one number at the Table 3 conditions, and reaching it
    // would need sin(...) = -0.167 where the printed equation gives +0.470
    // -- no sign, unit or angle convention makes that argument negative.
    // Scaling Eq. (31) by the 0.737 ratio would be an invented correction,
    // so the choice is left with whoever has evidence for it.
    const double b_0_printed = 0.8 - 0.014 * sD * sD
                     + (1.5 - 2.0 / std::sqrt(sD))
                       * std::sin(0.86 * a * (1.0 + 0.754 / (1.0 + 0.87 * sD * sD)));
    const double b_0 = std::isnan(b_0_override) ? b_0_printed : b_0_override;
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

// Warn for each input outside Baldauf's stated envelope, and report whether
// any was.
//
// The bounds have been `constexpr` in the header since the correlation
// landed, under a comment saying that outside them it is extrapolation --
// but nothing read them, so the correlation answered a 7.4 hole spacing as
// confidently as a 3.0 one. That matters beyond a user's own judgement: the
// validation harness derives its `extrapolated` column from this signal, so
// a silent correlation makes extrapolation indistinguishable from model
// error in the scorecard.
//
// `status` follows nusselt_dittus_boelter's contract -- when a caller passes
// one, it takes the flag and no warning is emitted; a null caller gets the
// warnings.
bool check_range(double M, double P, double alpha_deg, double sD, double Tu,
                 CorrelationStatus *status) {
    bool extrapolated = false;

    const auto flag = [&](bool outside, const char *name, double value,
                          double lo, double hi) {
        if (!outside) {
            return;
        }
        extrapolated = true;
        if (!status) {
            warn("film_effectiveness_baldauf_2002: " + std::string(name) +
                 " = " + std::to_string(value) + " is outside validated range ["
                 + std::to_string(lo) + ", " + std::to_string(hi) +
                 "]. Extrapolating; check results.");
        }
    };

    flag(M < baldauf2002::M_min || M > baldauf2002::M_max, "M", M,
         baldauf2002::M_min, baldauf2002::M_max);
    flag(P < baldauf2002::P_min || P > baldauf2002::P_max, "P", P,
         baldauf2002::P_min, baldauf2002::P_max);
    flag(alpha_deg < baldauf2002::alpha_deg_min
             || alpha_deg > baldauf2002::alpha_deg_max,
         "alpha_deg", alpha_deg, baldauf2002::alpha_deg_min,
         baldauf2002::alpha_deg_max);
    flag(sD < baldauf2002::s_over_D_min || sD > baldauf2002::s_over_D_max,
         "s_over_D", sD, baldauf2002::s_over_D_min,
         baldauf2002::s_over_D_max);
    flag(Tu < baldauf2002::Tu_min || Tu > baldauf2002::Tu_max, "Tu", Tu,
         baldauf2002::Tu_min, baldauf2002::Tu_max);

    if (status) {
        *status = extrapolated ? CorrelationStatus::Extrapolated
                               : CorrelationStatus::Valid;
    }
    return extrapolated;
}

}  // anonymous namespace

double film_effectiveness_baldauf_2002(double x_over_D, double M, double P,
                                       double alpha_deg, double s_over_D,
                                       double Tu,
                                       CorrelationStatus *status,
                                       double b_0_override) {
    check_inputs(x_over_D, M, P, alpha_deg, s_over_D, Tu);
    check_range(M, P, alpha_deg, s_over_D, Tu, status);
    return baldauf_eta_impl<double>(x_over_D, M, P, alpha_deg, s_over_D, Tu,
                                    b_0_override);
}

std::tuple<double, double, double> film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status, double b_0_override) {
    check_inputs(x_over_D, M, P, alpha_deg, s_over_D, Tu);
    check_range(M, P, alpha_deg, s_over_D, Tu, status);
    using D = DualN<2>;
    const D dM = D::seed(M, 0);
    const D dP = D::seed(P, 1);
    const D r = baldauf_eta_impl<D>(x_over_D, dM, dP, alpha_deg, s_over_D, Tu,
                                    b_0_override);
    return {r.v, r.d[0], r.d[1]};
}

// -------------------------------------------------------------
// Multi-row film superposition
// -------------------------------------------------------------

namespace {

void check_rows(const std::vector<double>& eta_rows) {
    if (eta_rows.empty()) {
        throw std::invalid_argument(
            "film superposition: at least one row is required");
    }
    for (double e : eta_rows) {
        if (!(e >= 0.0 && e <= 1.0)) {
            throw std::invalid_argument(
                "film superposition: each row effectiveness must lie in "
                "[0, 1]; eta is a normalised temperature difference");
        }
    }
}

void check_alphas(const std::vector<double>& eta_rows,
                  const std::vector<double>& alphas) {
    if (alphas.size() + 1 != eta_rows.size()) {
        throw std::invalid_argument(
            "film superposition: alpha_between_rows must have exactly one "
            "fewer entry than eta_rows -- alpha_j is the correction applied "
            "between row j and row j+1");
    }
    for (double a : alphas) {
        if (!(a >= 0.0 && a <= 1.0)) {
            throw std::invalid_argument(
                "film superposition: each alpha must lie in [0, 1]. It is "
                "the fraction of the film's temperature deficit surviving "
                "mixing to the next row; 1 is uncorrected Sellers, and a "
                "value above 1 would create coolant.");
        }
    }
}

}  // anonymous namespace

EffusionPanelFilm effusion_panel_film_effectiveness(int n_rows,
                                                    double pitch_x_over_D,
                                                    double s_over_D, double M,
                                                    double density_ratio,
                                                    double alpha_deg, double Tu) {
  if (n_rows < 1 || !(pitch_x_over_D > 0.0)) {
    throw std::invalid_argument(
        "effusion_panel_film_effectiveness: n_rows >= 1 and pitch_x_over_D > 0");
  }
  EffusionPanelFilm out;
  double total = 0.0;
  std::vector<double> rows;
  rows.reserve(static_cast<std::size_t>(n_rows));
  for (int n = 1; n <= n_rows; ++n) {
    rows.clear();
    for (int j = 1; j <= n; ++j) {
      CorrelationStatus status = CorrelationStatus::Valid;
      rows.push_back(film_effectiveness_baldauf_2002(j * pitch_x_over_D, M, density_ratio,
                                                     alpha_deg, s_over_D, Tu, &status));
      if (status != CorrelationStatus::Valid) {
        out.extrapolated = true;
      }
    }
    total += film_superposition_sellers(rows);
  }
  out.eta = total / n_rows;
  return out;
}

double film_superposition_sellers(const std::vector<double>& eta_rows) {
    check_rows(eta_rows);
    double remaining = 1.0;
    for (double e : eta_rows) {
        remaining *= (1.0 - e);
    }
    return 1.0 - remaining;
}

std::pair<double, std::vector<double>> film_superposition_sellers_and_gradient(
    const std::vector<double>& eta_rows) {
    check_rows(eta_rows);
    const std::size_t n = eta_rows.size();

    // Prefix and suffix products of (1 - eta), so the "product excluding i"
    // comes out without dividing -- a fully effective row would otherwise
    // divide by zero.
    std::vector<double> pre(n + 1, 1.0), suf(n + 1, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
        pre[i + 1] = pre[i] * (1.0 - eta_rows[i]);
    }
    for (std::size_t i = n; i-- > 0;) {
        suf[i] = suf[i + 1] * (1.0 - eta_rows[i]);
    }

    std::vector<double> grad(n);
    for (std::size_t i = 0; i < n; ++i) {
        grad[i] = pre[i] * suf[i + 1];
    }
    return {1.0 - pre[n], grad};
}

double film_superposition_corrected(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows) {
    check_rows(eta_rows);
    check_alphas(eta_rows, alpha_between_rows);
    const std::size_t n = eta_rows.size();

    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        // prod_{j=i..n-2} alpha_j  (zero-based: alpha_j sits between rows
        // j and j+1, so the film from row i passes alphas i .. n-2)
        double a = 1.0;
        for (std::size_t j = i; j + 1 < n; ++j) {
            a *= alpha_between_rows.at(j);
        }
        // prod_{k=i+1..n-1} (1 - eta_k)
        double b = 1.0;
        for (std::size_t k = i + 1; k < n; ++k) {
            b *= (1.0 - eta_rows[k]);
        }
        total += eta_rows[i] * a * b;
    }
    return total;
}

std::pair<double, std::vector<double>> film_superposition_corrected_and_gradient(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows) {
    check_rows(eta_rows);
    check_alphas(eta_rows, alpha_between_rows);
    const std::size_t n = eta_rows.size();

    std::vector<double> A(n, 1.0);   // prod of alphas from row i onward
    for (std::size_t i = 0; i < n; ++i) {
        double a = 1.0;
        for (std::size_t j = i; j + 1 < n; ++j) {
            a *= alpha_between_rows.at(j);
        }
        A[i] = a;
    }
    std::vector<double> suf(n + 1, 1.0);   // prod_{k>=i} (1 - eta_k)
    for (std::size_t i = n; i-- > 0;) {
        suf[i] = suf[i + 1] * (1.0 - eta_rows[i]);
    }

    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        total += eta_rows[i] * A[i] * suf[i + 1];
    }

    // d eta / d eta_m = A_m suf_{m+1}
    //                 - sum_{i<m} eta_i A_i prod_{k=i+1..n-1, k!=m} (1-eta_k)
    // The inner product is formed by an explicit loop rather than dividing
    // suf by (1 - eta_m), which would blow up for a fully effective row.
    std::vector<double> grad(n, 0.0);
    for (std::size_t m = 0; m < n; ++m) {
        double g = A[m] * suf[m + 1];
        for (std::size_t i = 0; i < m; ++i) {
            double prod = 1.0;
            for (std::size_t k = i + 1; k < n; ++k) {
                if (k != m) {
                    prod *= (1.0 - eta_rows[k]);
                }
            }
            g -= eta_rows[i] * A[i] * prod;
        }
        grad[m] = g;
    }
    return {total, grad};
}

double mainstream_temperature_correction(double mass_flow_ratio, double a,
                                         double b) {
    if (mass_flow_ratio < 0.0) {
        throw std::invalid_argument(
            "mainstream_temperature_correction: mass_flow_ratio must be "
            "non-negative");
    }
    const double ar = a * mass_flow_ratio;
    if (ar + 1.0 <= 0.0) {
        throw std::invalid_argument(
            "mainstream_temperature_correction: a * mass_flow_ratio + 1 must "
            "be positive");
    }
    return ar / (ar + 1.0) + b;
}

double equivalent_slot_width(double hole_area, double pitch) {
    if (hole_area <= 0.0 || pitch <= 0.0) {
        throw std::invalid_argument(
            "equivalent_slot_width: hole_area and pitch must be positive");
    }
    return hole_area / pitch;
}

double equivalent_blowing_ratio(double M_baseline, double area_baseline,
                                double area_equivalent) {
    if (area_equivalent <= 0.0 || area_baseline <= 0.0) {
        throw std::invalid_argument(
            "equivalent_blowing_ratio: areas must be positive");
    }
    return M_baseline * area_baseline / area_equivalent;
}

// -------------------------------------------------------------
// Effusion plate internal heat transfer -- Andrews 86-GT-225
// -------------------------------------------------------------

namespace {

void check_effusion_inputs(double Re, double Pr) {
    if (!(Re > 0.0)) {
        throw std::invalid_argument(
            "effusion internal Nusselt: Re must be positive");
    }
    if (!(Pr > 0.0)) {
        throw std::invalid_argument(
            "effusion internal Nusselt: Pr must be positive");
    }
}

}  // namespace

double mills_entry_length_factor(double L_over_D) {
    if (!(L_over_D > 0.0)) {
        throw std::invalid_argument(
            "mills_entry_length_factor: L/D must be positive");
    }
    namespace A = andrews1986;
    if (L_over_D <= A::mills_branch_L_over_D) {
        const double z = L_over_D;                                    // Eq. 13
        return 0.13 * z * z * z - 0.75 * z * z + 1.04 * z + 2.24;
    }
    const double w = 1.0 / L_over_D;                                  // Eq. 14
    return 1.0 - 49.3 * w * w * w * w + 58.6 * w * w * w
         - 26.5 * w * w + 7.48 * w;
}

double effusion_approach_nusselt(double Re, double Pr, double X_over_L) {
    check_effusion_inputs(Re, Pr);
    if (!(X_over_L > 0.0)) {
        throw std::invalid_argument(
            "effusion_approach_nusselt: X/L must be positive");
    }
    namespace A = andrews1986;
    return A::sparrow_coefficient * std::pow(Re, A::sparrow_re_exponent)
         * std::cbrt(Pr) * X_over_L / M_PI;                           // Eq. 18
}

double effusion_throat_nusselt(double Re, double Pr, double L_over_D) {
    check_effusion_inputs(Re, Pr);
    namespace A = andrews1986;
    return A::mills_coefficient * std::pow(Re, A::mills_re_exponent)
         * std::cbrt(Pr) * mills_entry_length_factor(L_over_D);       // Eq. 12
}

double effusion_internal_nusselt(double Re, double Pr, double X_over_L,
                                 double L_over_D) {
    return effusion_approach_nusselt(Re, Pr, X_over_L)
         + effusion_throat_nusselt(Re, Pr, L_over_D);                 // Eq. 19
}

}  // namespace combaero::cooling
