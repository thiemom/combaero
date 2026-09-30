// Multi-row film superposition: Sellers, and Gao et al. (2025) Eq. (7).
//
// Ground truth is the published algebra -- Gao Eqs. (1), (5), (7), (9), (10)
// -- and the physical claim the correction exists to fix: Sellers
// overestimates, and worsens as rows accumulate.
// See validation/cooling/extractions/film_superposition.md.

#include "cooling_correlations.h"
#include "math_constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace cc = combaero::cooling;

namespace {

// Gao Eq. (1) written out as the SUM the paper prints, rather than the
// product form the implementation uses. Two algebraically identical routes,
// so agreement checks the implementation against the paper's own statement.
double sellers_sum_form(const std::vector<double>& e) {
    double total = e.at(0);
    for (std::size_t i = 1; i < e.size(); ++i) {
        double prod = 1.0;
        for (std::size_t j = 0; j < i; ++j) {
            prod *= (1.0 - e[j]);
        }
        total += e[i] * prod;
    }
    return total;
}

}  // namespace

TEST(FilmSuperposition, ProductFormMatchesThePapersSumForm) {
    // Gao Eq. (1) is printed as a sum; 1 - prod(1 - eta_i) is the same
    // thing. The implementation uses the product because it is O(n) and
    // does not accumulate pairwise rounding.
    for (const std::vector<double>& e :
         {std::vector<double>{0.3},
          std::vector<double>{0.3, 0.22},
          std::vector<double>{0.30, 0.22, 0.18, 0.12},
          std::vector<double>{0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05}}) {
        EXPECT_NEAR(cc::film_superposition_sellers(e), sellers_sum_form(e), 1e-14);
    }
}

TEST(FilmSuperposition, CorrectedFormReducesExactlyToSellersAtAlphaOne) {
    // Eq. (7) must collapse onto Eq. (1) when nothing is corrected. If it
    // does not, the indexing of the two products is wrong.
    const std::vector<double> e{0.30, 0.22, 0.18, 0.12};
    const std::vector<double> ones(e.size() - 1, 1.0);
    EXPECT_NEAR(cc::film_superposition_corrected(e, ones),
                cc::film_superposition_sellers(e), 1e-14);

    // And for a single row there is nothing to correct at all.
    EXPECT_NEAR(cc::film_superposition_corrected({0.37}, {}), 0.37, 1e-14);
}

TEST(FilmSuperposition, CorrectionAlwaysReducesTheAccumulatedEffectiveness) {
    // The whole point. Gao: "the Sellers method accumulates prediction
    // errors as the number of hole rows increases, leading to an
    // OVERESTIMATION of the cooling efficiency." So alpha < 1 must pull the
    // total down, monotonically.
    const std::vector<double> e{0.30, 0.22, 0.18, 0.12};
    double prev = cc::film_superposition_sellers(e);
    for (double a : {0.98, 0.95, 0.90, 0.80, 0.60}) {
        const std::vector<double> al(e.size() - 1, a);
        const double v = cc::film_superposition_corrected(e, al);
        EXPECT_LT(v, prev) << "alpha = " << a;
        EXPECT_GT(v, 0.0);
        prev = v;
    }
}

TEST(FilmSuperposition, SellersOverestimateGrowsWithRowCount) {
    // The reason a film module built for a few rows cannot simply be reused
    // for effusion. With identical rows, Sellers marches towards 1 while the
    // corrected form is held back, and the GAP widens with every row.
    const double eta_row = 0.08, alpha = 0.9;
    double prev_gap = 0.0;
    for (std::size_t n : {2u, 4u, 8u, 16u, 30u}) {
        const std::vector<double> e(n, eta_row);
        const std::vector<double> al(n - 1, alpha);
        const double gap = cc::film_superposition_sellers(e)
                         - cc::film_superposition_corrected(e, al);
        EXPECT_GT(gap, prev_gap) << "n = " << n;
        prev_gap = gap;
    }
    // At 30 rows the uncorrected model is close to saturation, which is
    // exactly the regime effusion operates in.
    const std::vector<double> e30(30, eta_row);
    EXPECT_GT(cc::film_superposition_sellers(e30), 0.90);
    EXPECT_LT(cc::film_superposition_corrected(e30, std::vector<double>(29, alpha)),
              0.60);
}

TEST(FilmSuperposition, StaysWithinPhysicalBoundsForManyRows) {
    for (std::size_t n : {1u, 2u, 5u, 20u, 50u}) {
        for (double eta_row : {0.0, 0.05, 0.3, 0.9, 1.0}) {
            const std::vector<double> e(n, eta_row);
            const double s = cc::film_superposition_sellers(e);
            EXPECT_GE(s, 0.0);
            EXPECT_LE(s, 1.0);
            if (n > 1) {
                const double c =
                    cc::film_superposition_corrected(e, std::vector<double>(n - 1, 0.85));
                EXPECT_GE(c, 0.0);
                EXPECT_LE(c, 1.0);
                EXPECT_LE(c, s) << "the correction may only reduce";
            }
        }
    }
}

TEST(FilmSuperposition, GradientsAgreeWithFiniteDifferences) {
    const std::vector<double> e{0.30, 0.22, 0.18, 0.12, 0.07};
    const std::vector<double> al{0.95, 0.92, 0.90, 0.88};

    const auto [vs, gs] = cc::film_superposition_sellers_and_gradient(e);
    const auto [vc, gc] = cc::film_superposition_corrected_and_gradient(e, al);
    EXPECT_NEAR(vs, cc::film_superposition_sellers(e), 1e-14);
    EXPECT_NEAR(vc, cc::film_superposition_corrected(e, al), 1e-14);

    const double h = 1e-7;
    for (std::size_t i = 0; i < e.size(); ++i) {
        std::vector<double> up = e, dn = e;
        up[i] += h;
        dn[i] -= h;
        EXPECT_NEAR(gs[i],
                    (cc::film_superposition_sellers(up)
                     - cc::film_superposition_sellers(dn)) / (2.0 * h),
                    1e-6) << "Sellers d/d eta_" << i;
        EXPECT_NEAR(gc[i],
                    (cc::film_superposition_corrected(up, al)
                     - cc::film_superposition_corrected(dn, al)) / (2.0 * h),
                    1e-6) << "corrected d/d eta_" << i;
    }
}

TEST(FilmSuperposition, GradientSurvivesAFullyEffectiveRow) {
    // eta = 1 for some row makes prod(1 - eta) zero, so a gradient computed
    // by dividing the total product would be 0/0. Both gradients are built
    // from explicit partial products instead.
    const std::vector<double> e{0.3, 1.0, 0.2};
    const auto [vs, gs] = cc::film_superposition_sellers_and_gradient(e);
    EXPECT_NEAR(vs, 1.0, 1e-14);
    for (double g : gs) {
        EXPECT_TRUE(std::isfinite(g));
    }
    // Only the saturated row still moves the total; the others are masked
    // by it, which is physically right.
    EXPECT_NEAR(gs[0], 0.0, 1e-14);
    EXPECT_NEAR(gs[2], 0.0, 1e-14);
    EXPECT_GT(gs[1], 0.0);

    const auto [vc, gc] = cc::film_superposition_corrected_and_gradient(e, {0.9, 0.9});
    EXPECT_TRUE(std::isfinite(vc));
    for (double g : gc) {
        EXPECT_TRUE(std::isfinite(g));
    }
}

TEST(FilmSuperposition, CorrectionFactorFollowsTheSaturatingFormOfEquationFive) {
    // alpha = a r/(a r + 1) + b. At no coolant it is b; it rises
    // monotonically and saturates at 1 + b.
    const double a = 8.0, b = 0.1;
    EXPECT_NEAR(cc::mainstream_temperature_correction(0.0, a, b), b, 1e-14);
    double prev = b;
    for (double r : {0.001, 0.01, 0.05, 0.2, 1.0, 100.0}) {
        const double v = cc::mainstream_temperature_correction(r, a, b);
        EXPECT_GT(v, prev) << "r = " << r;
        prev = v;
    }
    EXPECT_NEAR(cc::mainstream_temperature_correction(1e9, a, b), 1.0 + b, 1e-6);
    EXPECT_THROW(cc::mainstream_temperature_correction(-0.1, a, b),
                 std::invalid_argument);
}

TEST(FilmSuperposition, EquivalentSlotAndBlowingRatio) {
    // Eq. (9): s = A_hole / pitch. A 1 mm hole on a 3 mm pitch.
    const double d = 1.0e-3, pitch = 3.0e-3;
    const double area = M_PI * d * d / 4.0;
    EXPECT_NEAR(cc::equivalent_slot_width(area, pitch), area / pitch, 1e-18);
    // Closer holes give a wider equivalent slot at fixed hole size.
    EXPECT_GT(cc::equivalent_slot_width(area, 2.0e-3),
              cc::equivalent_slot_width(area, 5.0e-3));

    // Eq. (10): M_e = M_0 A_0 / A_e. Spreading the same blowing over more
    // area lowers the equivalent blowing ratio.
    EXPECT_NEAR(cc::equivalent_blowing_ratio(1.5, 2.0e-6, 4.0e-6), 0.75, 1e-14);
    EXPECT_NEAR(cc::equivalent_blowing_ratio(1.5, 2.0e-6, 2.0e-6), 1.5, 1e-14);

    EXPECT_THROW(cc::equivalent_slot_width(0.0, pitch), std::invalid_argument);
    EXPECT_THROW(cc::equivalent_blowing_ratio(1.0, 1.0, 0.0), std::invalid_argument);
}

TEST(FilmSuperposition, RefusesMalformedInput) {
    EXPECT_THROW(cc::film_superposition_sellers({}), std::invalid_argument);
    EXPECT_THROW(cc::film_superposition_sellers({0.3, 1.2}), std::invalid_argument);
    EXPECT_THROW(cc::film_superposition_sellers({0.3, -0.1}), std::invalid_argument);

    // alpha must have exactly one fewer entry than eta -- it lives BETWEEN
    // rows, and an off-by-one here would silently shift every correction.
    EXPECT_THROW(cc::film_superposition_corrected({0.3, 0.2}, {}),
                 std::invalid_argument);
    EXPECT_THROW(cc::film_superposition_corrected({0.3, 0.2}, {0.9, 0.9}),
                 std::invalid_argument);
    // alpha above 1 would create coolant.
    EXPECT_THROW(cc::film_superposition_corrected({0.3, 0.2}, {1.4}),
                 std::invalid_argument);
    EXPECT_THROW(cc::film_superposition_corrected({0.3, 0.2}, {-0.1}),
                 std::invalid_argument);
}

TEST(FilmSuperposition, GaosPublishedCoefficientsAreRecordedButNotUsable) {
    // Gao section 4.3.1: "The empirical coefficients in Equation (5), a and
    // b, were determined to be 12 and 0.9465, respectively." An earlier note
    // in this repo claimed they were never published -- a grep looked for
    // "a = 12" while the paper states the value in prose.
    namespace FS = combaero::cooling::film_superposition;
    EXPECT_DOUBLE_EQ(FS::gao_a_case1, 12.0);
    EXPECT_DOUBLE_EQ(FS::gao_b_case1, 0.9465);

    // Why they are recorded rather than defaulted. With b = 0.9465 the
    // whole reachable range of alpha is [b, 1 + b), so the correction
    // saturates at 5.35% per row however much coolant is added.
    EXPECT_NEAR(cc::mainstream_temperature_correction(0.0, FS::gao_a_case1,
                                                      FS::gao_b_case1),
                FS::gao_b_case1, 1e-12);
    const double reachable = 1.0 - FS::gao_b_case1;
    EXPECT_LT(reachable, 0.06) << "the published pair can barely correct";

    // The floor is what rules them out, and it needs no rig arithmetic:
    // a r/(a r + 1) >= 0, so alpha >= b for EVERY r. Murray's plate needs
    // about 0.85 at low blowing and 0.69 near M = 1, both below b.
    //
    // Asserted this way ON PURPOSE. An r-dependent claim -- "alpha exceeds
    // 1 above M = 0.20 on that rig" -- was written first and withdrawn:
    // Gao's test section dimensions are not stated in the paper, so r
    // cannot be put on their scale and the claim was unverifiable.
    for (double r : {0.0, 1e-6, 1e-4, 0.001, 0.01, 0.1, 1.0, 100.0}) {
        EXPECT_GE(cc::mainstream_temperature_correction(r, FS::gao_a_case1,
                                                        FS::gao_b_case1),
                  FS::gao_b_case1) << "r = " << r;
    }
    EXPECT_LT(0.85, FS::gao_b_case1) << "Murray's low-blowing alpha is reachable";
    EXPECT_LT(0.69, FS::gao_b_case1) << "Murray's high-blowing alpha is reachable";

    // Large r takes alpha past 1, which the superposition must refuse.
    const double big = cc::mainstream_temperature_correction(
        1.0, FS::gao_a_case1, FS::gao_b_case1);
    EXPECT_GT(big, 1.0);
    EXPECT_THROW(cc::film_superposition_corrected({0.3, 0.2}, {big}),
                 std::invalid_argument);
}

TEST(FilmSuperposition, CouplingFormWithCOneIsExactlySellers) {
    // The "third category" correction from Gao's own survey:
    // eta = eta1 + eta2 - C eta1 eta2, applied recursively. C = 1 must be
    // Sellers exactly, the same identity property alpha = 1 has -- that is
    // what makes it a correction rather than a free curve.
    //
    // Recorded here because it fits Murray better than alpha does and is
    // parameterised on blowing ratio, the axis that data varies. Not
    // implemented as an API: see
    // validation/cooling/extractions/gao_alpha_fit_on_murray.md.
    auto coupled = [](const std::vector<double>& e, double C) {
        double total = e.at(0);
        for (std::size_t i = 1; i < e.size(); ++i) {
            total = total + e[i] - C * total * e[i];
        }
        return total;
    };
    for (const std::vector<double>& e :
         {std::vector<double>{0.3, 0.2},
          std::vector<double>{0.30, 0.22, 0.18, 0.12},
          std::vector<double>(10, 0.08)}) {
        EXPECT_NEAR(coupled(e, 1.0), cc::film_superposition_sellers(e), 1e-12);
    }
    // C > 1 damps, which is the direction the data requires.
    const std::vector<double> rows(10, 0.08);
    EXPECT_LT(coupled(rows, 2.15), cc::film_superposition_sellers(rows));
    EXPECT_LT(coupled(rows, 4.29), coupled(rows, 2.15));
}
