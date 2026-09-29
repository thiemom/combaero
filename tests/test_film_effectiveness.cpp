// Baldauf et al. (2002) laterally averaged film-cooling effectiveness.
//
// Ground truth is the paper's own worked example (Tables 3 and 4) and its
// Fig. 14 measurement-versus-correlation plot, not the code's output.
// See validation/cooling/extractions/baldauf_2002_film_effectiveness.md.

#include "cooling_correlations.h"
#include "correlation_status.h"
#include "math_constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

namespace cc = combaero::cooling;
namespace B = combaero::cooling::baldauf2002;

namespace {

// Table 3: the paper's explicit application example, which is also the case
// plotted in Fig. 14.
constexpr double T3_M = 2.0;
constexpr double T3_ALPHA = 30.0;
constexpr double T3_SD = 3.0;
constexpr double T3_P = 1.2;
constexpr double T3_TU = 0.015;

double eta(double x_over_D, double M = T3_M, double P = T3_P,
           double alpha = T3_ALPHA, double sD = T3_SD, double Tu = T3_TU) {
    return cc::film_effectiveness_baldauf_2002(x_over_D, M, P, alpha, sD, Tu);
}

// Same call, but takes the extrapolation flag instead of emitting a warning.
// Passing a status is also what keeps the suite's own out-of-envelope probes
// from writing to stderr.
bool extrapolates(double M, double P, double alpha, double sD, double Tu) {
    combaero::CorrelationStatus status = combaero::CorrelationStatus::Valid;
    cc::film_effectiveness_baldauf_2002(20.0, M, P, alpha, sD, Tu, &status);
    return status == combaero::CorrelationStatus::Extrapolated;
}

}  // namespace

TEST(BaldaufFilmEffectiveness, ReproducesFigure14ForThePapersOwnExample) {
    // Fig. 14 plots measurement against correlation for exactly the Table 3
    // conditions over x/D = 0 to 80. The correlation curve rises from the
    // ejection point and settles onto a plateau at eta ~ 0.11 to 0.12.
    EXPECT_LT(eta(1.0), 0.01);          // still climbing out of the ejection
    EXPECT_GT(eta(10.0), 0.05);
    for (double x : {40.0, 50.0, 60.0, 70.0, 80.0}) {
        const double e = eta(x);
        EXPECT_GT(e, 0.105) << "x/D = " << x;
        EXPECT_LT(e, 0.130) << "x/D = " << x;
    }
    // The apex sits well upstream of the plateau, which is the whole point
    // of the base curve: a rise to a peak, then a slow decay.
    EXPECT_GT(eta(30.0), eta(80.0));
}

TEST(BaldaufFilmEffectiveness, RisesFromZeroAtTheEjectionPoint) {
    // Older correlations extrapolate an exponential decay backwards and
    // produce an unphysical maximum at x = 0. This one is built to be valid
    // "from the point of the ejection", so eta must climb from ~0.
    double prev = eta(0.5);
    EXPECT_GE(prev, 0.0);
    EXPECT_LT(prev, 0.005);
    for (double x = 1.0; x <= 20.0; x += 0.5) {
        const double e = eta(x);
        EXPECT_GT(e, prev) << "not rising at x/D = " << x;
        prev = e;
    }
}

TEST(BaldaufFilmEffectiveness, PeaksThenFallsWithBlowingRate) {
    // Jet lift-off. The paper's Fig. 8 shows the peak effectiveness rising
    // at low blowing rate and falling at moderate to high, and Eq. (10) is
    // the two-branch fit that encodes it. So eta must be non-monotone in M.
    double best = -1.0, best_M = 0.0;
    for (int i = 0; i <= 60; ++i) {
        const double M = 0.2 + (2.5 - 0.2) * i / 60.0;
        const double e = eta(20.0, M);
        if (e > best) { best = e; best_M = M; }
    }
    EXPECT_GT(best_M, 0.2);
    EXPECT_LT(best_M, 2.5) << "peak must be interior, not at an endpoint";
    EXPECT_LT(eta(20.0, 2.5), best) << "lift-off must reduce eta at high M";
    EXPECT_LT(eta(20.0, 0.2), best);
}

TEST(BaldaufFilmEffectiveness, CloserHolesCoolBetter) {
    // More coolant per unit span at fixed M. Monotone in s/D across the
    // paper's range.
    double prev = eta(20.0, 1.0, T3_P, T3_ALPHA, 2.0);
    for (double sD : {2.5, 3.0, 4.0, 5.0}) {
        const double e = eta(20.0, 1.0, T3_P, T3_ALPHA, sD);
        EXPECT_LT(e, prev) << "s/D = " << sD;
        prev = e;
    }
}

TEST(BaldaufFilmEffectiveness, StaysWithinPhysicalBounds) {
    // eta is a normalised temperature difference; it cannot leave [0, 1].
    for (double M : {0.2, 0.5, 1.0, 1.8, 2.5}) {
        for (double sD : {2.0, 3.0, 5.0}) {
            for (double al : {30.0, 60.0, 90.0}) {
                for (double x : {0.5, 5.0, 25.0, 100.0}) {
                    const double e = eta(x, M, 1.5, al, sD);
                    EXPECT_TRUE(std::isfinite(e));
                    EXPECT_GE(e, 0.0) << "M=" << M << " s/D=" << sD
                                      << " a=" << al << " x/D=" << x;
                    EXPECT_LE(e, 1.0) << "M=" << M << " s/D=" << sD
                                      << " a=" << al << " x/D=" << x;
                }
            }
        }
    }
}

TEST(BaldaufFilmEffectiveness, AnalyticDerivativesAgreeWithFiniteDifferences) {
    // The (f, J) rule. One templated equation chain serves both the value
    // and the dual-number path, so this also confirms the two have not
    // drifted apart.
    for (double x : {5.0, 20.0, 60.0}) {
        for (double M : {0.4, 1.0, 2.0}) {
            for (double P : {1.2, 1.7}) {
                const auto [e, dM, dP] =
                    cc::film_effectiveness_baldauf_2002_and_derivatives(
                        x, M, P, T3_ALPHA, T3_SD, T3_TU);
                EXPECT_NEAR(e, eta(x, M, P), 1e-12);

                const double hM = M * 1e-6, hP = P * 1e-6;
                const double fdM =
                    (eta(x, M + hM, P) - eta(x, M - hM, P)) / (2.0 * hM);
                const double fdP =
                    (eta(x, M, P + hP) - eta(x, M, P - hP)) / (2.0 * hP);
                EXPECT_NEAR(dM, fdM, std::max(1e-9, std::abs(fdM) * 1e-5))
                    << "d/dM at x/D=" << x << " M=" << M << " P=" << P;
                EXPECT_NEAR(dP, fdP, std::max(1e-9, std::abs(fdP) * 1e-5))
                    << "d/dP at x/D=" << x << " M=" << M << " P=" << P;
            }
        }
    }
}

TEST(BaldaufFilmEffectiveness, DerivativeChangesSignThroughLiftOff) {
    // Not a curiosity: d eta/dM must be positive below the lift-off peak and
    // negative above it, or the solver would push blowing the wrong way.
    const auto lo = cc::film_effectiveness_baldauf_2002_and_derivatives(
        20.0, 0.4, T3_P, T3_ALPHA, T3_SD, T3_TU);
    const auto hi = cc::film_effectiveness_baldauf_2002_and_derivatives(
        20.0, 2.0, T3_P, T3_ALPHA, T3_SD, T3_TU);
    EXPECT_GT(std::get<1>(lo), 0.0);
    EXPECT_LT(std::get<1>(hi), 0.0);
}

TEST(BaldaufFilmEffectiveness, RefusesInputsTheCorrelationCannotRepresent) {
    EXPECT_THROW(eta(0.0), std::invalid_argument);       // at/upstream of ejection
    EXPECT_THROW(eta(-1.0), std::invalid_argument);
    EXPECT_THROW(eta(10.0, 0.0), std::invalid_argument); // no coolant
    EXPECT_THROW(eta(10.0, 1.0, 0.0), std::invalid_argument);
    // Eq. (34) carries exp(-0.0012/Tu^2), singular at zero turbulence.
    EXPECT_THROW(eta(10.0, 1.0, 1.2, 30.0, 3.0, 0.0), std::invalid_argument);
    // alpha is to the SURFACE, so 0 and >90 are not ejection angles.
    EXPECT_THROW(eta(10.0, 1.0, 1.2, 0.0), std::invalid_argument);
    EXPECT_THROW(eta(10.0, 1.0, 1.2, 120.0), std::invalid_argument);
}

TEST(BaldaufFilmEffectiveness, PublishedEnvelopeIsReportedNotEnforced) {
    // Still NOT enforced, and for the original reason: a network solve
    // transits odd states during Newton iteration, and refusing there would
    // break convergence rather than protect anyone. Evaluating outside the
    // envelope must keep answering, finitely.
    EXPECT_DOUBLE_EQ(B::M_min, 0.2);
    EXPECT_DOUBLE_EQ(B::M_max, 2.5);
    EXPECT_DOUBLE_EQ(B::s_over_D_min, 2.0);
    EXPECT_DOUBLE_EQ(B::s_over_D_max, 5.0);
    EXPECT_DOUBLE_EQ(B::rms_deviation, 0.055);   // the paper's own stated RMS
    EXPECT_TRUE(std::isfinite(eta(20.0, 3.0, 1.2, 30.0, 6.0, 0.10)));

    // What changed: it is now REPORTED. The constants were declared with a
    // comment saying that outside them the correlation is extrapolation, and
    // for two releases nothing read them -- so the correlation answered a
    // 7.4 hole spacing as confidently as a 3.0 one. The validation harness
    // derives its `extrapolated` column from this signal, so silence there
    // makes extrapolation indistinguishable from model error.
    EXPECT_TRUE(extrapolates(3.0, 1.2, 30.0, 6.0, 0.10));
}

TEST(BaldaufFilmEffectiveness, EachEnvelopeBoundIsCheckedIndependently) {
    // Every bound, low side and high side, one parameter at a time -- so a
    // missing check cannot hide behind a neighbouring one that fires.
    EXPECT_FALSE(extrapolates(T3_M, T3_P, T3_ALPHA, T3_SD, T3_TU))
        << "the paper's own Table 3 case must be inside its own envelope";

    struct Case { const char* what; double M, P, alpha, sD, Tu; };
    const Case outside[] = {
        {"M below",        B::M_min * 0.5,  T3_P, T3_ALPHA, T3_SD, T3_TU},
        {"M above",        B::M_max * 1.5,  T3_P, T3_ALPHA, T3_SD, T3_TU},
        {"P below",        T3_M, B::P_min * 0.5,  T3_ALPHA, T3_SD, T3_TU},
        {"P above",        T3_M, B::P_max * 1.5,  T3_ALPHA, T3_SD, T3_TU},
        {"alpha below",    T3_M, T3_P, B::alpha_deg_min * 0.5, T3_SD, T3_TU},
        {"s/D below",      T3_M, T3_P, T3_ALPHA, B::s_over_D_min * 0.5, T3_TU},
        {"s/D above",      T3_M, T3_P, T3_ALPHA, B::s_over_D_max * 1.5, T3_TU},
        {"Tu below",       T3_M, T3_P, T3_ALPHA, T3_SD, B::Tu_min * 0.5},
        {"Tu above",       T3_M, T3_P, T3_ALPHA, T3_SD, B::Tu_max * 1.5},
    };
    for (const Case& c : outside) {
        EXPECT_TRUE(extrapolates(c.M, c.P, c.alpha, c.sD, c.Tu)) << c.what;
    }

    // alpha_deg_max is 90, which check_inputs already refuses above, so the
    // high side is covered by the throw rather than by this flag.
    EXPECT_DOUBLE_EQ(B::alpha_deg_max, 90.0);

    // The bounds themselves are inclusive -- sitting exactly on a published
    // limit is inside the published range, not outside it.
    EXPECT_FALSE(extrapolates(B::M_min, B::P_min, B::alpha_deg_min,
                              B::s_over_D_min, B::Tu_min));
    EXPECT_FALSE(extrapolates(B::M_max, B::P_max, B::alpha_deg_max,
                              B::s_over_D_max, B::Tu_max));
}

TEST(BaldaufFilmEffectiveness, AndreiTwentyFourteenRigIsOutsideTheEnvelope) {
    // The rig behind validation/cooling/data/andrei2014: d = 1.5 mm holes at
    // 30 deg, spanwise pitch s/d = 7.37, blowing 1 to 3, density ratio 1.0
    // and 1.5. Recorded as a test because it decides how that dataset may be
    // read: NO combination of its conditions sits inside Baldauf's envelope,
    // so scoring against it is a cross-source ACCURACY check on an
    // extrapolated model, never a fidelity check.
    const double sD = 7.37, alpha = 30.0, Tu = 0.05;
    for (double M : {1.0, 2.0, 3.0}) {
        for (double P : {1.0, 1.5}) {
            EXPECT_TRUE(extrapolates(M, P, alpha, sD, Tu))
                << "BR = " << M << ", DR = " << P;
        }
    }
    // s/D alone is enough: 7.37 against a published maximum of 5.
    EXPECT_GT(sD, B::s_over_D_max);
}

TEST(BaldaufFilmEffectiveness, DerivativesReportExtrapolationToo) {
    // The solver-facing form shares the check, so a Newton step cannot
    // wander outside the envelope unreported while the value form would
    // have said so.
    combaero::CorrelationStatus status = combaero::CorrelationStatus::Valid;
    cc::film_effectiveness_baldauf_2002_and_derivatives(
        20.0, T3_M, T3_P, T3_ALPHA, B::s_over_D_max * 1.5, T3_TU, &status);
    EXPECT_EQ(status, combaero::CorrelationStatus::Extrapolated);

    status = combaero::CorrelationStatus::Extrapolated;
    cc::film_effectiveness_baldauf_2002_and_derivatives(
        20.0, T3_M, T3_P, T3_ALPHA, T3_SD, T3_TU, &status);
    EXPECT_EQ(status, combaero::CorrelationStatus::Valid);
}

TEST(BaldaufFilmEffectiveness, ReportingDoesNotChangeTheValue) {
    // The whole point of reporting rather than enforcing: eta outside the
    // envelope must be bit-identical to what it was before the check
    // existed, whether or not a status is requested.
    combaero::CorrelationStatus status = combaero::CorrelationStatus::Valid;
    const double with_status = cc::film_effectiveness_baldauf_2002(
        20.0, 3.0, 1.0, 30.0, 7.37, 0.05, &status);
    const double without = cc::film_effectiveness_baldauf_2002(
        20.0, 3.0, 1.0, 30.0, 7.37, 0.05, nullptr);
    EXPECT_EQ(with_status, without);
    EXPECT_EQ(status, combaero::CorrelationStatus::Extrapolated);
}

TEST(BaldaufFilmEffectiveness, PeakEffectivenessMatchesFigure8Magnitudes) {
    // Fig. 8 plots the PEAK laterally averaged effectiveness against M, one
    // panel per hole spacing, with axis bands of roughly 0.35-0.6 at s/D=2,
    // 0.15-0.45 at s/D=3 and 0.05-0.3 at s/D=5. Those bands span all nine
    // alpha/P combinations in each panel, so this is a magnitude and
    // ordering check, not a tight one -- the tight guard is the pinned
    // values below.
    auto peak = [](double M, double sD) {
        double best = 0.0;
        for (int i = 1; i <= 400; ++i) {
            best = std::max(best, cc::film_effectiveness_baldauf_2002(
                                      0.5 * i, M, 1.2, 30.0, sD, 0.015));
        }
        return best;
    };
    for (double M : {0.5, 1.0, 2.0, 2.5}) {
        // Closer spacing always cools better at fixed M.
        EXPECT_GT(peak(M, 2.0), peak(M, 3.0)) << "M = " << M;
        EXPECT_GT(peak(M, 3.0), peak(M, 5.0)) << "M = " << M;
        // And every peak sits inside its panel's plotted band, generously.
        EXPECT_GT(peak(M, 2.0), 0.25);  EXPECT_LT(peak(M, 2.0), 0.62);
        EXPECT_GT(peak(M, 3.0), 0.10);  EXPECT_LT(peak(M, 3.0), 0.47);
        EXPECT_GT(peak(M, 5.0), 0.02);  EXPECT_LT(peak(M, 5.0), 0.32);
    }
}

TEST(BaldaufFilmEffectiveness, PinnedValuesAcrossTheEnvelope) {
    // REGRESSION PINS, and labelled as such: these are this implementation's
    // own output, not values the paper prints. Their job is to catch drift
    // in any single coefficient, which the figure-based checks above are too
    // loose to see -- falsification showed that perturbing Eq. (15)'s 0.048,
    // dropping Eq. (39)'s density backscale, or collapsing Eq. (32) to
    // b_1 = b_0 all passed the physical tests untouched.
    //
    // Their literature anchor is the Fig. 14 and Fig. 8 tests above: those
    // say the curve is in the right place, these say it has not moved.
    struct Case { double x, M, P, alpha, sD, Tu, eta; };
    const Case cases[] = {
        {20.0, 2.0, 1.2, 30.0, 3.0, 0.0150, 0.121878791320},
        {60.0, 2.0, 1.2, 30.0, 3.0, 0.0150, 0.117546338033},
        {10.0, 0.5, 1.2, 30.0, 3.0, 0.0150, 0.279616731336},
        {30.0, 1.0, 1.8, 60.0, 2.0, 0.0400, 0.273521417342},
        {15.0, 1.5, 1.5, 90.0, 5.0, 0.0035, 0.069395466122},
        {50.0, 2.5, 1.2, 45.0, 2.5, 0.0750, 0.214573682614},
        { 5.0, 0.2, 1.8, 30.0, 4.0, 0.0200, 0.185276331495},
        {80.0, 1.2, 1.4, 75.0, 3.5, 0.0100, 0.101217046875},
    };
    for (const auto& c : cases) {
        EXPECT_NEAR(cc::film_effectiveness_baldauf_2002(c.x, c.M, c.P, c.alpha,
                                                        c.sD, c.Tu),
                    c.eta, 1e-9)
            << "x/D=" << c.x << " M=" << c.M << " P=" << c.P
            << " alpha=" << c.alpha << " s/D=" << c.sD << " Tu=" << c.Tu;
    }
}

// -------------------------------------------------------------------------
// Transcription pins
// -------------------------------------------------------------------------
//
// These recompute the paper's coefficient functions from the equations AS
// TRANSCRIBED and compare against the numbers Table 4 prints. They pin the
// TRANSCRIPTION -- they do not exercise the library -- and they are here so
// that the one place the paper contradicts itself cannot be quietly
// rediscovered or quietly "fixed" later.

TEST(BaldaufTable4, SeventeenCoefficientsReproduceThePapersWorkedExample) {
    const double sD = T3_SD, M = T3_M, P = T3_P, Tu = T3_TU;
    const double a = T3_ALPHA * M_PI / 180.0;
    const double U = M / P;

    EXPECT_NEAR(0.6 + 0.4 * (2 - std::cos(a)) / (1 + std::pow((sD - 1) / 3.3, 6.0)),
                1.0321731, 1e-7);                                  // xi_c,  Eq. 18
    EXPECT_NEAR(0.465 / (1 + 0.048 * sD * sD), 0.32472067, 1e-8);  // eta_c0, Eq. 15
    EXPECT_NEAR(U * std::pow(P, 0.8) * (1 - (0.03 + 0.11 * (5 - sD)) * std::cos(a)),
                1.5108774, 1e-7);                                  // mu,     Eq. 9
    EXPECT_NEAR(std::exp(1.92 - 7.5 * std::pow(sD, -1.5)), 1.6106283, 1e-7);  // b, Eq. 12
    EXPECT_NEAR(0.125 + 0.063 * std::pow(sD, 1.8), 0.58015447, 1e-8);        // mu_0, Eq. 14
    EXPECT_NEAR(0.7 + 336 * std::exp(-1.85 * sD), 2.0061856, 1e-7);          // c, Eq. 13

    const double xi_hat = 1.17 * (1 - (sD - 1) / (1 + 0.2 * (sD - 1) * (sD - 1)))
                        * (std::cos(2.3 * a) + 2.45);
    EXPECT_NEAR(xi_hat, -0.36508783, 1e-8);                                  // Eq. 26
    const double eta_hat = 0.022 * (sD + 1) * (0.9 - std::sin(2 * a))
                         - (0.08 + 0.46 / (1 + (sD - 3.2) * (sD - 3.2)));
    EXPECT_NEAR(eta_hat, -0.51931793, 1e-8);                                 // Eq. 27
    const double g = 0.75 * (1 - std::exp(-0.8 * (sD - 1)));
    EXPECT_NEAR(g, 0.59857761, 1e-8);                                        // Eq. 24
    const double k = 2 * (1 - std::exp(0.57 * (1 - sD)))
                   + 0.91 * std::pow(std::cos(a), 0.65);
    EXPECT_NEAR(k, 2.1891363, 1e-7);                                         // Eq. 25
    const double tr = 1.0 / (1 + std::pow(U * std::pow(P, g) / k, -5.0));
    EXPECT_NEAR(1 + xi_hat * tr, 0.88819413, 1e-8);                          // Eq. 22
    EXPECT_NEAR(1 + eta_hat * tr, 0.84096213, 1e-8);                         // Eq. 23

    const double a_1 = 0.04 + 0.23 * sD + (0.95 - 0.19 * sD) * std::cos(1.5 * a);
    EXPECT_NEAR(a_1, 0.99870058, 1e-8);                                      // Eq. 30
    EXPECT_NEAR(65 / std::pow(M / 2.5, a_1), 81.226444, 1e-6);               // Eq. 29
    EXPECT_NEAR(7.5 + sD, 10.5, 1e-12);                                      // Eq. 33

    const double b_T = 0.7 * (1 + (1.22 / (1 + 7 * std::pow(sD - 1, -7.0))
                                   + 0.87 + std::cos(2.5 * a))
                              * std::exp(2.6 * Tu - 0.0012 / (Tu * Tu) - 1.76));
    EXPECT_NEAR(b_T, 0.70138176, 1e-8);                                      // Eq. 34
    EXPECT_NEAR(2.5 * std::pow(5.8 / 2.5, b_T / 0.7), 5.809643, 1e-6);       // Eq. 35
}

TEST(BaldaufTable4, EquationThirtyOneContradictsThePapersOwnTable) {
    // THE ONE DISCREPANCY, pinned deliberately.
    //
    // Eq. (31) as printed gives 0.83612467 where Table 4 prints 0.61626073,
    // a 36% gap, while all seventeen other coefficients above land within
    // 2e-6. The printed equation was confirmed from four independent
    // channels, no angle convention or unit choice reaches Table 4's value,
    // and Table 4 is self-consistent with Eq. (32) to the digit.
    //
    // We implement Eq. (31) AS PRINTED. This test exists so that a future
    // reader meets the contradiction immediately instead of rediscovering
    // it, and so that anyone who "fixes" the constant has to delete an
    // assertion that says why not to.
    const double sD = T3_SD;
    const double a = T3_ALPHA * M_PI / 180.0;
    const double b_0_printed =
        0.8 - 0.014 * sD * sD
        + (1.5 - 2.0 / std::sqrt(sD))
          * std::sin(0.86 * a * (1 + 0.754 / (1 + 0.87 * sD * sD)));

    EXPECT_NEAR(b_0_printed, 0.83612467, 1e-8);
    EXPECT_NEAR(B::b0_table4_at_table3, 0.61626073, 1e-12);
    EXPECT_GT(std::abs(b_0_printed - B::b0_table4_at_table3)
              / B::b0_table4_at_table3, 0.30);

    // Table 4's own b_0 and b_1 are mutually consistent with Eq. (32),
    // which is what rules out a typo in the table's b_1 alone.
    EXPECT_NEAR(0.54778731 / B::b0_table4_at_table3,
                1.0 / (1.0 + std::pow(T3_M, -3.0)), 1e-7);
}
