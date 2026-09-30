// Effusion plate internal heat transfer -- Andrews 86-GT-225.
//
// Ground truth is the paper's own algebra and its own worked geometry, not
// this code's output. See
// validation/cooling/extractions/andrews_effusion_internal_h.md.

#include "cooling_correlations.h"
#include "math_constants.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

namespace cc = combaero::cooling;
namespace A = combaero::cooling::andrews1986;

namespace {

// Andrews' own comparison geometry for Eq. (19): X = 6.11 mm, D = 0.64 mm,
// L = 6.35 mm.
constexpr double EX_X = 6.11e-3;
constexpr double EX_L = 6.35e-3;

}  // namespace

TEST(EffusionInternal, MillsBranchesAgreeWhereTheyMeet) {
    // Eqs. (13) and (14) are separate curve fits joined at L/D = 2. That
    // they meet is not something the paper asserts -- it is the check that
    // caught whether the coefficients had been transcribed correctly from a
    // poor scan, and both branches landing on 2.36 from different
    // polynomials is not a coincidence a misreading would reproduce.
    const double below = cc::mills_entry_length_factor(2.0 - 1e-12);
    const double above = cc::mills_entry_length_factor(2.0 + 1e-12);
    EXPECT_NEAR(below, above, 2e-3);
    EXPECT_NEAR(below, 2.36, 5e-3);
}

TEST(EffusionInternal, EntryLengthFactorDecaysToFullyDeveloped) {
    // The physical meaning of R_Nu: a short hole never develops, so it is
    // enhanced; a long one must recover the plain Dittus-Boelter value.
    EXPECT_GT(cc::mills_entry_length_factor(0.5), 2.0);
    EXPECT_GT(cc::mills_entry_length_factor(2.0), 2.0);

    double previous = cc::mills_entry_length_factor(3.0);
    for (double LD : {5.0, 10.0, 25.0, 100.0, 500.0}) {
        const double v = cc::mills_entry_length_factor(LD);
        EXPECT_LT(v, previous) << "L/D = " << LD;
        EXPECT_GT(v, 1.0);
        previous = v;
    }
    EXPECT_NEAR(cc::mills_entry_length_factor(5000.0), 1.0, 2e-3);

    EXPECT_THROW(cc::mills_entry_length_factor(0.0), std::invalid_argument);
    EXPECT_THROW(cc::mills_entry_length_factor(-1.0), std::invalid_argument);
}

TEST(EffusionInternal, ReproducesEquationNineteensPrintedCoefficient) {
    // Eq. (19) prints the summed form for one geometry:
    //
    //   Nu = (0.27 Re^0.476 + 0.023 Re^0.8 R_Nu) Pr^(1/3)
    //
    // 0.27 is the ONLY place the paper evaluates Eq. (18)'s geometry factor
    // numerically, so it is the one arithmetic check available on the
    // approach term. 0.881 * X/(pi L) with X = 6.11 mm and L = 6.35 mm.
    const double X_over_L = EX_X / EX_L;
    const double implied =
        A::sparrow_coefficient * X_over_L / M_PI;
    EXPECT_NEAR(implied, 0.27, 5e-4);

    // And the function must carry that same factor. Compared RELATIVELY,
    // because 0.881 X/(pi L) is 0.26983 and the paper prints it rounded to
    // 0.27 -- a 0.06% difference that an absolute tolerance on a Nusselt
    // number of order 14 would flag as an error when it is the source's own
    // rounding.
    const double Re = 5000.0, Pr = 0.72;
    const double printed = 0.27 * std::pow(Re, 0.476) * std::cbrt(Pr);
    EXPECT_NEAR(cc::effusion_approach_nusselt(Re, Pr, X_over_L) / printed,
                1.0, 1e-3);
}

TEST(EffusionInternal, TheSumIsTheTwoTermsAndNothingElse) {
    // Eq. (19) is a summation -- "the authors have treated the wall heat
    // transfer as the summation of Equations 13 or 14 and 15" -- not a
    // blend, not a maximum. Both terms are on the hole diameter and hole
    // internal area, which is what makes adding them legitimate.
    const double Re = 4000.0, Pr = 0.72, X_over_L = 2.42, L_over_D = 1.93;
    EXPECT_DOUBLE_EQ(cc::effusion_internal_nusselt(Re, Pr, X_over_L, L_over_D),
                     cc::effusion_approach_nusselt(Re, Pr, X_over_L)
                         + cc::effusion_throat_nusselt(Re, Pr, L_over_D));
}

TEST(EffusionInternal, TheApproachTermDominatesAtLowReynolds) {
    // The paper's headline: "the hole approach flow heat transfer is much
    // larger than the internal hole heat transfer". It is a claim about
    // exponents -- 0.476 against 0.8 -- so the approach must lead at low Re
    // and be overtaken as Re grows. Andrews' plate C geometry.
    const double X_over_L = 15.24 / 6.3, L_over_D = 6.3 / 3.27, Pr = 0.727;

    const double lo = 500.0;
    EXPECT_GT(cc::effusion_approach_nusselt(lo, Pr, X_over_L),
              cc::effusion_throat_nusselt(lo, Pr, L_over_D));

    const double hi = 50000.0;
    EXPECT_LT(cc::effusion_approach_nusselt(hi, Pr, X_over_L),
              cc::effusion_throat_nusselt(hi, Pr, L_over_D));

    // The crossover sits inside the rig's own range, which is why both
    // terms are needed rather than whichever is nominally "dominant".
    double crossover = 0.0;
    for (double Re = 100.0; Re < 100000.0; Re *= 1.02) {
        if (cc::effusion_throat_nusselt(Re, Pr, L_over_D)
            > cc::effusion_approach_nusselt(Re, Pr, X_over_L)) {
            crossover = Re;
            break;
        }
    }
    EXPECT_GT(crossover, 1000.0);
    EXPECT_LT(crossover, 10000.0);
}

TEST(EffusionInternal, PlateAreaConversionIsNotOptional) {
    // Nu is on the hole internal area; an effusion element wants a
    // coefficient per unit PLATE area. For Andrews' plate C the two differ
    // by A/A_h = 3.46, so leaving the conversion out is not a refinement.
    const double D = 3.27e-3, X = 15.24e-3, L = 6.3e-3;
    const double A_plate = X * X - M_PI * D * D / 4.0;
    const double A_hole = M_PI * D * L;
    EXPECT_NEAR(A_plate / A_hole, 3.4, 0.1);   // Table 1 prints 3.4

    // X is 0.6 inch; the table's rounded 15.2 mm does not reproduce its own
    // hole count N = 4306 per square metre.
    EXPECT_NEAR(1.0 / (X * X), 4306.0, 1.0);
    EXPECT_GT(std::abs(1.0 / (15.2e-3 * 15.2e-3) - 4306.0), 20.0);
}

TEST(EffusionInternal, RefusesInputsTheCorrelationsCannotRepresent) {
    EXPECT_THROW(cc::effusion_approach_nusselt(0.0, 0.72, 2.4),
                 std::invalid_argument);
    EXPECT_THROW(cc::effusion_approach_nusselt(1000.0, 0.0, 2.4),
                 std::invalid_argument);
    EXPECT_THROW(cc::effusion_approach_nusselt(1000.0, 0.72, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(cc::effusion_throat_nusselt(-1.0, 0.72, 2.0),
                 std::invalid_argument);
    EXPECT_THROW(cc::effusion_throat_nusselt(1000.0, 0.72, -2.0),
                 std::invalid_argument);
}
