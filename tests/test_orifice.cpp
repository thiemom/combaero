#include <gtest/gtest.h>
#include "../include/orifice.h"
#include <cmath>

// Test fixture for orifice tests
class OrificeTest : public ::testing::Test {
protected:
    OrificeGeometry geom;
    OrificeState state;

    void SetUp() override {
        // Standard test case: 50mm orifice in 100mm pipe
        geom.d = 0.050;  // 50 mm
        geom.D = 0.100;  // 100 mm

        // Typical flow conditions
        state.Re_D = 100000.0;
        state.dP = 10000.0;   // 10 kPa
        state.rho = 1.2;      // kg/m³ (air at ~20°C)
        state.mu = 1.8e-5;    // Pa·s
    }
};

// -------------------------------------------------------------
// Geometry tests
// -------------------------------------------------------------

TEST_F(OrificeTest, GeometryBeta) {
    EXPECT_DOUBLE_EQ(geom.beta(), 0.5);
}

TEST_F(OrificeTest, GeometryArea) {
    const double expected = 3.14159265358979323846 * 0.050 * 0.050 / 4.0;
    EXPECT_NEAR(geom.area(), expected, 1e-10);
}

TEST_F(OrificeTest, GeometryValid) {
    EXPECT_TRUE(geom.is_valid());

    OrificeGeometry invalid;
    invalid.d = 0.0;
    invalid.D = 0.1;
    EXPECT_FALSE(invalid.is_valid());

    invalid.d = 0.15;  // d > D
    invalid.D = 0.1;
    EXPECT_FALSE(invalid.is_valid());
}

// -------------------------------------------------------------
// Sharp thin-plate Cd tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CdSharpThinPlateRange) {
    // Cd should be in reasonable range for sharp-edged orifice
    double Cd = Cd_sharp_thin_plate(geom, state);
    EXPECT_GT(Cd, 0.58);
    EXPECT_LT(Cd, 0.68);
}

TEST_F(OrificeTest, CdSharpThinPlateBetaDependence) {
    // Cd should increase slightly with beta
    OrificeGeometry geom_low, geom_high;
    geom_low.d = 0.020;
    geom_low.D = 0.100;  // beta = 0.2
    geom_high.d = 0.070;
    geom_high.D = 0.100; // beta = 0.7

    double Cd_low = Cd_sharp_thin_plate(geom_low, state);
    double Cd_high = Cd_sharp_thin_plate(geom_high, state);

    // Both should be in valid range
    EXPECT_GT(Cd_low, 0.58);
    EXPECT_LT(Cd_low, 0.68);
    EXPECT_GT(Cd_high, 0.58);
    EXPECT_LT(Cd_high, 0.68);
}

TEST_F(OrificeTest, CdSharpThinPlateReynoldsDependence) {
    // Cd should approach asymptotic value at high Re
    OrificeState state_low, state_high;
    state_low.Re_D = 10000.0;
    state_high.Re_D = 1000000.0;

    double Cd_low = Cd_sharp_thin_plate(geom, state_low);
    double Cd_high = Cd_sharp_thin_plate(geom, state_high);

    // Both should be reasonable
    EXPECT_GT(Cd_low, 0.58);
    EXPECT_GT(Cd_high, 0.58);

    // Difference should be small at high Re
    EXPECT_LT(std::abs(Cd_high - Cd_low), 0.02);
}

// -------------------------------------------------------------
// Thick-plate Cd tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CdThickPlateCorrection) {
    // Thick plate should have higher Cd than thin plate
    geom.t = 0.010;  // 10 mm thickness, t/d = 0.2

    double Cd_thin = Cd_sharp_thin_plate(geom, state);
    double Cd_thick = Cd_thick_plate(geom, state);

    EXPECT_GT(Cd_thick, Cd_thin);
    EXPECT_LT(Cd_thick, Cd_thin * 1.3);  // Correction should be bounded
}

TEST_F(OrificeTest, CdThickPlateVeryThin) {
    // Very thin plate should have negligible correction
    geom.t = 0.0005;  // 0.5 mm, t/d = 0.01

    double Cd_thin = Cd_sharp_thin_plate(geom, state);
    double Cd_thick = Cd_thick_plate(geom, state);

    EXPECT_NEAR(Cd_thick, Cd_thin, 0.001);
}

// -------------------------------------------------------------
// Rounded-entry Cd tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CdRoundedEntryHigher) {
    // Rounded entry should have higher Cd than sharp edge
    geom.r = 0.010;  // 10 mm radius, r/d = 0.2

    double Cd_sharp = Cd_sharp_thin_plate(geom, state);
    double Cd_round = Cd_rounded_entry(geom, state);

    EXPECT_GT(Cd_round, Cd_sharp);
    EXPECT_LT(Cd_round, 1.0);  // Cannot exceed 1.0
}

TEST_F(OrificeTest, CdRoundedEntryWellRounded) {
    // Well-rounded entry (r/d >= 0.15) should approach ~0.98
    geom.r = 0.010;  // r/d = 0.2

    double Cd = Cd_rounded_entry(geom, state);
    EXPECT_GT(Cd, 0.95);
    EXPECT_LT(Cd, 1.0);
}

// -------------------------------------------------------------
// Auto-select Cd tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CdAutoSelectThinPlate) {
    // Default geometry should use thin-plate correlation
    double Cd_auto = Cd(geom, state);
    double Cd_thin = Cd_sharp_thin_plate(geom, state);
    EXPECT_DOUBLE_EQ(Cd_auto, Cd_thin);
}

TEST_F(OrificeTest, CdAutoSelectThickPlate) {
    geom.t = 0.010;  // t/d = 0.2
    double Cd_auto = Cd(geom, state);
    double Cd_thick = Cd_thick_plate(geom, state);
    EXPECT_DOUBLE_EQ(Cd_auto, Cd_thick);
}

TEST_F(OrificeTest, CdAutoSelectRounded) {
    geom.r = 0.005;  // r/d = 0.1
    double Cd_auto = Cd(geom, state);
    double Cd_round = Cd_rounded_entry(geom, state);
    EXPECT_DOUBLE_EQ(Cd_auto, Cd_round);
}

// -------------------------------------------------------------
// Flow calculation tests
// -------------------------------------------------------------

TEST_F(OrificeTest, OrificeMdot) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);

    // Expected: Cd * E * A * sqrt(2 * rho * dP)
    const double beta = geom.beta();
    const double E = 1.0 / std::sqrt(1.0 - std::pow(beta, 4.0));
    double expected = Cd_val * E * geom.area() * std::sqrt(2.0 * state.rho * state.dP);
    EXPECT_NEAR(mdot, expected, 1e-10);
}

TEST_F(OrificeTest, OrificeDpRoundTrip) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);
    double dP_calc = orifice_dP(geom, Cd_val, mdot, state.rho);
    EXPECT_NEAR(dP_calc, state.dP, 1e-6);
}

TEST_F(OrificeTest, OrificeCdFromMeasurement) {
    double Cd_val = 0.61;
    double mdot = orifice_mdot(geom, Cd_val, state.dP, state.rho);
    double Cd_calc = orifice_Cd_from_measurement(geom, mdot, state.dP, state.rho);
    EXPECT_NEAR(Cd_calc, Cd_val, 1e-10);
}

// -------------------------------------------------------------
// Correlation class tests
// -------------------------------------------------------------

TEST_F(OrificeTest, CorrelationFactory) {
    auto corr = make_correlation(CdCorrelation::ReaderHarrisGallagher);
    ASSERT_NE(corr, nullptr);
    EXPECT_FALSE(corr->name().empty());

    double Cd_class = corr->Cd(geom, state);
    double Cd_func = Cd_sharp_thin_plate(geom, state);
    EXPECT_DOUBLE_EQ(Cd_class, Cd_func);
}

TEST_F(OrificeTest, UserCorrelation) {
    auto corr = make_user_correlation(
        [](const OrificeGeometry&, const OrificeState&) { return 0.65; },
        "TestConstant");

    ASSERT_NE(corr, nullptr);
    EXPECT_EQ(corr->name(), "TestConstant");
    EXPECT_DOUBLE_EQ(corr->Cd(geom, state), 0.65);
}

// -------------------------------------------------------------
// Namespace function tests
// -------------------------------------------------------------

TEST_F(OrificeTest, KFromCd) {
    double Cd_val = 0.61;
    double beta = geom.beta();
    double K = orifice::K_from_Cd(Cd_val, beta);

    // K should be positive for Cd < 1
    EXPECT_GT(K, 0.0);

    // Round-trip
    double Cd_back = orifice::Cd_from_K(K, beta);
    EXPECT_NEAR(Cd_back, Cd_val, 1e-10);
}

TEST_F(OrificeTest, ThicknessCorrection) {
    const double Re_d = 1e5;  // Typical Reynolds number

    // No correction for thin plate
    EXPECT_DOUBLE_EQ(orifice::thickness_correction(0.01, 0.5, Re_d), 1.0);

    // Small thickness: reattachment benefit
    double corr_small = orifice::thickness_correction(0.2, 0.5, Re_d);
    EXPECT_GT(corr_small, 1.0);

    // Peak around t/d ~ 0.3 (calibrated to Idelchik data)
    double corr_peak = orifice::thickness_correction(0.3, 0.5, Re_d);
    EXPECT_GT(corr_peak, corr_small);

    // Long tube: friction dominates, k_t < 1.0
    double corr_long = orifice::thickness_correction(3.0, 0.5, Re_d);
    EXPECT_LT(corr_long, corr_peak);  // Falls at large t/d
    EXPECT_LT(corr_long, 1.0);  // Long-tube behavior
}
// -------------------------------------------------------------
// Hardening and Numerical Stability tests
// -------------------------------------------------------------

TEST_F(OrificeTest, HardeningHighBeta) {
    // Test that Cd remains finite for beta -> 1.0 (approaching pipe diameter)
    OrificeGeometry geom_extreme;
    geom_extreme.D = 0.100;
    geom_extreme.d = 0.0999999; // beta = 0.9999999

    double Cd_val = Cd_sharp_thin_plate(geom_extreme, state);
    // Should be capped at 1.5, not 500,000
    EXPECT_LE(Cd_val, 1.5);
    EXPECT_GT(Cd_val, 0.5);
}

TEST_F(OrificeTest, KFromCdHighBeta) {
    // Test that K_from_Cd and Cd_from_K are stable and self-consistent near beta = 1.0
    double beta_extreme = 0.99999;
    double Cd_val = 0.95;

    double K = orifice::K_from_Cd(Cd_val, beta_extreme);
    EXPECT_TRUE(std::isfinite(K));
    EXPECT_GE(K, 0.0);

    // Round-trip should be exact even with clamping because both use the same internal clamp
    double Cd_back = orifice::Cd_from_K(K, beta_extreme);
    EXPECT_NEAR(Cd_back, Cd_val, 1e-10);
}

TEST_F(OrificeTest, CdRoundedHighBetaStability) {
    // Test that Rounded-entry correlation doesn't divide by zero at d=D
    double beta_extreme = 0.99999;
    geom.r = 0.005; // Rounded

    double Cd_val = orifice::Cd_rounded(geom.r / (geom.D * beta_extreme), beta_extreme, state.Re_D);
    EXPECT_TRUE(std::isfinite(Cd_val));
    EXPECT_LE(Cd_val, 1.5);
}

// -------------------------------------------------------------
// McGreehan and Schotsch (1988)
//
// Every expectation below is anchored OUTSIDE the implementation: a value
// printed elsewhere in the paper, a different author's correlation, a curve
// drawn on the page, a physical identity, or another paper entirely. None is
// a golden value captured from this code. Provenance and the measured
// agreement for each are in
// validation/cooling/extractions/orifice_discharge_coefficient.md.
// -------------------------------------------------------------

namespace ms = orifice::mcgreehan_schotsch;

// The paper states a baseline sharp-edged Cd of 0.60 at Re = 3.2e4 in running
// text on p.213, independently of Eq. (8) on p.214. Eq. (8) must reproduce it.
TEST(McGreehanSchotsch, ReynoldsBaselineReproducesStatedReferencePoint) {
    EXPECT_NEAR(ms::reynolds_baseline(3.2e4), ms::cd_reference, 5e-4);
}

// Eq. (12) against Schoder and Dawson's Eq. (10), dCd/Cd = 3.10 (r/d), which
// the paper cites as "a good approximation for 0 < r/d < 0.1" but does not use.
// A different author's correlation, so this could genuinely have disagreed.
TEST(McGreehanSchotsch, CornerFactorAgreesWithSchoderDawson) {
    for (const double rd : {0.02, 0.05, 0.10}) {
        const double schoder = ms::cd_reference * (1.0 + 3.10 * rd);
        const double got     = 1.0 - ms::corner_factor(rd) * (1.0 - ms::cd_reference);
        EXPECT_NEAR(got, schoder, 0.015 * schoder) << "at r/d = " << rd;
    }
}

// p.214: an ASME nozzle Cd is the limiting value, reached at r/d = 0.82.
// Extrapolating the corner-radius fit that far must land on Eq. (9)'s own
// high-Re nozzle asymptote -- a different equation fitted to different data.
TEST(McGreehanSchotsch, CornerRadiusExtrapolatesToNozzleAsymptote) {
    const double fully_rounded =
        1.0 - ms::corner_factor(ms::r_over_d_nozzle_limit) * (1.0 - ms::re_c0);
    EXPECT_NEAR(fully_rounded, ms::nozzle_c0, 0.005 * ms::nozzle_c0);
}

// Identities the fitted constants must satisfy for the chain to be consistent.
// g(0) = 1 is a property of 1.3 and 0.435 conspiring, not something imposed.
TEST(McGreehanSchotsch, FactorsAreIdentitiesAtZero) {
    EXPECT_NEAR(ms::length_factor(0.0), 1.0, 1e-3);
    EXPECT_DOUBLE_EQ(ms::corner_factor(0.0), 1.0);
}

// Eq. (17) must reduce exactly to its input at zero crossflow.
TEST(McGreehanSchotsch, CrossflowIsIdentityAtZero) {
    const double base = ms::cd_with_corner_and_length(1.0e4, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(ms::cd(1.0e4, 0.0, 1.0, 0.0), base);
}

// Eq. (13)/(14) against the curve drawn in Fig. 3, read at ~0.01 in Cd.
// Fig. 3 is drawn for a basic Cd of 0.60, which Eq. (8) gives at Re = 3.2e4.
// L/d = 0.5 is excluded: it sits on the steepest part of the curve, where a
// 0.1 error in reading L/d moves Cd by 0.017 (see check E of the extraction).
TEST(McGreehanSchotsch, LengthEffectMatchesFigure3) {
    const struct { double L_over_d; double Cd_drawn; } points[] = {
        {1.0, 0.77}, {2.0, 0.81}, {3.0, 0.81}, {5.0, 0.78}, {8.0, 0.76}, {10.0, 0.74},
    };
    for (const auto& p : points) {
        const double got = ms::cd_with_corner_and_length(3.2e4, 0.0, p.L_over_d);
        EXPECT_NEAR(got, p.Cd_drawn, 0.02 * p.Cd_drawn) << "at L/d = " << p.L_over_d;
    }
}

// Eq. (17) against the Cd:r,L = 0.6 curve drawn in Figs. 4 and 7. A basic Cd
// of 0.60 is reached at Re = 3.2e4 with r/d = L/d = 0.
TEST(McGreehanSchotsch, CrossflowEffectMatchesFigure4) {
    const struct { double U1_over_Vi; double Cd_drawn; } points[] = {
        {1.0, 0.40}, {2.0, 0.24}, {4.0, 0.12},
    };
    for (const auto& p : points) {
        const double got = ms::cd(3.2e4, 0.0, 0.0, p.U1_over_Vi);
        EXPECT_NEAR(got, p.Cd_drawn, 0.05 * p.Cd_drawn) << "at U1/Vi = " << p.U1_over_Vi;
    }
}

// The rise before the fall is drawn on Fig. 4 and supported by Rohde's own
// result 2 (slanting an orifice into the flow increases Cd). It is deliberate,
// so it gets a test: anyone who "fixes" Eq. (17) into a monotone decay, or
// clamps the factor at 1, breaks this.
TEST(McGreehanSchotsch, CrossflowIsNonMonotonicNearTheOrigin) {
    const double base = ms::cd(3.2e4, 0.0, 0.0, 0.0);

    double peak = base;
    double peak_at = 0.0;
    for (int i = 1; i <= 400; ++i) {
        const double u  = 0.005 * i;
        const double cd = ms::cd(3.2e4, 0.0, 0.0, u);
        if (cd > peak) { peak = cd; peak_at = u; }
    }

    // Fig. 4's drawn peak is ~0.635 at U1/Vi ~ 0.1, against a 0.60 baseline.
    EXPECT_GT(peak, base) << "the rise before the fall has been lost";
    EXPECT_NEAR(peak, 0.635, 0.02 * 0.635);
    EXPECT_NEAR(peak_at, 0.1, 0.05);
    // ...and it really does fall away afterwards.
    EXPECT_LT(ms::cd(3.2e4, 0.0, 0.0, 4.0), 0.5 * base);
}

// The #375 anchor. A bare plenum-fed jet plate -- sharp holes, t/d = 1, no
// supply-side crossflow -- against Florschuetz, Truman and Metzger (1981),
// who recommend C_D = 0.79 and measure 0.73-0.85 (their Table 1). A different
// paper, rig and decade, and the value combaero currently hard-codes.
TEST(McGreehanSchotsch, BareJetPlateAgreesWithFlorschuetzDefault) {
    EXPECT_NEAR(ms::cd(1.0e4, 0.0, 1.0, 0.0), 0.79, 0.01 * 0.79);
}

TEST(McGreehanSchotsch, JetPlateStaysInsideFlorschuetzMeasuredBand) {
    for (const double Re : {5.0e3, 1.0e4, 3.0e4, 7.0e4}) {
        for (const double t_over_d : {1.0, 1.5, 2.0, 3.0}) {
            const double got = ms::cd(Re, 0.0, t_over_d, 0.0);
            EXPECT_GE(got, 0.73) << "Re = " << Re << ", t/d = " << t_over_d;
            EXPECT_LE(got, 0.85) << "Re = " << Re << ", t/d = " << t_over_d;
        }
    }
}

// Eqs. (15)/(16) only engage when a corner radius is present, so the chain has
// a branch at r/d = 0. g(0) = 1.0005 rather than exactly 1, so the branch is
// not perfectly continuous. Bound the step, so that if anyone widens it the
// test says so.
TEST(McGreehanSchotsch, CombinedCorrectionBranchIsEffectivelyContinuous) {
    const double at_zero = ms::cd_with_corner_and_length(3.2e4, 0.0, 2.0);
    const double just_above = ms::cd_with_corner_and_length(3.2e4, 1e-9, 2.0);
    EXPECT_NEAR(just_above, at_zero, 5e-4);
}

// Eqs. (13)+(15)+(16) collapse algebraically to a single product form,
//   Cd = 1 - g(L/d - r/d) * g(r/d) * (1 - Cd:r),
// derived by hand from the two sequential steps the implementation actually
// performs. An independent expression of the same thing: it catches a dropped
// Eq. (16) subtraction, a g() evaluated at the wrong argument, or the printed
// "L/D" being taken literally.
TEST(McGreehanSchotsch, CombinedCorrectionMatchesItsClosedProductForm) {
    const struct { double r_over_d; double L_over_d; } cases[] = {
        {0.1, 1.0}, {0.3, 2.0}, {0.5, 3.0}, {0.05, 0.5},
    };
    for (const auto& c : cases) {
        const double cd_r = ms::cd_with_corner(3.2e4, c.r_over_d);
        const double expected =
            1.0 - ms::length_factor(c.L_over_d - c.r_over_d) *
                      ms::length_factor(c.r_over_d) * (1.0 - cd_r);
        EXPECT_NEAR(ms::cd_with_corner_and_length(3.2e4, c.r_over_d, c.L_over_d),
                    expected, 1e-12)
            << "at r/d = " << c.r_over_d << ", L/d = " << c.L_over_d;
    }
}

// Rounding the inlet must raise Cd monotonically, and at the r/d = 0.82 knee
// the paper says an ASME nozzle is reached -- so a long, well-rounded orifice
// must sit close to 1.
TEST(McGreehanSchotsch, CornerRadiusRaisesCdTowardTheNozzleLimit) {
    double previous = ms::cd_with_corner_and_length(3.2e4, 0.0, 2.0);
    for (const double rd : {0.05, 0.1, 0.2, 0.3, 0.5, 0.82}) {
        const double got = ms::cd_with_corner_and_length(3.2e4, rd, 2.0);
        EXPECT_GT(got, previous) << "not monotone at r/d = " << rd;
        previous = got;
    }
    EXPECT_GT(previous, 0.98);
    EXPECT_LE(previous, 1.0);
}

TEST(McGreehanSchotsch, WrapperMatchesNamespacedChain) {
    EXPECT_DOUBLE_EQ(orifice::Cd_McGreehanSchotsch(2.0e4, 0.1, 1.5, 0.3),
                     ms::cd(2.0e4, 0.1, 1.5, 0.3));
}

TEST(McGreehanSchotsch, StaysFiniteOnSolverExcursions) {
    for (const double Re : {1.0e-3, 1.0, 1.0e3, 1.0e9}) {
        for (const double u : {0.0, 1.0e-12, 50.0}) {
            const double got = ms::cd(Re, 0.0, 1.0, u);
            EXPECT_TRUE(std::isfinite(got)) << "Re = " << Re << ", U1/Vi = " << u;
            EXPECT_GT(got, 0.0);
            EXPECT_LE(got, 1.5);
        }
    }
    // Negative inputs are clamped rather than producing NaN from pow().
    EXPECT_TRUE(std::isfinite(ms::cd(1.0e4, -1.0, -1.0, -1.0)));
}

// A constant Cd must be settable, not just selectable: make_correlation() can
// only hand back the default.
TEST(McGreehanSchotsch, ConstantCorrelationTakesAnExplicitValue) {
    OrificeGeometry g;
    g.d = 0.001;
    g.D = 0.010;
    OrificeState s;
    s.Re_D = 1.0e4;

    auto fixed = make_constant_correlation(0.79);
    ASSERT_NE(fixed, nullptr);
    EXPECT_DOUBLE_EQ(fixed->Cd(g, s), 0.79);
}

// Eq. (8) diverges below its stated floor (it reaches Cd = 1.0 at Re = 904),
// so the chain holds Re at re_min rather than extrapolating into nonsense.
// Constant below the floor is a deliberate, documented choice.
TEST(McGreehanSchotsch, ReynoldsBaselineIsHeldAtItsValidityFloor) {
    // Eq. (8) at its own stated floor: 0.5885 + 372/1e4.
    const double at_floor = ms::reynolds_baseline(ms::re_min);
    EXPECT_NEAR(at_floor, 0.6257, 1e-4);

    for (const double Re : {1.0e-3, 1.0, 9.0e2, 5.0e3}) {
        EXPECT_DOUBLE_EQ(ms::reynolds_baseline(Re), at_floor)
            << "not held at the floor for Re = " << Re;
    }
    // Above the floor it must still respond to Re.
    EXPECT_LT(ms::reynolds_baseline(1.0e5), at_floor);
    EXPECT_NEAR(ms::reynolds_baseline(3.2e4), ms::cd_reference, 5e-4);
}

// Eq. (17) applied to a caller-supplied baseline. This is how the source uses
// it in its own validation, so the chain must decompose exactly this way.
TEST(McGreehanSchotsch, ChainDecomposesIntoBaselineAndCrossflow) {
    for (const double u : {0.0, 0.1, 0.5, 1.0, 3.0}) {
        const double base = ms::cd_with_corner_and_length(3.2e4, 0.1, 2.0);
        EXPECT_DOUBLE_EQ(ms::cd(3.2e4, 0.1, 2.0, u),
                         ms::cd_with_crossflow(base, u))
            << "at U1/Vi = " << u;
    }
}

// Eq. (17) is otherwise validated against a DRAWN curve by
// CrossflowEffectMatchesFigure4, and against Rohde's digitised data by
// python/tests/test_orifice_validation.py. No figure-read test is duplicated
// here for the Fig. 6 baseline curves.

// -------------------------------------------------------------
// Expansion factor Y, Eqs. (4)-(7)
// -------------------------------------------------------------

TEST(McGreehanSchotschY, NozzleFormIsTheIsentropicExpansionFactor) {
    // Eq. (5) IS the isentropic nozzle expansion factor. Checked here against
    // the closed-form isentropic mass flux ratio, derived independently:
    //   Y_n = [mdot_isentropic / mdot_incompressible] at the same dP.
    // (The stronger check, against combaero's own nozzle_flow with real gas
    // properties, lives in python/tests/test_orifice_expansion.py.)
    const double g = 1.4;
    for (const double S : {0.99, 0.95, 0.90, 0.80, 0.70, 0.60}) {
        // Isentropic mass flux, normalised the same way Y is defined.
        const double num = std::sqrt(
            (2.0 * g / (g - 1.0)) *
            (std::pow(S, 2.0 / g) - std::pow(S, (g + 1.0) / g)));
        const double den = std::sqrt(2.0 * (1.0 - S));
        EXPECT_NEAR(ms::expansion_nozzle(S, g), num / den, 1e-10)
            << "at S = " << S;
    }
}

TEST(McGreehanSchotschY, BothFormsTendToOneAtZeroPressureDrop) {
    for (const double g : {1.3, 1.4, 1.67}) {
        EXPECT_NEAR(ms::expansion_orifice(1.0, g), 1.0, 1e-12);
        EXPECT_NEAR(ms::expansion_nozzle(1.0, g), 1.0, 1e-9);
        EXPECT_NEAR(ms::expansion_nozzle(1.0 - 1e-10, g), 1.0, 1e-6);
    }
}

TEST(McGreehanSchotschY, OrificeFormMatchesThePrintedCoefficient) {
    // Eq. (4) printed as 1 - 0.41 (P_t1 - P_s2)/(gamma P_t1). Checked against
    // that literal form rather than the (1-S)/gamma rewrite used internally.
    const double g = 1.4, P_t1 = 4.0e5, P_s2 = 3.0e5;
    const double printed = 1.0 - 0.41 * (P_t1 - P_s2) / (g * P_t1);
    EXPECT_NEAR(ms::expansion_orifice(P_s2 / P_t1, g), printed, 1e-12);
}

TEST(McGreehanSchotschY, BlendWeightReducesToThePaperAtZeroSmoothing) {
    // eps = 0 must recover Eq. (7) clamped hard, exactly.
    for (const double cd : {0.70, 0.80, 0.82, 0.86, 0.90, 0.94, 0.99}) {
        const double x = ms::x_blend_slope * (cd - ms::x_blend_cd0);
        const double hard = std::max(0.0, std::min(1.0, x));
        EXPECT_NEAR(ms::expansion_blend_weight(cd, 0.0), hard, 1e-12)
            << "at Cd = " << cd;
    }
    // X reaches 1 at Cd = 0.94, as the paper's constants imply.
    EXPECT_NEAR(ms::expansion_blend_weight(0.94, 0.0), 1.0, 1e-3);
}

// The reviewer's concern, made a test: a hard clamp puts an exactly-zero
// derivative outside [0.82, 0.94] and a discontinuous jump at both knees.
// The smoothed default must have neither.
TEST(McGreehanSchotschY, SmoothedBlendHasNoDeadZoneAndNoJump) {
    const double g = 1.4, S = 0.7;
    const double h = 1e-6;
    auto dY = [&](double cd, double eps) {
        return (ms::expansion_factor(cd + h, S, g, eps)
                - ms::expansion_factor(cd - h, S, g, eps)) / (2.0 * h);
    };

    // Hard clamp: dead on both sides. This documents what is being avoided.
    EXPECT_DOUBLE_EQ(dY(0.78, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(dY(0.99, 0.0), 0.0);

    // Smoothed: not merely non-zero but with ROOM in it. An earlier version
    // of this assertion tested != 0, which passes at eps = 0.005 where the
    // smallest slope is 1.4e-6 -- numerically a floor. 5e-5 is below the
    // 1.36e-4 the default achieves and above the 2.2e-5 that eps = 0.02 gets,
    // so this assertion is itself a lower wall on eps.
    for (const double cd : {0.70, 0.78, 0.818, 0.822, 0.90, 0.936, 0.944,
                            0.97, 0.99}) {
        EXPECT_GT(std::abs(dY(cd, ms::x_smooth_eps)), 5.0e-5)
            << "derivative is effectively floored at Cd = " << cd;
    }

    // ...and CONTINUOUS across both knees. Tested by how the change scales
    // with the sampling interval rather than against a threshold: a genuine
    // discontinuity holds the same jump however closely you sample, while a
    // continuous derivative's change shrinks with the interval. Measured, the
    // hard clamp sits at 0.73353 for every h; the smoothed one halves each
    // time h halves.
    for (const double knee : {0.82, 0.94}) {
        auto jump = [&](double half, double eps) {
            return std::abs(dY(knee + half, eps) - dY(knee - half, eps));
        };

        // Hard clamp: the jump does not shrink. This is what is being avoided.
        EXPECT_NEAR(jump(0.004, 0.0), jump(0.00025, 0.0), 1e-6);
        EXPECT_NEAR(jump(0.004, 0.0), 0.7335, 1e-3);

        // Smoothed: halving the interval halves the change, to within 5%.
        const double wide = jump(0.002, ms::x_smooth_eps);
        const double half = jump(0.001, ms::x_smooth_eps);
        EXPECT_NEAR(half / wide, 0.5, 0.05)
            << "derivative is not continuous at the knee Cd = " << knee;
    }
}

TEST(McGreehanSchotschY, SmoothingCostIsBounded) {
    // What the smoothing buys must be paid for in a quantified amount of Y.
    double worst = 0.0;
    for (int i = 0; i <= 300; ++i) {
        const double cd = 0.70 + 0.30 * i / 300.0;
        for (const double S : {0.6, 0.7, 0.8, 0.9, 0.99}) {
            worst = std::max(worst,
                             std::abs(ms::expansion_factor(cd, S, 1.4, ms::x_smooth_eps)
                                      - ms::expansion_factor(cd, S, 1.4, 0.0)));
        }
    }
    EXPECT_LT(worst, 0.006) << "smoothing deviation grew: " << worst;
}

// Below the critical pressure ratio the isentropic form turns over and
// predicts decreasing flow. Y must saturate instead.
// The evidence for x_smooth_eps, rather than an assurance that it is sensible.
//
// Two independent properties bound it from opposite directions, and the
// default has to sit between them. Measured walls: continuity is satisfied
// for eps >= 0.0304, fidelity for eps <= 0.0976. The default sits 29% across
// that window. If someone moves the constant, this says which wall they hit.
TEST(McGreehanSchotschY, DefaultSmoothingSitsInsideItsAdmissibleWindow) {
    const double g = 1.4, S = 0.7, h = 1e-6;
    auto dY = [&](double cd, double eps) {
        return (ms::expansion_factor(cd + h, S, g, eps)
                - ms::expansion_factor(cd - h, S, g, eps)) / (2.0 * h);
    };
    // Continuity, as the ratio of the derivative change at two sampling
    // intervals. 0.5 means continuous; 1.0 means a kink.
    auto continuity = [&](double eps) {
        auto jump = [&](double half) {
            return std::abs(dY(0.94 + half, eps) - dY(0.94 - half, eps));
        };
        const double wide = jump(0.002);
        return wide > 0.0 ? jump(0.001) / wide : 1.0;
    };
    // Departure from the paper's exact Eq. (7).
    auto cost = [&](double eps) {
        double worst = 0.0;
        for (int i = 0; i <= 300; ++i) {
            const double cd = 0.70 + 0.30 * i / 300.0;
            for (const double s : {0.6, 0.7, 0.8, 0.9, 0.99}) {
                worst = std::max(worst,
                                 std::abs(ms::expansion_factor(cd, s, g, eps)
                                          - ms::expansion_factor(cd, s, g, 0.0)));
            }
        }
        return worst;
    };

    // The default satisfies both.
    EXPECT_NEAR(continuity(ms::x_smooth_eps), 0.5, 0.05);
    EXPECT_LT(cost(ms::x_smooth_eps), 0.006);

    // Below the window, continuity goes -- the transition becomes narrower
    // than anything a Newton step can resolve.
    EXPECT_GT(std::abs(continuity(0.02) - 0.5), 0.05)
        << "the lower wall has moved; re-measure the window";

    // Above it, fidelity goes.
    EXPECT_GT(cost(0.15), 0.006)
        << "the upper wall has moved; re-measure the window";

    // And the window really is narrow, so the default is not free to drift.
    EXPECT_LT(cost(ms::x_smooth_eps) / 0.006, 0.9);
}

TEST(McGreehanSchotschY, SaturatesAtChoking) {
    const double g = 1.4;
    const double s_crit = ms::critical_pressure_ratio(g);
    EXPECT_NEAR(s_crit, 0.5283, 1e-3);

    auto flux = [&](double S) {
        return ms::expansion_factor(0.99, S, g) * std::sqrt(std::max(1.0 - S, 0.0));
    };
    // The proxy mass flux must not fall away below the choke point.
    const double at_crit = flux(s_crit);
    for (const double S : {0.45, 0.35, 0.20, 0.05}) {
        EXPECT_GT(flux(S), 0.97 * at_crit) << "flow collapses below S* at S = " << S;
    }
}

// Both bounds on the pressure ratio must be soft, not just the choke point.
// A scan for zero-derivative regions and derivative discontinuities caught an
// earlier revision that smoothed S* but left a hard min(S, 1): it put a kink
// of ~0.455 in dY/dS at S = 1 and a floor above it. S > 1 is reverse flow,
// which a solver iterate can reach.
TEST(McGreehanSchotschY, PressureRatioIsSmoothAtBothBounds) {
    const double g = 1.4, cd = 0.90, h = 1e-7;
    auto dS = [&](double S) {
        return (ms::expansion_factor(cd, S + h, g)
                - ms::expansion_factor(cd, S - h, g)) / (2.0 * h);
    };
    auto continuous_at = [&](double S0) {
        auto jump = [&](double half) {
            return std::abs(dS(S0 + half) - dS(S0 - half));
        };
        const double wide = jump(0.004);
        if (wide < 1e-9) return true;   // already flat on both sides
        return jump(0.002) / wide < 0.75;  // shrinks with the interval
    };

    EXPECT_TRUE(continuous_at(1.0)) << "kink at the no-flow bound S = 1";
    EXPECT_TRUE(continuous_at(ms::critical_pressure_ratio(g)))
        << "kink at the choke point";

    // No floor anywhere a solver can wander, including into reverse flow.
    for (const double S : {0.2, 0.4, ms::critical_pressure_ratio(g), 0.7,
                           0.95, 0.999, 1.0, 1.05, 1.2}) {
        EXPECT_NE(dS(S), 0.0) << "derivative floored at S = " << S;
    }
}

TEST(McGreehanSchotschY, BlendSelectsOrificeWhenSharpAndNozzleWhenRounded) {
    const double g = 1.4, S = 0.7;
    // A sharp hole (Cd well under 0.82) is pure orifice.
    EXPECT_NEAR(ms::expansion_factor(0.70, S, g),
                ms::expansion_orifice(S, g), 5e-3);
    // One rounded enough to act like a nozzle is pure nozzle.
    EXPECT_NEAR(ms::expansion_factor(0.995, S, g),
                ms::expansion_nozzle(S, g), 5e-3);
    // They genuinely differ, so the blend is doing work.
    EXPECT_GT(std::abs(ms::expansion_orifice(S, g) - ms::expansion_nozzle(S, g)),
              0.08);
}

TEST(McGreehanSchotschY, StaysFiniteAndBounded) {
    for (const double g : {1.1, 1.4, 1.67}) {
        for (const double S : {1.5, 1.0, 0.5, 0.1, 0.0, -0.1}) {
            for (const double cd : {0.3, 0.8, 1.2}) {
                const double y = ms::expansion_factor(cd, S, g);
                EXPECT_TRUE(std::isfinite(y)) << "g=" << g << " S=" << S;
                EXPECT_GT(y, 0.0);
                EXPECT_LE(y, 1.01);
            }
        }
    }
    EXPECT_DOUBLE_EQ(ms::expansion_factor(0.8, 0.7, 1.0), 1.0);  // degenerate gamma
}
