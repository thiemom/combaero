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

    // Expected: Cd * A * sqrt(2 * rho * dP)
    double expected = Cd_val * geom.area() * std::sqrt(2.0 * state.rho * state.dP);
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
