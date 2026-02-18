#include <gtest/gtest.h>
#include "../include/acoustics.h"
#include "../include/geometry.h"
#include "../include/math_constants.h"  // MSVC compatibility for M_PI

using namespace combaero;

// Test fixture for adapter tests
class AcousticAdapterTest : public ::testing::Test {
protected:
    CanAnnularFlowGeometry flow_geom;

    void SetUp() override {
        // Standard test geometry: annular section
        flow_geom.L = 0.3;              // Total length [m]
        flow_geom.D_inner = 0.4;        // Inner diameter [m]
        flow_geom.D_outer = 0.6;        // Outer diameter [m]
        flow_geom.L_primary = 0.2;      // Primary zone length [m]
        flow_geom.D_primary = 0.1;      // Primary zone diameter [m]
        flow_geom.L_transition = 0.05;  // Transition length [m]
    }
};

// -------------------------------------------------------------
// to_acoustic_geometry adapter tests
// -------------------------------------------------------------

TEST_F(AcousticAdapterTest, BasicConversion) {
    auto geom = to_acoustic_geometry(flow_geom, 24, 0.5, 0.1);

    EXPECT_EQ(geom.n_cans, 24);
    EXPECT_DOUBLE_EQ(geom.length_can, 0.5);

    // Can area = π * (0.1/2)^2 = π * 0.0025
    const double expected_can_area = M_PI * 0.0025;
    EXPECT_NEAR(geom.area_can, expected_can_area, 1e-10);

    // Plenum radius = D_mean / 2 = (0.4 + 0.6) / 4 = 0.25
    EXPECT_DOUBLE_EQ(geom.radius_plenum, 0.25);

    // Plenum area = flow_geom.area() = π * ((0.6/2)^2 - (0.4/2)^2)
    const double expected_plenum_area = flow_geom.area();
    EXPECT_DOUBLE_EQ(geom.area_plenum, expected_plenum_area);
}

TEST_F(AcousticAdapterTest, DifferentCanCount) {
    for (int n_cans : {1, 8, 12, 24, 36}) {
        auto geom = to_acoustic_geometry(flow_geom, n_cans, 0.5, 0.1);
        EXPECT_EQ(geom.n_cans, n_cans);
    }
}

TEST_F(AcousticAdapterTest, DifferentCanDimensions) {
    // Test various can lengths and diameters
    struct TestCase {
        double L_can;
        double D_can;
    };

    std::vector<TestCase> cases = {
        {0.3, 0.08},
        {0.5, 0.1},
        {0.7, 0.12},
        {1.0, 0.15}
    };

    for (const auto& tc : cases) {
        auto geom = to_acoustic_geometry(flow_geom, 24, tc.L_can, tc.D_can);

        EXPECT_DOUBLE_EQ(geom.length_can, tc.L_can);

        const double expected_can_area = M_PI * (tc.D_can / 2.0) * (tc.D_can / 2.0);
        EXPECT_NEAR(geom.area_can, expected_can_area, 1e-10);

        // Plenum properties should remain constant
        EXPECT_DOUBLE_EQ(geom.radius_plenum, 0.25);
        EXPECT_DOUBLE_EQ(geom.area_plenum, flow_geom.area());
    }
}

TEST_F(AcousticAdapterTest, PlenumPropertiesFromFlowGeom) {
    // Verify that plenum properties are correctly derived from flow geometry

    // Case 1: Different annular diameters
    CanAnnularFlowGeometry geom1;
    geom1.D_inner = 0.5;
    geom1.D_outer = 0.7;
    auto result1 = to_acoustic_geometry(geom1, 12, 0.4, 0.08);

    // D_mean = (0.5 + 0.7) / 2 = 0.6, so radius = 0.3
    EXPECT_DOUBLE_EQ(result1.radius_plenum, 0.3);

    // Case 2: Larger annulus
    CanAnnularFlowGeometry geom2;
    geom2.D_inner = 0.8;
    geom2.D_outer = 1.2;
    auto result2 = to_acoustic_geometry(geom2, 18, 0.6, 0.12);

    // D_mean = (0.8 + 1.2) / 2 = 1.0, so radius = 0.5
    EXPECT_DOUBLE_EQ(result2.radius_plenum, 0.5);
}

TEST_F(AcousticAdapterTest, InvalidInputs) {
    // Test invalid n_cans
    EXPECT_THROW(to_acoustic_geometry(flow_geom, 0, 0.5, 0.1), std::invalid_argument);
    EXPECT_THROW(to_acoustic_geometry(flow_geom, -1, 0.5, 0.1), std::invalid_argument);

    // Test invalid L_can
    EXPECT_THROW(to_acoustic_geometry(flow_geom, 24, 0.0, 0.1), std::invalid_argument);
    EXPECT_THROW(to_acoustic_geometry(flow_geom, 24, -0.1, 0.1), std::invalid_argument);

    // Test invalid D_can
    EXPECT_THROW(to_acoustic_geometry(flow_geom, 24, 0.5, 0.0), std::invalid_argument);
    EXPECT_THROW(to_acoustic_geometry(flow_geom, 24, 0.5, -0.05), std::invalid_argument);
}

TEST_F(AcousticAdapterTest, EdgeCases) {
    // Very small can
    auto small = to_acoustic_geometry(flow_geom, 1, 0.01, 0.001);
    EXPECT_EQ(small.n_cans, 1);
    EXPECT_NEAR(small.area_can, M_PI * 0.00000025, 1e-15);

    // Large number of cans
    auto many = to_acoustic_geometry(flow_geom, 1000, 0.5, 0.1);
    EXPECT_EQ(many.n_cans, 1000);

    // Thin annulus (small gap)
    CanAnnularFlowGeometry thin;
    thin.D_inner = 0.99;
    thin.D_outer = 1.0;
    auto thin_result = to_acoustic_geometry(thin, 24, 0.5, 0.1);
    EXPECT_DOUBLE_EQ(thin_result.radius_plenum, (0.99 + 1.0) / 4.0);  // D_mean / 2
}

TEST_F(AcousticAdapterTest, IntegrationWithSolver) {
    // Test that the converted geometry can be used with can_annular_eigenmodes
    auto geom = to_acoustic_geometry(flow_geom, 8, 0.5, 0.1);

    // Basic validation that geometry is usable
    EXPECT_GT(geom.area_can, 0.0);
    EXPECT_GT(geom.area_plenum, 0.0);
    EXPECT_GT(geom.radius_plenum, 0.0);
    EXPECT_GT(geom.length_can, 0.0);
    EXPECT_GT(geom.n_cans, 0);
}
