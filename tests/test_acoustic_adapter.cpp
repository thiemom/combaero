#include <gtest/gtest.h>
#include "../include/acoustics.h"
#include "../include/can_annular_solvers.h"
#include "../include/geometry.h"
#include "../include/math_constants.h"  // MSVC compatibility for M_PI
#include <cmath>
#include <vector>

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

// -------------------------------------------------------------
// Both can-annular solver methods against the straight-duct identity
// -------------------------------------------------------------
// A_can = 2 A_plenum, L_can = s/2 (s = 2 pi R / N) and one gas: for m = 0 the
// can and the two closed half-sector stubs form one closed-closed duct of
// length 2 L seen from its middle: below 1200 Hz its m = 0 modes are exactly
// c / (4 L) = 500 Hz (pressure node at the junction) and c / (2 L) = 1000 Hz.
// With the plenum term of the wrong sign the admittance sum vanishes
// identically and every scan minimum is "a mode".
namespace {
CanAnnularGeometry straight_duct_geometry() {
    CanAnnularGeometry g;
    g.n_cans = 12;
    g.length_can = 0.25;
    g.area_can = 0.02;
    g.area_plenum = 0.01;
    g.radius_plenum = 2.0 * g.length_can * g.n_cans / (2.0 * M_PI);
    return g;
}
}  // namespace

TEST(CanAnnularSolverMethods, StraightDuctIdentityBothMethods) {
    const auto g = straight_duct_geometry();
    const double c = 500.0;
    for (auto method : {CanAnnularSolverMethod::MagnitudeMinimization,
                        CanAnnularSolverMethod::ArgumentPrinciple}) {
        const auto modes = can_annular_eigenmodes_with_method(
            g, c, c, 1.0, 1.0, 1200.0, BoundaryCondition::Closed, method);
        std::vector<double> m0;
        for (const auto& mode : modes) {
            if (mode.m_azimuthal == 0) m0.push_back(mode.frequency);
        }
        ASSERT_EQ(m0.size(), 2u) << "method " << static_cast<int>(method);
        EXPECT_NEAR(m0[0], c / (4.0 * g.length_can), 1e-6);
        // 1000 Hz is a multiple root of the pole-free product (sin(kL) and
        // sin(2kL) vanish together), so Muller converges only linearly there.
        EXPECT_NEAR(m0[1], c / (2.0 * g.length_can), 1e-3);
    }
}

// The Argument Principle method returns the same mode set as magnitude
// minimization, which python/tests/test_can_annular_acoustics.py checks root
// by root against an independently derived dispersion relation. With D as a
// bare admittance sum the winding number was zeros minus poles and the method
// dropped 4 of the 21 modes of the first case below; it also skipped roots
// next to a box edge and accepted non-converged Muller iterates.
TEST(CanAnnularSolverMethods, ArgumentPrincipleMatchesMagnitudeMinimization) {
    struct Case {
        int n;
        double L, Ac, R, Ap, c_c, c_p, r_c, r_p, f_max;
        BoundaryCondition bc;
    };
    const Case cases[] = {
        {12, 0.5, 0.01, 0.5, 0.02, 500.0, 600.0, 1.0, 1.2, 1200.0, BoundaryCondition::Closed},
        {24, 0.5, 0.01, 0.5, 0.02, 500.0, 600.0, 1.0, 1.2, 2000.0, BoundaryCondition::Closed},
        {16, 0.8, 0.03, 1.2, 0.05, 800.0, 450.0, 0.4, 1.8, 1500.0, BoundaryCondition::Closed},
        {2, 0.4, 0.01, 0.3, 0.01, 340.0, 340.0, 1.2, 1.2, 2000.0, BoundaryCondition::Closed},
        {12, 0.5, 0.01, 0.5, 0.02, 500.0, 600.0, 1.0, 1.2, 1500.0, BoundaryCondition::Open},
    };
    for (const auto& k : cases) {
        CanAnnularGeometry g;
        g.n_cans = k.n;
        g.length_can = k.L;
        g.area_can = k.Ac;
        g.radius_plenum = k.R;
        g.area_plenum = k.Ap;
        const auto ap = can_annular_eigenmodes_with_method(
            g, k.c_c, k.c_p, k.r_c, k.r_p, k.f_max, k.bc,
            CanAnnularSolverMethod::ArgumentPrinciple);
        const auto mm = can_annular_eigenmodes_with_method(
            g, k.c_c, k.c_p, k.r_c, k.r_p, k.f_max, k.bc,
            CanAnnularSolverMethod::MagnitudeMinimization);
        EXPECT_EQ(ap.size(), mm.size()) << "N=" << k.n;
        auto contains = [](const std::vector<BlochMode>& set, const BlochMode& x) {
            for (const auto& y : set) {
                if (x.m_azimuthal == y.m_azimuthal && std::abs(x.frequency - y.frequency) < 1e-3) {
                    return true;
                }
            }
            return false;
        };
        for (const auto& a : ap) {
            EXPECT_TRUE(contains(mm, a)) << "N=" << k.n << " extra m=" << a.m_azimuthal
                                         << " f=" << a.frequency;
        }
        for (const auto& b : mm) {
            EXPECT_TRUE(contains(ap, b)) << "N=" << k.n << " missed m=" << b.m_azimuthal
                                         << " f=" << b.frequency;
        }
    }
}
