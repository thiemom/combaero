#include <gtest/gtest.h>
#include "../include/incompressible.h"
#include <cmath>
#include <stdexcept>

// -------------------------------------------------------------
// pressure_loss / velocity_from_pressure_loss
// -------------------------------------------------------------

TEST(PressureLoss, RoundTrip) {
    const double v = 10.0;
    const double rho = 1.2;
    const double zeta = 2.5;
    double dP = pressure_loss(v, rho, zeta);
    EXPECT_DOUBLE_EQ(dP, zeta * 0.5 * rho * v * v);
    EXPECT_NEAR(velocity_from_pressure_loss(dP, rho, zeta), v, 1e-12);
}

TEST(PressureLoss, ZetaOne_EqualsOneDynamicPressure) {
    const double v = 5.0;
    const double rho = 1.2;
    EXPECT_DOUBLE_EQ(pressure_loss(v, rho, 1.0), dynamic_pressure(v, rho));
}

TEST(PressureLoss, ZetaZero_NoPressureDrop) {
    EXPECT_DOUBLE_EQ(pressure_loss(20.0, 1.2, 0.0), 0.0);
}

TEST(PressureLoss, InvalidRho) {
    EXPECT_THROW(pressure_loss(10.0, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(pressure_loss(10.0, -1.0, 1.0), std::invalid_argument);
}

TEST(PressureLoss, InvalidZeta) {
    EXPECT_THROW(pressure_loss(10.0, 1.2, -0.1), std::invalid_argument);
}

TEST(VelocityFromPressureLoss, InvalidArgs) {
    EXPECT_THROW(velocity_from_pressure_loss(-1.0, 1.2, 1.0), std::invalid_argument);
    EXPECT_THROW(velocity_from_pressure_loss(100.0, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(velocity_from_pressure_loss(100.0, 1.2, 0.0), std::invalid_argument);
}

// -------------------------------------------------------------
// zeta_from_Cd / Cd_from_zeta
// -------------------------------------------------------------

TEST(ZetaCd, RoundTrip) {
    for (double Cd : {0.3, 0.6, 0.61, 0.8, 1.0}) {
        EXPECT_NEAR(Cd_from_zeta(zeta_from_Cd(Cd)), Cd, 1e-12);
    }
}

TEST(ZetaCd, KnownValues) {
    EXPECT_DOUBLE_EQ(zeta_from_Cd(1.0), 1.0);
    EXPECT_NEAR(zeta_from_Cd(0.5), 4.0, 1e-12);
    EXPECT_NEAR(Cd_from_zeta(4.0), 0.5, 1e-12);
}

TEST(ZetaCd, SharpEdgeOrifice) {
    // Typical sharp-edge Cd ≈ 0.61 → ζ ≈ 2.69
    const double Cd = 0.61;
    const double zeta = zeta_from_Cd(Cd);
    EXPECT_NEAR(zeta, 1.0 / (Cd * Cd), 1e-12);
    EXPECT_NEAR(Cd_from_zeta(zeta), Cd, 1e-12);
}

TEST(ZetaCd, InvalidArgs) {
    EXPECT_THROW(zeta_from_Cd(0.0), std::invalid_argument);
    EXPECT_THROW(zeta_from_Cd(-0.5), std::invalid_argument);
    EXPECT_THROW(Cd_from_zeta(0.0), std::invalid_argument);
    EXPECT_THROW(Cd_from_zeta(-1.0), std::invalid_argument);
}

// -------------------------------------------------------------
// Consistency: zeta-based dP matches Cd-based orifice_dP
// -------------------------------------------------------------

TEST(ZetaCd, ConsistentWithOrificeDP) {
    // For a restriction where reference area = throat area:
    //   orifice_dP(mdot, A, Cd, rho) == pressure_loss(v_throat, rho, zeta_from_Cd(Cd))
    const double A = 1e-4;   // 1 cm²
    const double Cd = 0.61;
    const double rho = 1.2;
    const double mdot = 0.05;

    double dP_orifice = orifice_dP(mdot, A, Cd, rho);
    double v_throat = mdot / (rho * A);
    double dP_zeta = pressure_loss(v_throat, rho, zeta_from_Cd(Cd));
    EXPECT_NEAR(dP_zeta, dP_orifice, 1e-9);
}
