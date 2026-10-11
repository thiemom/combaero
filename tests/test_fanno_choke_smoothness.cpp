#include <gtest/gtest.h>

#include "../include/compressible.h"
#include "../include/humidair.h"
#include "../include/math_constants.h"
#include "../include/solver_interface.h"
#include "../include/thermo.h"

#include <cmath>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// Smoothness of the compressible channel drop through choke onset (issue #356)
//
// The network solver differentiates this drop by central finite differences
// with a 1e-6 relative step. Anything that makes the drop a staircase, or that
// kinks its gradient, is handed straight to Newton as a wrong Jacobian and
// shows up as "the iteration is not making good progress". These tests pin the
// two defects that did exactly that.
// ---------------------------------------------------------------------------

namespace {

constexpr double kL = 1.0;
constexpr double kD = 0.10;

}  // namespace

// The network's compressible channel no longer differences this march (it
// uses the Mach-quadrature Fanno flow, fanno_mach.h, #481); the march remains
// the profile API and an independent reference, so its choke station must
// stay grid independent.
TEST(FannoChokeSmoothnessTest, ChokedOutletIsGridIndependent) {
    const std::vector<double> X = dry_air();
    const double T = 1700.0;
    // Chokes well before the duct end: L* = 0.20 m against L = 1 m with the
    // corrected gradient (#362). Earlier values here (0.91, then 0.70) were
    // picked against a march that never choked -- 0.70 now has L* = 1.15 m and
    // traverses the duct, which is the right answer, not a regression.
    const double M_in = 0.85;
    const double a = speed_of_sound(T, X);
    const double u = M_in * a;
    const double A = 0.25 * M_PI * kD * kD;
    const double rho = 1.03 / (u * A);
    const double P = rho * specific_gas_constant(X) * T;

    const auto coarse = fanno_channel(T, P, u, kL, kD, 0.02, X, 100);
    const auto fine = fanno_channel(T, P, u, kL, kD, 0.02, X, 6400);
    ASSERT_TRUE(coarse.choked);
    ASSERT_TRUE(fine.choked);

    // Both report the choke station, so refining must not move the outlet by
    // more than the march's own truncation error.
    EXPECT_NEAR(coarse.outlet.P, fine.outlet.P, 20.0);
    EXPECT_NEAR(coarse.L_choke, fine.L_choke, 1e-3);
}
