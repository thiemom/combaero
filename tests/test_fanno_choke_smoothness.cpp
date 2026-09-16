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
constexpr double kRough = 1e-5;
constexpr double kT = 1700.0;
constexpr double kPt = 130000.0;

double dP_at(double m_dot, const std::vector<double>& Y) {
    const auto res = solver::channel_compressible_residuals_and_jacobian(
        m_dot, kPt, kT, Y, 101325.0, kL, kD, kRough, "haaland", 1.0);
    return res.dP_calc;
}

}  // namespace

// The Fanno march truncates at a step boundary when it chokes. L_choke is
// interpolated within the breaking step; if the outlet state is not, the
// reported drop snaps to the march grid and the gradient develops a sawtooth
// -- measured at 22% peak-to-peak before the fix.
TEST(FannoChokeSmoothnessTest, ChokedBranchGradientHasNoMarchGridSawtooth) {
    const std::vector<double> Y = mole_to_mass(dry_air());

    // Inside the choked band for this duct.
    const double m_lo = 1.560;
    const double step = 5e-4;
    const int n = 30;

    std::vector<double> slopes;
    double prev = dP_at(m_lo, Y);
    for (int i = 1; i <= n; ++i) {
        const double m = m_lo + i * step;
        const double cur = dP_at(m, Y);
        slopes.push_back((cur - prev) / step);
        prev = cur;
    }

    // Consecutive slopes must differ by only the smooth physical trend. The
    // grid staircase alternated by ~22%; require well under that.
    double worst = 0.0;
    for (std::size_t i = 1; i < slopes.size(); ++i) {
        const double rel = std::abs(slopes[i] - slopes[i - 1]) /
                           std::max(std::abs(slopes[i - 1]), 1.0);
        worst = std::max(worst, rel);
    }
    EXPECT_LT(worst, 0.05) << "gradient sawtooth in the choked band: "
                           << worst * 100.0 << "% between consecutive samples";
}

// The drop must not fall as m_dot rises. It used to: the truncated drop covers
// the marched part only, so it collapses as choking moves upstream, and with a
// weak barrier the composite went non-monotone -- which is precisely the
// spurious-root condition the barrier exists to prevent.
TEST(FannoChokeSmoothnessTest, DropIsMonotoneThroughChokeOnset) {
    const std::vector<double> Y = mole_to_mass(dry_air());

    double prev = dP_at(1.500, Y);
    for (double m = 1.502; m <= 1.660; m += 0.002) {
        const double cur = dP_at(m, Y);
        EXPECT_GE(cur, prev) << "drop decreased with rising m_dot at m = " << m;
        prev = cur;
    }
}

// When the march chokes, the reported outlet has to sit on the same station as
// L_choke. If it lags at the overshooting step boundary, refining the march
// moves the reported outlet -- the signature of grid snapping.
TEST(FannoChokeSmoothnessTest, ChokedOutletIsGridIndependent) {
    const std::vector<double> X = dry_air();
    const double T = 1700.0;
    const double M_in = 0.91;  // chokes before the duct end
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
