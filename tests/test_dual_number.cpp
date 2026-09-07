// Forward-mode dual numbers: every operation's derivative checked against a
// finite difference of its own value.
//
// The repo's solver rule requires an analytic derivative, and requires it to
// be cross-checked independently -- two code paths must never lean on the same
// derivation. Here the independent check is a central difference of the SAME
// routine's primal, which is legitimate precisely because the primal and the
// partials are computed by different code inside each operator.
//
// This header was extracted from src/ejector.cpp for the junction (f, J) port
// (issue #271) and gained dexp, dsin, dcos, datan2, dabs and the two angle
// wraps, which the Mynard closure needs. The pre-existing operators are
// covered too: the extraction must not have changed them, and the ejector's
// own Jacobian tests are the other half of that gate.

#include <cmath>
#include <functional>

#include <gtest/gtest.h>

#include "dual_number.h"
#include "math_constants.h"

using combaero::solver::DualN;

namespace {

using D1 = DualN<1>;
using D2 = DualN<2>;

constexpr double kStep = 1e-6;
constexpr double kTol = 1e-6;

// Central difference of a scalar function, against which the dual partial is
// compared. Relative where the derivative is large, absolute where it is small.
void ExpectDerivative(const std::function<double(double)>& f,
                      const std::function<D1(D1)>& f_dual, double x,
                      double tol = kTol) {
  D1 out = f_dual(D1::seed(x, 0));
  double fd = (f(x + kStep) - f(x - kStep)) / (2.0 * kStep);
  EXPECT_NEAR(out.v, f(x), 1e-12) << "primal disagrees at x=" << x;
  EXPECT_NEAR(out.d[0], fd, tol * std::max(1.0, std::abs(fd)))
      << "derivative disagrees at x=" << x;
}

} // namespace

// -----------------------------------------------------------------------------
// Seeding and constants
// -----------------------------------------------------------------------------

TEST(DualNumber, SeedCarriesUnitPartialAndConstantCarriesNone) {
  D2 x = D2::seed(3.0, 0);
  D2 c = D2::constant(3.0);

  EXPECT_DOUBLE_EQ(x.v, 3.0);
  EXPECT_DOUBLE_EQ(x.d[0], 1.0);
  EXPECT_DOUBLE_EQ(x.d[1], 0.0);
  EXPECT_DOUBLE_EQ(c.d[0], 0.0);
  EXPECT_DOUBLE_EQ(c.d[1], 0.0);
}

TEST(DualNumber, PartialsStayIndependentAcrossSeeds) {
  // The whole point of N seeds: a product's partial in each direction is that
  // direction's alone.
  D2 x = D2::seed(3.0, 0);
  D2 y = D2::seed(5.0, 1);
  D2 p = x * y;

  EXPECT_DOUBLE_EQ(p.v, 15.0);
  EXPECT_DOUBLE_EQ(p.d[0], 5.0);
  EXPECT_DOUBLE_EQ(p.d[1], 3.0);
}

// -----------------------------------------------------------------------------
// Arithmetic
// -----------------------------------------------------------------------------

TEST(DualNumber, ArithmeticMatchesFiniteDifference) {
  for (double x : {-2.7, -0.4, 0.3, 1.0, 4.2}) {
    ExpectDerivative([](double v) { return v + 2.5; },
                     [](D1 v) { return v + 2.5; }, x);
    ExpectDerivative([](double v) { return 2.5 - v; },
                     [](D1 v) { return 2.5 - v; }, x);
    ExpectDerivative([](double v) { return -v; }, [](D1 v) { return -v; }, x);
    ExpectDerivative([](double v) { return v * v; },
                     [](D1 v) { return v * v; }, x);
    ExpectDerivative([](double v) { return v * 3.0; },
                     [](D1 v) { return v * 3.0; }, x);
    ExpectDerivative([](double v) { return (v * v + 1.0) / (v * v + 2.0); },
                     [](D1 v) { return (v * v + 1.0) / (v * v + 2.0); }, x);
    ExpectDerivative([](double v) { return 6.0 / (v * v + 1.0); },
                     [](D1 v) { return 6.0 / (v * v + 1.0); }, x);
  }
}

// -----------------------------------------------------------------------------
// Elementary functions
// -----------------------------------------------------------------------------

TEST(DualNumber, SqrtAndPowMatchFiniteDifference) {
  for (double x : {0.05, 0.5, 1.0, 3.3, 120.0}) {
    ExpectDerivative([](double v) { return std::sqrt(v); },
                     [](D1 v) { return dsqrt(v); }, x);
    ExpectDerivative([](double v) { return std::pow(v, 1.4); },
                     [](D1 v) { return dpow(v, 1.4); }, x);
    ExpectDerivative([](double v) { return std::pow(v, -0.7); },
                     [](D1 v) { return dpow(v, -0.7); }, x);
  }
}

TEST(DualNumber, ExpAndLogMatchFiniteDifference) {
  for (double x : {-3.0, -0.2, 0.4, 2.0}) {
    ExpectDerivative([](double v) { return std::exp(v); },
                     [](D1 v) { return dexp(v); }, x);
  }
  for (double x : {0.05, 0.9, 4.0, 900.0}) {
    ExpectDerivative([](double v) { return std::log(v); },
                     [](D1 v) { return dlog(v); }, x);
  }
}

TEST(DualNumber, ExpIsCorrectAtTheDampingRegulariser) {
  // The closure's one use of exp: 1 - exp(-flow_ratio / 0.02), whose argument
  // is large and negative for any ordinary split. Checked here because a naive
  // FD on the composite loses all its digits there, while the dual does not.
  constexpr double kTau = 0.02;
  for (double q : {0.01, 0.05, 0.2, 0.5, 0.9}) {
    D1 ratio = D1::seed(q, 0);
    D1 damping = 1.0 - dexp(0.0 - ratio / kTau);
    EXPECT_NEAR(damping.v, 1.0 - std::exp(-q / kTau), 1e-14);
    EXPECT_NEAR(damping.d[0], std::exp(-q / kTau) / kTau,
                1e-9 * std::max(1.0, std::exp(-q / kTau) / kTau));
  }
}

TEST(DualNumber, SinAndCosMatchFiniteDifference) {
  for (double x : {-3.0, -0.5, 0.0, 0.7, 2.4, 5.9}) {
    ExpectDerivative([](double v) { return std::sin(v); },
                     [](D1 v) { return dsin(v); }, x);
    ExpectDerivative([](double v) { return std::cos(v); },
                     [](D1 v) { return dcos(v); }, x);
  }
}

TEST(DualNumber, AbsTakesTheOneSidedDerivativeAtZero) {
  // Deliberate: a zero column there would drop the unknown from the Newton
  // step entirely, which is worse than a one-sided slope.
  D1 at_zero = dabs(D1::seed(0.0, 0));
  EXPECT_DOUBLE_EQ(at_zero.v, 0.0);
  EXPECT_DOUBLE_EQ(at_zero.d[0], 1.0);

  EXPECT_DOUBLE_EQ(dabs(D1::seed(-2.0, 0)).d[0], -1.0);
  EXPECT_DOUBLE_EQ(dabs(D1::seed(2.0, 0)).d[0], 1.0);
}

// -----------------------------------------------------------------------------
// atan2 -- the one with two dual arguments
// -----------------------------------------------------------------------------

TEST(DualNumber, Atan2MatchesFiniteDifferenceInEveryQuadrant) {
  const double pts[][2] = {{1.0, 2.0},  {2.0, -1.0}, {-1.5, 0.7},
                           {-2.0, -3.0}, {0.0, 4.0},  {4.0, 0.0}};
  for (const auto& p : pts) {
    double y0 = p[0];
    double x0 = p[1];
    D2 out = datan2(D2::seed(y0, 0), D2::seed(x0, 1));

    EXPECT_NEAR(out.v, std::atan2(y0, x0), 1e-14);
    double fd_y = (std::atan2(y0 + kStep, x0) - std::atan2(y0 - kStep, x0)) / (2.0 * kStep);
    double fd_x = (std::atan2(y0, x0 + kStep) - std::atan2(y0, x0 - kStep)) / (2.0 * kStep);
    EXPECT_NEAR(out.d[0], fd_y, kTol * std::max(1.0, std::abs(fd_y)))
        << "d/dy at (" << y0 << "," << x0 << ")";
    EXPECT_NEAR(out.d[1], fd_x, kTol * std::max(1.0, std::abs(fd_x)))
        << "d/dx at (" << y0 << "," << x0 << ")";
  }
}

TEST(DualNumber, Atan2DerivativeSurvivesTheBranchCut) {
  // atan2 jumps by 2*pi across the negative-x axis, but the jump is a
  // CONSTANT, so the derivative is continuous there. Approaching from both
  // sides must give the same partials -- a formula written as an arctangent
  // of y/x would not.
  const double x0 = -3.0;
  D2 above = datan2(D2::seed(1e-9, 0), D2::seed(x0, 1));
  D2 below = datan2(D2::seed(-1e-9, 0), D2::seed(x0, 1));

  EXPECT_NEAR(above.d[0], below.d[0], 1e-9);
  EXPECT_NEAR(above.d[1], below.d[1], 1e-9);
  // And the values are on opposite sides of the cut, i.e. this really is it.
  EXPECT_GT(std::abs(above.v - below.v), 6.0);
}

// -----------------------------------------------------------------------------
// Angle wrapping
// -----------------------------------------------------------------------------

// Expected values produced by the Python reference
// (combaero.network._mynard2010._wrap_to_pi / _wrap_to_2pi), NOT by Matlab's
// documented range. The two differ at exactly +pi: Matlab's wrapToPi is
// (-pi, pi] and returns +pi, while the Python expression
// `(x + pi) % (2*pi) - pi` returns -pi, so its real range is [-pi, pi). The
// port has to match the reference implementation, boundary included, or the
// equivalence gate fails on the one input where a lateral is exactly
// anti-parallel to the main duct. The Python docstring claimed Matlab's range
// and has been corrected.
TEST(DualNumber, WrapToPiMatchesThePythonReference) {
  struct { double in; double out; } cases[] = {
      {0.0, 0.0},         {M_PI, -M_PI},             {-M_PI, -M_PI},
      {3.0 * M_PI, -M_PI}, {1.5 * M_PI, -0.5 * M_PI}, {-1.5 * M_PI, 0.5 * M_PI},
      {2.0 * M_PI, 0.0},   {-0.5 * M_PI, -0.5 * M_PI},
  };
  for (const auto& c : cases) {
    EXPECT_NEAR(dwrap_to_pi(D1::seed(c.in, 0)).v, c.out, 1e-12) << "in=" << c.in;
  }
}

TEST(DualNumber, WrapTo2PiMatchesThePythonReference) {
  struct { double in; double out; } cases[] = {
      {0.0, 0.0},          {M_PI, M_PI},         {2.0 * M_PI, 0.0},
      {-0.5 * M_PI, 1.5 * M_PI}, {3.0 * M_PI, M_PI}, {-M_PI, M_PI},
      {-1.5 * M_PI, 0.5 * M_PI},
  };
  for (const auto& c : cases) {
    EXPECT_NEAR(dwrap_to_2pi(D1::seed(c.in, 0)).v, c.out, 1e-12) << "in=" << c.in;
  }
}

TEST(DualNumber, WrappingPassesPartialsThroughUnchanged) {
  // Both wraps shift by an integer multiple of 2*pi, so within a branch the
  // derivative is the identity. Pinned because zeroing the partials here would
  // silently drop the junction's phi angle from the Jacobian.
  for (double x : {-7.0, -1.0, 0.4, 3.5, 9.2}) {
    D2 a;
    a.v = x;
    a.d = {2.5, -0.75};
    for (const D2& wrapped : {dwrap_to_pi(a), dwrap_to_2pi(a)}) {
      EXPECT_DOUBLE_EQ(wrapped.d[0], 2.5) << "at x=" << x;
      EXPECT_DOUBLE_EQ(wrapped.d[1], -0.75) << "at x=" << x;
    }
  }
}

TEST(DualNumber, WrappingIsAShiftOfTheValueOnly) {
  for (double x : {-7.0, -1.0, 0.4, 3.5, 9.2}) {
    D1 a = D1::seed(x, 0);
    double two_pi = 2.0 * M_PI;
    double k_pi = (dwrap_to_pi(a).v - x) / two_pi;
    double k_2pi = (dwrap_to_2pi(a).v - x) / two_pi;
    EXPECT_NEAR(k_pi, std::round(k_pi), 1e-12) << "wrapToPi shifted by a non-integer at " << x;
    EXPECT_NEAR(k_2pi, std::round(k_2pi), 1e-12) << "wrapTo2Pi shifted by a non-integer at " << x;
  }
}
