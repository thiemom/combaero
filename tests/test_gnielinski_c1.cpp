#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "friction.h"
#include "heat_transfer.h"

using combaero::nusselt_gnielinski;
using combaero::nusselt_gnielinski_smooth_with_derivative;
using combaero::nusselt_gnielinski_with_derivative;

// #446. The clamp is the old max(Re, 3000) everywhere outside its C1 band,
// so nothing moves outside Re 2500-3500.
TEST(GnielinskiC1, ClampEqualsTheOldMaxOutsideItsBand) {
  for (double Re : {100.0, 1000.0, 2300.0, 2500.0, 3500.0, 4000.0, 1.0e5, 1.0e6}) {
    EXPECT_DOUBLE_EQ(friction_petukhov_clamped(Re),
                     friction_petukhov(std::max(Re, 3000.0)))
        << "Re=" << Re;
  }
}

TEST(GnielinskiC1, ClampSlopeIsExactAndContinuous) {
  for (double Re : {2600.0, 2999.0, 3000.0, 3001.0, 3400.0, 5000.0}) {
    const double h = 1e-3;
    const double fd = (friction_petukhov_clamped(Re + h) -
                       friction_petukhov_clamped(Re - h)) / (2.0 * h);
    EXPECT_NEAR(friction_petukhov_clamped_dRe(Re), fd, 1e-6 * std::abs(fd) + 1e-14)
        << "Re=" << Re;
  }
  for (double edge : {2500.0, 3500.0}) {
    EXPECT_NEAR(friction_petukhov_clamped_dRe(edge - 1e-9),
                friction_petukhov_clamped_dRe(edge + 1e-9), 1e-12)
        << "edge " << edge;
  }
}

// The analytic form IS nusselt_gnielinski, smooth-pipe and fixed-f alike.
TEST(GnielinskiC1, AnalyticValueMatchesTheLibraryFunction) {
  combaero::CorrelationStatus st;  // silences the below-range warning
  for (double Re : {-50.0, 0.0, 500.0, 1000.0, 1700.0, 2300.0, 2700.0, 3000.0,
                    3300.0, 1.0e4, 1.0e6}) {
    for (double Pr : {0.7, 5.0}) {
      const double smooth = nusselt_gnielinski(Re, Pr, &st);
      EXPECT_NEAR(nusselt_gnielinski_smooth_with_derivative(Re, Pr).Nu, smooth,
                  1e-12 * smooth) << "Re=" << Re;
      const double fixed = nusselt_gnielinski(Re, Pr, 0.05, &st);
      EXPECT_NEAR(nusselt_gnielinski_with_derivative(Re, Pr, 0.05).Nu, fixed,
                  1e-12 * fixed) << "Re=" << Re;
    }
  }
}

// The derivative everywhere, Re 3000 INCLUDED -- the point the old clamp
// broke -- and through the Hermite blend.
TEST(GnielinskiC1, AnalyticDerivativeMatchesCentralDifferences) {
  combaero::CorrelationStatus st;
  for (double Re : {1200.0, 2000.0, 2400.0, 2600.0, 2999.0, 3000.0, 3200.0,
                    3499.0, 5000.0, 1.0e5}) {
    const double h = Re * 1e-6;
    const double fd = (nusselt_gnielinski(Re + h, 0.7, &st) -
                       nusselt_gnielinski(Re - h, 0.7, &st)) / (2.0 * h);
    const double an = nusselt_gnielinski_smooth_with_derivative(Re, 0.7).dNu_dRe;
    EXPECT_LT(std::abs(an - fd) / std::max(std::abs(fd), 1e-12), 1e-5) << "Re=" << Re;
  }
}

// The kink #446 reported: the one-sided slopes at Re 3000 now agree.
TEST(GnielinskiC1, NoSlopeJumpAtRe3000) {
  const auto lo = nusselt_gnielinski_smooth_with_derivative(3000.0 - 1e-6, 0.7);
  const auto hi = nusselt_gnielinski_smooth_with_derivative(3000.0 + 1e-6, 0.7);
  EXPECT_NEAR(lo.dNu_dRe, hi.dNu_dRe, 1e-6 * hi.dNu_dRe);
}
