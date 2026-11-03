// Film-side h downstream of a slot (#962): Seban (1960) Fig. 12 as drawn,
// and the Ballal & Lefebvre (1973) wall-jet option. The curves against
// their closed forms, the B/C blend's continuity (value and slope), the
// element mean against a fine quadrature, the hand-back, the book's printed
// 1926 W/m^2K, and every analytic derivative against central differences.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "humidair.h"
#include "slot_film_htc.h"

using combaero::slot_film_htc_ballal_lefebvre_1973;
using combaero::slot_film_htc_seban_1960;
using combaero::slot_film_nusselt_ballal_lefebvre_1973;
using combaero::slot_film_stanton_seban_1960;
using combaero::SlotFilmHtc;
namespace sb = combaero::seban1960;
namespace bl = combaero::ballal_lefebvre1973;

namespace {

double central(const std::function<double(double)> &f, double x, double rel) {
  const double h = rel * std::max(std::abs(x), 1e-12);
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

void expect_fd(double analytic, double fd, double rtol, const char *what) {
  EXPECT_NEAR(analytic, fd, rtol * std::max(std::abs(fd), 1e-12)) << what;
}

double point_St(double X, double Re) {
  return slot_film_stanton_seban_1960(X, X, Re).St;
}

} // namespace

TEST(SlotFilmHtc, SebanCurvesAsDrawn) {
  const double Re = 5000.0;
  // Away from the blend: the closed forms exactly.
  for (double X : {4.0, 8.0, 20.0, 29.0}) {
    EXPECT_NEAR(point_St(X, Re),
                sb::C_near * std::pow(Re, -0.3) * std::pow(X, -0.3),
                1e-14 * point_St(X, Re));
  }
  for (double X : {48.0, 100.0, 300.0}) {
    EXPECT_NEAR(point_St(X, Re),
                sb::C_far * std::pow(Re, -0.3) * std::pow(X, -0.6),
                1e-14 * point_St(X, Re));
  }
  // The crossing: both curves meet at x/s 37.4, where the blend is exact.
  const double Xs = std::pow(sb::C_far / sb::C_near, 1.0 / 0.3);
  EXPECT_NEAR(Xs, 37.36, 0.01);
  EXPECT_NEAR(point_St(Xs, Re),
              sb::C_near * std::pow(Re, -0.3) * std::pow(Xs, -0.3),
              1e-12 * point_St(Xs, Re));
  // Within the blend it stays within 1.5 % of min(B, C).
  for (double X = 29.0; X < 48.0; X += 0.25) {
    const double B = sb::C_near * std::pow(X, -0.3);
    const double C = sb::C_far * std::pow(X, -0.6);
    const double f = point_St(X, Re) * std::pow(Re, 0.3);
    EXPECT_LE(std::abs(f / std::min(B, C) - 1.0), 0.015) << X;
  }
  // The printed constants reproduce Eqs. 13/14 as printed.
  EXPECT_NEAR(slot_film_stanton_seban_1960(10.0, 10.0, Re, sb::C_near_printed,
                                           sb::C_far_printed)
                  .St,
              0.10 * std::pow(Re, -0.3) * std::pow(10.0, -0.3), 1e-15);
}

TEST(SlotFilmHtc, SebanBlendIsC1) {
  // Value and slope of ln f against ln X continuous at both blend edges.
  const double Xs = std::pow(sb::C_far / sb::C_near, 1.0 / 0.3);
  const double d = sb::blend_half_width_decades * std::log(10.0);
  const auto lnf = [](double u) { return std::log(point_St(std::exp(u), 1.0)); };
  for (double ue : {std::log(Xs) - d, std::log(Xs) + d}) {
    const double e = 1e-7;
    EXPECT_NEAR(lnf(ue - e), lnf(ue + e), 1e-6);
    const double sl = (lnf(ue - e) - lnf(ue - 3 * e)) / (2 * e);
    const double sr = (lnf(ue + 3 * e) - lnf(ue + e)) / (2 * e);
    EXPECT_NEAR(sl, sr, 1e-5) << ue;
  }
}

TEST(SlotFilmHtc, SpanMeanIsTheIntegral) {
  const double Re = 3000.0;
  for (auto [a, b] : std::vector<std::pair<double, double>>{
           {0.0, 20.0}, {5.0, 33.0}, {20.0, 45.0}, {30.0, 40.0},
           {40.0, 200.0}, {10.0, 120.0}}) {
    // Fine midpoint rule in ln X (a > 0) or X.
    const int n = 200000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double X = a + (b - a) * (i + 0.5) / n;
      sum += point_St(X, Re);
    }
    const double mean = sum / n;
    EXPECT_NEAR(slot_film_stanton_seban_1960(a, b, Re).St, mean,
                (a == 0.0 ? 2e-4 : 1e-8) * mean)
        << a << " " << b;
  }
}

TEST(SlotFilmHtc, FarEndHandsBackToTheDuct) {
  const double Re = 3000.0;
  const auto r = slot_film_stanton_seban_1960(200.0, 400.0, Re);
  EXPECT_DOUBLE_EQ(r.fraction_in_box, 0.5);
  EXPECT_FALSE(r.x_in_box);
  EXPECT_NEAR(r.St, slot_film_stanton_seban_1960(200.0, 300.0, Re).St, 1e-15);
  const auto beyond = slot_film_stanton_seban_1960(310.0, 400.0, Re);
  EXPECT_EQ(beyond.fraction_in_box, 0.0);
  EXPECT_EQ(beyond.branch, combaero::SLOT_FILM_BRANCH_NONE);
  // Dimensional: past the far end the side's h is the duct's.
  const auto X = combaero::dry_air();
  const double s = 0.002;
  const SlotFilmHtc h = slot_film_htc_seban_1960(0.7, 0.8, s, 60.0, 2.0, 150.0,
                                                 600.0, 1.0e6, X);
  EXPECT_DOUBLE_EQ(h.h, 150.0);
  EXPECT_DOUBLE_EQ(h.dh_dh_duct, 1.0);
  // Half the span in the box: the length-weighted mix.
  const SlotFilmHtc g = slot_film_htc_seban_1960(0.5, 0.7, s, 60.0, 2.0, 150.0,
                                                 600.0, 1.0e6, X);
  EXPECT_NEAR(g.h, 0.5 * g.h_film + 0.5 * 150.0, 1e-10 * g.h);
}

TEST(SlotFilmHtc, HandsBackBelowUnitMassVelocityRatio) {
  const auto X = combaero::dry_air();
  const double s = 0.002, h_duct = 150.0;
  const auto at = [&](double M) {
    return slot_film_htc_seban_1960(0.0, 0.04, s, 60.0, M, h_duct, 600.0,
                                    1.0e6, X);
  };
  EXPECT_DOUBLE_EQ(at(0.9).h, h_duct);
  EXPECT_DOUBLE_EQ(at(1.0).h, h_duct);
  EXPECT_FALSE(at(0.9).M_in_box);
  const SlotFilmHtc full = at(1.2);
  EXPECT_DOUBLE_EQ(full.h, full.h_film);
  EXPECT_DOUBLE_EQ(full.weight_M, 1.0);
  EXPECT_TRUE(full.M_in_box);
  EXPECT_NEAR(at(1.05).h, 0.5 * (full.h_film + h_duct), 1e-10 * full.h);
  EXPECT_FALSE(at(9.5).M_in_box);
  // The wall jet hands back below 1.3.
  const auto wj = slot_film_htc_ballal_lefebvre_1973(0.0, 0.04, s, 60.0, 1.25,
                                                     h_duct, 600.0, 1.0e6, X);
  EXPECT_DOUBLE_EQ(wj.h, h_duct);
}

TEST(SlotFilmHtc, BookWallJetCoefficient) {
  // Lefebvre & Ballal p. 340: C1 = 0.10 (k_a/x) Re_x^0.8 (x/s)^-0.36 with
  // the book's k_a 0.0553, Re_s 1.81e4, x/s 18, s 1.45 mm: 1926 W/m^2K.
  const double s = 0.00145;
  const auto n = slot_film_nusselt_ballal_lefebvre_1973(18.0, 18.0, 1.81e4);
  EXPECT_NEAR(n.Nu_s * 0.0553 / s, 1926.0, 1.0);
  EXPECT_EQ(n.branch, combaero::SLOT_FILM_BRANCH_WALL_JET);
  // Its span mean is the integral of the point values.
  const int N = 200000;
  double sum = 0.0;
  for (int i = 0; i < N; ++i) {
    const double Xi = 5.0 + 25.0 * (i + 0.5) / N;
    sum += slot_film_nusselt_ballal_lefebvre_1973(Xi, Xi, 1.81e4).Nu_s;
  }
  EXPECT_NEAR(slot_film_nusselt_ballal_lefebvre_1973(5.0, 30.0, 1.81e4).Nu_s,
              sum / N, 1e-8 * sum / N);
}

TEST(SlotFilmHtc, DimensionalDerivatives) {
  const auto X = combaero::dry_air();
  const double s = 0.0015;
  using Kernel = SlotFilmHtc (*)(double, double, double, double, double,
                                 double, double, double,
                                 const std::vector<double> &);
  for (Kernel k : {Kernel(&slot_film_htc_seban_1960),
                   Kernel(&slot_film_htc_ballal_lefebvre_1973)}) {
    // M inside the hand-back band (1.05 / 1.35) and above it (2.5).
    for (double M : {1.05, 1.35, 2.5}) {
      for (const auto &ab : std::vector<std::pair<double, double>>{
               {0.0, 0.03}, {0.02, 0.09}, {0.05, 0.6}}) {
        const double a = ab.first, b = ab.second;
        const double G = 400.0, T = 820.0, hd = 600.0, P = 2.0e6;
        const SlotFilmHtc r = k(a, b, s, G, M, hd, T, P, X);
        expect_fd(r.dh_dG_s,
                  central([&](double g) { return k(a, b, s, g, M, hd, T, P, X).h; },
                          G, 1e-6),
                  1e-6, "dh_dG_s");
        expect_fd(r.dh_dT_s,
                  central([&](double t) { return k(a, b, s, G, M, hd, t, P, X).h; },
                          T, 1e-6),
                  1e-6, "dh_dT_s");
        expect_fd(r.dh_dM,
                  central([&](double m) { return k(a, b, s, G, m, hd, T, P, X).h; },
                          M, 1e-7),
                  1e-6, "dh_dM");
        // Dimensionless and often exactly 0 (fully filmed): absolute.
        EXPECT_NEAR(r.dh_dh_duct,
                    central([&](double d) { return k(a, b, s, G, M, d, T, P, X).h; },
                            hd, 1e-6),
                    1e-7)
            << "dh_dh_duct";
      }
    }
  }
}

TEST(SlotFilmHtc, Refusals) {
  const auto X = combaero::dry_air();
  EXPECT_THROW(slot_film_stanton_seban_1960(0.0, 0.0, 3000.0),
               std::invalid_argument);
  EXPECT_THROW(slot_film_stanton_seban_1960(10.0, 5.0, 3000.0),
               std::invalid_argument);
  EXPECT_THROW(slot_film_stanton_seban_1960(1.0, 5.0, 0.0),
               std::invalid_argument);
  EXPECT_THROW(slot_film_htc_seban_1960(0.0, 0.01, 0.0, 60.0, 2.0, 100.0,
                                        600.0, 1e6, X),
               std::invalid_argument);
  EXPECT_THROW(slot_film_htc_seban_1960(0.0, 0.01, 0.001, 0.0, 2.0, 100.0,
                                        600.0, 1e6, X),
               std::invalid_argument);
}
