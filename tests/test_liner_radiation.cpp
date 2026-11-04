// Liner radiation kernels (#793): Lefebvre & Ballal Ch. 8 Eqs. 8.7-8.17.
// Every analytic derivative against central differences, and the local wall
// balance with radiation (wall_radiation_sides) against its own definition.

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "liner_radiation.h"
#include "thermo_transport_data.h"

using combaero::elemental_fuel_air_ratio;
using combaero::gas_emissivity_luminous;
using combaero::liner_casing_exchange_factor;
using combaero::liner_casing_radiation;
using combaero::liner_gas_radiation;
using combaero::luminosity_factor_hydrogen;
using combaero::mean_beam_length;
using combaero::STEFAN_BOLTZMANN;
using combaero::wall_radiation_sides;
using combaero::wall_adiabatic_face;
using combaero::WALL_ADIABATIC_N_PARAMS;
using combaero::WALL_RAD_N_PARAMS;

namespace {

double central(const std::function<double(double)> &f, double x, double rel) {
  const double h = rel * std::max(std::abs(x), 1.0);
  return (f(x + h) - f(x - h)) / (2.0 * h);
}

void expect_fd(double analytic, double fd, double rtol, const char *what) {
  EXPECT_NEAR(analytic, fd, rtol * std::max(std::abs(fd), 1e-12)) << what;
}

} // namespace

TEST(LinerRadiation, GasEmissivityDerivatives) {
  for (double Tg : {900.0, 1600.0, 2280.0}) {
    for (double P : {1e5, 3.04e6}) {
      for (double L : {1.0, 1.7, 3.0}) {
        const double q = 0.0588;
        const double lb = 0.0792;
        const auto r = gas_emissivity_luminous(Tg, P, q, lb, L);
        expect_fd(r.deps_dT,
                  central([&](double x) { return gas_emissivity_luminous(x, P, q, lb, L).eps; },
                          Tg, 1e-6),
                  1e-6, "T");
        expect_fd(r.deps_dP,
                  central([&](double x) { return gas_emissivity_luminous(Tg, x, q, lb, L).eps; },
                          P, 1e-6),
                  1e-6, "P");
        expect_fd(r.deps_dq,
                  central([&](double x) { return gas_emissivity_luminous(Tg, P, x, lb, L).eps; },
                          q, 1e-8),
                  1e-6, "q");
        expect_fd(r.deps_dl_b,
                  central([&](double x) { return gas_emissivity_luminous(Tg, P, q, x, L).eps; },
                          lb, 1e-8),
                  1e-6, "l_b");
        expect_fd(r.deps_dL,
                  central([&](double x) { return gas_emissivity_luminous(Tg, P, q, lb, x).eps; },
                          L, 1e-6),
                  1e-6, "L");
      }
    }
  }
  // q l_b = 0: no emission, finite partials.
  const auto z = gas_emissivity_luminous(2000.0, 1e6, 0.0, 0.1, 1.0);
  EXPECT_EQ(z.eps, 0.0);
  EXPECT_EQ(z.deps_dq, 0.0);
  EXPECT_THROW(gas_emissivity_luminous(0.0, 1e6, 0.05, 0.1, 1.0), std::invalid_argument);
  EXPECT_THROW(gas_emissivity_luminous(2000.0, -1.0, 0.05, 0.1, 1.0), std::invalid_argument);
}

TEST(LinerRadiation, LuminosityAndBeamLength) {
  const auto [L, dL] = luminosity_factor_hydrogen(13.9);
  EXPECT_NEAR(L, 336.0 / (13.9 * 13.9), 1e-15);
  expect_fd(dL,
            central([](double H) { return std::get<0>(luminosity_factor_hydrogen(H)); }, 13.9,
                    1e-7),
            1e-7, "dL/dH");
  // A long cylinder: V/A -> D/4, 3.4 V/A = 0.85 D (Lefebvre's 0.6-0.9 D_L).
  EXPECT_NEAR(mean_beam_length(0.25, 1.0), 0.85, 1e-15);
  EXPECT_NEAR(mean_beam_length(0.25, 1.0, combaero::BEAM_LENGTH_COEFF_HOWELL), 0.9, 1e-15);
  EXPECT_THROW(luminosity_factor_hydrogen(0.0), std::invalid_argument);
  EXPECT_THROW(mean_beam_length(1.0, 0.0), std::invalid_argument);
}

TEST(LinerRadiation, R1AndR2Derivatives) {
  for (double Tw : {700.0, 1640.0, 2400.0}) {
    const double Tg = 2280.0;
    const double eg = 0.61;
    const double ew = 0.7;
    const auto r = liner_gas_radiation(Tg, Tw, eg, ew);
    expect_fd(r.dR1_dT_g,
              central([&](double x) { return liner_gas_radiation(x, Tw, eg, ew).R1; }, Tg, 1e-7),
              1e-6, "dR1/dTg");
    expect_fd(r.dR1_dT_w,
              central([&](double x) { return liner_gas_radiation(Tg, x, eg, ew).R1; }, Tw, 1e-7),
              1e-6, "dR1/dTw");
    expect_fd(r.dR1_deps_g,
              central([&](double x) { return liner_gas_radiation(Tg, Tw, x, ew).R1; }, eg, 1e-7),
              1e-6, "dR1/deps_g");
    expect_fd(r.dR1_deps_w,
              central([&](double x) { return liner_gas_radiation(Tg, Tw, eg, x).R1; }, ew, 1e-7),
              1e-6, "dR1/deps_w");
    const auto [R2, dTw, dTc] = liner_casing_radiation(Tw, 880.0, ew, 0.4, 0.7);
    (void)R2;
    expect_fd(dTw,
              central([&](double x) { return std::get<0>(liner_casing_radiation(x, 880.0, ew, 0.4, 0.7)); },
                      Tw, 1e-7),
              1e-6, "dR2/dTw");
    expect_fd(dTc,
              central([&](double x) { return std::get<0>(liner_casing_radiation(Tw, x, ew, 0.4, 0.7)); },
                      880.0, 1e-7),
              1e-6, "dR2/dTc");
  }
  // Equal temperatures exchange nothing; black parallel surfaces give sigma.
  // (to round-off: GCC may contract Tg15 Tg - Tw15 Tw into an FMA)
  EXPECT_NEAR(liner_gas_radiation(1500.0, 1500.0, 0.5, 0.7).R1, 0.0,
              1e-12 * liner_gas_radiation(1500.0, 0.0, 0.5, 0.7).R1);
  EXPECT_NEAR(liner_casing_exchange_factor(1.0, 1.0, 0.7), STEFAN_BOLTZMANN, 1e-22);
  EXPECT_THROW(liner_casing_radiation(1000.0, 800.0, 0.0, 0.4, 0.7), std::invalid_argument);
}

TEST(LinerRadiation, ElementalFuelAirRatio) {
  // Products of n-heptane burnt completely at q = 0.05 hold the same C and H:
  // build them by mass balance from the reactants and compare.
  const std::size_t n = combaero::molar_masses.size();
  const auto idx = [](const char *s) { return static_cast<std::size_t>(combaero::species_index.at(s)); };
  std::vector<double> Y(n, 0.0);
  const double q = 0.05;
  Y[idx("NC7H16")] = q / (1.0 + q);
  Y[idx("O2")] = 0.232 / (1.0 + q);
  Y[idx("N2")] = 0.768 / (1.0 + q);
  const auto [q_r, dq] = elemental_fuel_air_ratio(Y);
  EXPECT_NEAR(q_r, q, 1e-12);
  for (std::size_t k = 0; k < n; ++k) {
    auto f = [&](double x) {
      auto Yp = Y;
      Yp[k] = x;
      return std::get<0>(elemental_fuel_air_ratio(Yp));
    };
    expect_fd(dq[k], central(f, Y[k], 1e-7), 1e-6, "dq/dY");
  }
  EXPECT_THROW(elemental_fuel_air_ratio(std::vector<double>(3, 0.0)), std::invalid_argument);
}

namespace {

std::array<double, WALL_RAD_N_PARAMS> book_params(double T_aw_a, double h_a) {
  const double eg = gas_emissivity_luminous(2280.0, 3040e3, 0.0588, 0.0792, 1.7).eps;
  const double k_a = 0.5 * STEFAN_BOLTZMANN * 1.7 * eg;
  const double k_b = liner_casing_exchange_factor(0.7, 0.4, 0.7);
  return {T_aw_a, h_a, 880.0, 921.0, 0.0012 / 26.0, 2280.0, k_a, 880.0, k_b};
}

combaero::WallRadiationSides eval(const std::array<double, WALL_RAD_N_PARAMS> &p, int ka,
                                  int kb) {
  return wall_radiation_sides(p[0], p[1], p[2], p[3], p[4], ka, p[5], p[6], kb, p[7], p[8]);
}

} // namespace

TEST(LinerRadiation, WallSidesSatisfyTheBalance) {
  const auto p = book_params(2280.0, 565.0);
  const auto w = eval(p, 1, 2);
  // Each face's flux is its convection plus Eq. 8.7 / Eq. 8.17 exactly.
  const double eg = gas_emissivity_luminous(2280.0, 3040e3, 0.0588, 0.0792, 1.7).eps;
  const double R1 = liner_gas_radiation(2280.0, w.T_s_a, eg, 0.7).R1;
  const double R2 = std::get<0>(liner_casing_radiation(w.T_s_b, 880.0, 0.7, 0.4, 0.7));
  const double q_hot = 565.0 * (2280.0 - w.T_s_a) + R1;
  const double q_cold = 921.0 * (w.T_s_b - 880.0) + R2;
  const double q_wall = (w.T_s_a - w.T_s_b) / p[4];
  EXPECT_NEAR(q_hot, w.q, 1e-9 * w.q);
  EXPECT_NEAR(q_cold, w.q, 1e-9 * w.q);
  EXPECT_NEAR(q_wall, w.q, 1e-9 * w.q);
  EXPECT_NEAR(w.q_rad_a, R1, 1e-9 * R1);
  EXPECT_NEAR(w.q_rad_b, R2, 1e-9 * R2);
  // Effective sides: H (T* - T_s) is the face's total flux.
  EXPECT_NEAR(w.H_a * (w.T_a - w.T_s_a), w.q, 1e-9 * w.q);
  EXPECT_NEAR(w.H_b * (w.T_s_b - w.T_b), w.q, 1e-9 * w.q);
}

TEST(LinerRadiation, WallSidesDerivatives) {
  for (int ka : {0, 1, 2}) {
    for (int kb : {0, 1, 2}) {
      auto p = book_params(2280.0, 565.0);
      if (kb == 1) {
        p[7] = 1900.0;  // a second gas on side B, colder than side A's
      }
      const auto w = eval(p, ka, kb);
      // Scaled derivatives p_j dy/dp_j against relative central differences:
      // k is O(1e-8), so an absolute step drowns in roundoff.
      const std::array<double, 4> y0{w.H_a, w.T_a, w.H_b, w.T_b};
      for (int j = 0; j < WALL_RAD_N_PARAMS; ++j) {
        const double e = 1e-6;
        const double s = p[j] != 0.0 ? p[j] : 1.0;
        auto pp = p;
        auto pm = p;
        pp[j] += e * s;
        pm[j] -= e * s;
        const auto wp = eval(pp, ka, kb);
        const auto wm = eval(pm, ka, kb);
        const std::array<double, 4> yp{wp.H_a, wp.T_a, wp.H_b, wp.T_b};
        const std::array<double, 4> ym{wm.H_a, wm.T_a, wm.H_b, wm.T_b};
        for (int i = 0; i < 4; ++i) {
          const double fd = (yp[i] - ym[i]) / (2.0 * e);
          EXPECT_NEAR(w.d[i][j] * s, fd, 1e-6 * std::max(std::abs(fd), 1e-3 * std::abs(y0[i])))
              << "kinds " << ka << kb << " out " << i << " param " << j;
        }
        const double fa = (wp.q_rad_a - wm.q_rad_a) / (2.0 * e);
        const double fb = (wp.q_rad_b - wm.q_rad_b) / (2.0 * e);
        EXPECT_NEAR(w.dq_rad[0][j] * s, fa, 1e-6 * std::max(std::abs(fa), 1e-3 * w.q)) << j;
        EXPECT_NEAR(w.dq_rad[1][j] * s, fb, 1e-6 * std::max(std::abs(fb), 1e-3 * w.q)) << j;
      }
    }
  }
}

TEST(LinerRadiation, WallSidesWithoutRadiationAreTheConvectiveSides) {
  const auto p = book_params(2280.0, 565.0);
  const auto w = eval(p, 0, 0);
  EXPECT_EQ(w.H_a, 565.0);
  EXPECT_EQ(w.T_a, 2280.0);
  EXPECT_EQ(w.H_b, 921.0);
  EXPECT_EQ(w.T_b, 880.0);
  EXPECT_EQ(w.q_rad_a, 0.0);
  EXPECT_EQ(w.iterations, 0);
}

// Fidelity, printed quantities: Lefebvre & Ballal 3rd ed. pp. 325-328
// (uncooled), 338-340 (film-cooled, x/s = 18) and 350-351 (0.5 mm TBC),
// solved from the printed inputs with the printed convective coefficients:
// C1 = 565 (Tg - Tw1) (Eq. 8.28), C2 = 921 (Tw2 - 880) (Eq. 8.29, which
// prints 291 for 921), film C1 = 1926 (1176 - Tw1). The book rounds eps_g to
// 0.61 and C1 to 1 280 500 - 562 Tw1 before solving; from the inputs the
// uncooled liner comes out 1.6 K / 1.4 K hotter than printed.
TEST(LinerRadiation, LefebvreWorkedExamples) {
  const auto unc = eval(book_params(2280.0, 565.0), 1, 2);
  EXPECT_NEAR(unc.T_s_a, 1640.0, 2.0);
  EXPECT_NEAR(unc.T_s_b, 1603.0, 2.0);
  const auto film = eval(book_params(1176.0, 1926.0), 1, 2);
  EXPECT_NEAR(film.T_s_a, 1283.0, 1.0);
  EXPECT_NEAR(film.T_s_b, 1265.0, 1.0);
  auto p = book_params(1176.0, 1926.0);
  const double r_tbc = 0.0005 / 2.6;
  p[4] += r_tbc;
  const auto tbc = eval(p, 1, 2);
  EXPECT_NEAR(tbc.T_s_a, 1304.0, 1.0);
  EXPECT_NEAR(tbc.T_s_a - tbc.q * r_tbc, 1236.0, 1.0);
  EXPECT_NEAR(tbc.T_s_b, 1220.0, 1.0);
}

// A face whose other side is adiabatic (#970): the face's own balance.
TEST(LinerRadiation, AdiabaticFaceSatisfiesItsBalance) {
  const double eg = gas_emissivity_luminous(2280.0, 3040e3, 0.0588, 0.0792, 1.7).eps;
  const double k_gas = 0.5 * STEFAN_BOLTZMANN * 1.7 * eg;
  // Gas radiation into a face at T_aw 1500 K from a 2280 K gas: the face
  // runs hot of T_aw, and its convection hands back exactly what it absorbs.
  const auto g = wall_adiabatic_face(1500.0, 565.0, 1, 2280.0, k_gas);
  const double R1 = liner_gas_radiation(2280.0, g.T_s, eg, 0.7).R1;
  EXPECT_GT(g.T_s, 1500.0);
  EXPECT_LT(g.T_s, 2280.0);
  EXPECT_NEAR(565.0 * (1500.0 - g.T_s) + R1, 0.0, 1e-9 * R1);
  EXPECT_NEAR(g.q_rad, R1, 1e-9 * R1);
  // A surface face radiating to a colder casing runs cold of T_aw.
  const double F = liner_casing_exchange_factor(0.7, 0.4, 0.7);
  const auto c = wall_adiabatic_face(900.0, 300.0, 2, 600.0, F);
  const double R2 = std::get<0>(liner_casing_radiation(c.T_s, 600.0, 0.7, 0.4, 0.7));
  EXPECT_LT(c.T_s, 900.0);
  EXPECT_GT(c.T_s, 600.0);
  EXPECT_NEAR(300.0 * (900.0 - c.T_s) - R2, 0.0, 1e-9 * R2);
  // No radiation: the face sits at its adiabatic wall temperature.
  const auto n = wall_adiabatic_face(812.5, 40.0, 0, 0.0, 0.0);
  EXPECT_EQ(n.T_s, 812.5);
  EXPECT_EQ(n.q_rad, 0.0);
  EXPECT_EQ(n.dT_s[0], 1.0);
  EXPECT_THROW(wall_adiabatic_face(800.0, 0.0, 1, 2000.0, k_gas), std::invalid_argument);
  EXPECT_THROW(wall_adiabatic_face(800.0, 10.0, 3, 2000.0, k_gas), std::invalid_argument);
}

TEST(LinerRadiation, AdiabaticFaceDerivatives) {
  const double eg = gas_emissivity_luminous(2280.0, 3040e3, 0.0588, 0.0792, 1.7).eps;
  const double k_gas = 0.5 * STEFAN_BOLTZMANN * 1.7 * eg;
  const double F = liner_casing_exchange_factor(0.7, 0.4, 0.7);
  for (int kind : {1, 2}) {
    const std::array<double, WALL_ADIABATIC_N_PARAMS> p =
        kind == 1 ? std::array<double, WALL_ADIABATIC_N_PARAMS>{1500.0, 565.0, 2280.0, k_gas}
                  : std::array<double, WALL_ADIABATIC_N_PARAMS>{900.0, 300.0, 600.0, F};
    const auto w = wall_adiabatic_face(p[0], p[1], kind, p[2], p[3]);
    for (int j = 0; j < WALL_ADIABATIC_N_PARAMS; ++j) {
      // p_j dT_s/dp_j against a relative central difference (k is O(1e-8)).
      const double e = 1e-6;
      auto pp = p;
      auto pm = p;
      pp[static_cast<std::size_t>(j)] *= 1.0 + e;
      pm[static_cast<std::size_t>(j)] *= 1.0 - e;
      const double fd = (wall_adiabatic_face(pp[0], pp[1], kind, pp[2], pp[3]).T_s -
                         wall_adiabatic_face(pm[0], pm[1], kind, pm[2], pm[3]).T_s) /
                        (2.0 * e);
      const double an = w.dT_s[static_cast<std::size_t>(j)] * p[static_cast<std::size_t>(j)];
      EXPECT_NEAR(an, fd, 1e-6 * std::max(std::abs(fd), 1e-3 * w.T_s))
          << "kind " << kind << " param " << j;
    }
  }
}
