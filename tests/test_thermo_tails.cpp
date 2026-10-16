// The species thermo beyond its NASA-9 data (#606): the shipped data is what
// it claims (NASA-9, the pinned ranges), cp stays positive and bounded at
// every T, the tails join the data C2, h and s are cp's integrals on them, and
// the inverses recover T from h, s, u anywhere -- under both thermo models.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "thermo.h"
#include "thermo_transport_data.h"

using combaero::ThermoModel;

namespace {

struct ModelGuard {
  explicit ModelGuard(ThermoModel m) { combaero::set_thermo_model(m); }
  ~ModelGuard() { combaero::set_thermo_model(ThermoModel::NASA9); }
};

const ThermoModel kModels[] = {ThermoModel::NASA9, ThermoModel::NASA9_C2};

double H_R(std::size_t i, double T) { return combaero::h_species(i, T) / combaero::thermo::R_GAS; }

double T_lo(std::size_t i) { return combaero::nasa_coeffs[i].intervals.front().T_min; }
double T_hi(std::size_t i) { return combaero::nasa_coeffs[i].intervals.back().T_max; }

// Log-spaced from 1 K to 100 000 K.
std::vector<double> sweep() {
  std::vector<double> Ts;
  for (int k = 0; k <= 500; ++k) {
    Ts.push_back(std::pow(10.0, 5.0 * k / 500.0));
  }
  return Ts;
}

std::vector<double> pure(std::size_t i) {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[i] = 1.0;
  return X;
}

}  // namespace

// The data must be the CEA NASA-9 set: every species' ranges as pinned, and at
// least one interval using the T^-2 or T^-1 term NASA-7 does not have (a
// NASA-7 species written into the NASA-9 struct has both zero everywhere).
// AR's 200-1000 K interval is genuinely constant cp in CEA.
TEST(ThermoTails, TheDataIsNasa9OverThePinnedRanges) {
  struct Pin {
    const char* name;
    std::size_t n;
    double T_hi;
  };
  const Pin pins[] = {
      {"N2", 3, 20000.0},    {"O2", 3, 20000.0},     {"AR", 3, 20000.0},
      {"CO2", 3, 20000.0},   {"H2O", 2, 6000.0},     {"CH4", 2, 6000.0},
      {"C2H6", 2, 6000.0},   {"C3H8", 2, 6000.0},    {"IC4H10", 2, 6000.0},
      {"NC5H12", 2, 6000.0}, {"NC6H14", 2, 6000.0},  {"NC7H16", 2, 6000.0},
      {"CO", 3, 20000.0},    {"H2", 3, 20000.0},     {"NH3", 2, 6000.0},
  };
  ASSERT_EQ(combaero::num_species(), std::size(pins));
  for (const Pin& p : pins) {
    const std::size_t i = combaero::species_index_from_name(p.name);
    const auto& ivs = combaero::nasa_coeffs[i].intervals;
    ASSERT_EQ(ivs.size(), p.n) << p.name;
    EXPECT_EQ(ivs.front().T_min, 200.0) << p.name;
    EXPECT_EQ(ivs[0].T_max, 1000.0) << p.name;
    EXPECT_EQ(ivs[1].T_max, 6000.0) << p.name;
    EXPECT_EQ(ivs.back().T_max, p.T_hi) << p.name;
    for (std::size_t k = 1; k < ivs.size(); ++k) {
      EXPECT_EQ(ivs[k].T_min, ivs[k - 1].T_max) << p.name;
    }
    const bool nasa9_terms = std::any_of(ivs.begin(), ivs.end(), [](const auto& iv) {
      return iv.coeffs[0] != 0.0 || iv.coeffs[1] != 0.0;
    });
    EXPECT_TRUE(nasa9_terms) << p.name << " has no T^-2 / T^-1 term: NASA-7 data?";
  }
}

// cp > 0 and finite at every T, and on the tails between its value at the
// data end and its constant: never the polynomial's runaway (cp/R -25 for
// CH4 at 50 K, -241 for H2O at 20 000 K before #606).
TEST(ThermoTails, CpIsPositiveAndBoundedAtEveryT) {
  for (ThermoModel m : kModels) {
    ModelGuard g(m);
    for (std::size_t i = 0; i < combaero::num_species(); ++i) {
      const double cp_lo = combaero::cp_R(i, T_lo(i));
      const double cp_hi = combaero::cp_R(i, T_hi(i));
      const double cp_lo_inf = combaero::cp_R(i, 1.0);
      const double cp_hi_inf = combaero::cp_R(i, 1e5);
      for (double T : sweep()) {
        const double cp = combaero::cp_R(i, T);
        ASSERT_TRUE(std::isfinite(cp)) << i << " " << T;
        EXPECT_GT(cp, 2.0) << i << " " << T;
        if (T < T_lo(i)) {
          EXPECT_GE(cp, std::min(cp_lo, cp_lo_inf) - 1e-12) << i << " " << T;
          EXPECT_LE(cp, std::max(cp_lo, cp_lo_inf) + 1e-12) << i << " " << T;
        } else if (T > T_hi(i)) {
          EXPECT_GE(cp, std::min(cp_hi, cp_hi_inf) - 1e-12) << i << " " << T;
          EXPECT_LE(cp, std::max(cp_hi, cp_hi_inf) + 1e-12) << i << " " << T;
        }
      }
    }
  }
}

// C2 where the tails meet the data and where their fades end: cp, dcp/dT
// continuous, d2cp/dT2 continuous to its one-sided difference's truncation.
TEST(ThermoTails, TheTailsJoinTheDataC2) {
  using combaero::thermo::TAIL_FRACTION_HIGH;
  using combaero::thermo::TAIL_FRACTION_LOW;
  for (ThermoModel m : kModels) {
    ModelGuard g(m);
    for (std::size_t i = 0; i < combaero::num_species(); ++i) {
      const double edges[] = {T_lo(i), T_lo(i) * (1.0 - TAIL_FRACTION_LOW), T_hi(i),
                              T_hi(i) * (1.0 + TAIL_FRACTION_HIGH)};
      for (double E : edges) {
        const double e = 1e-9 * E;
        const double cp = combaero::cp_R(i, E);
        const double d_in = combaero::dcp_R_dT(i, E - e);
        const double d_out = combaero::dcp_R_dT(i, E + e);
        // The probes' own separation, 2e, taken out with the slope.
        EXPECT_NEAR(combaero::cp_R(i, E + e) - combaero::cp_R(i, E - e), (d_in + d_out) * e,
                    1e-12 * cp)
            << i << " " << E;
        // dcp/dT as its effect on cp over 1 % of T.
        EXPECT_LT(std::abs(d_out - d_in) * (0.01 * E) / cp, 1e-8) << i << " " << E;
        const double h = 1e-6 * E;
        const double c2_in = (combaero::dcp_R_dT(i, E - h) - combaero::dcp_R_dT(i, E - 2 * h)) / h;
        const double c2_out = (combaero::dcp_R_dT(i, E + 2 * h) - combaero::dcp_R_dT(i, E + h)) / h;
        EXPECT_LT(std::abs(c2_out - c2_in) * (0.01 * E) * (0.01 * E) / cp, 1e-6) << i << " " << E;
      }
    }
  }
}

// h' = cp and T s' = cp on the tails (the fade and beyond), and h, s
// continuous at the data ends: the integrals are exact, not re-based.
TEST(ThermoTails, EnthalpyAndEntropyAreTheIntegralsOfCpOnTheTails) {
  for (ThermoModel m : kModels) {
    ModelGuard g(m);
    for (std::size_t i = 0; i < combaero::num_species(); ++i) {
      for (double T : {20.0, 100.0, 160.0, 190.0, 199.9}) {
        const double e = 1e-4 * T;
        const double cp = combaero::cp_R(i, T);
        EXPECT_NEAR((H_R(i, T + e) - H_R(i, T - e)) / (2 * e), cp, 1e-6 * cp) << i << " " << T;
        EXPECT_NEAR((combaero::s_R(i, T + e) - combaero::s_R(i, T - e)) / (2 * e) * T, cp,
                    1e-6 * cp)
            << i << " " << T;
      }
      for (double f : {1.0001, 1.05, 1.099, 1.2, 3.0}) {
        const double T = f * T_hi(i);
        const double e = 1e-5 * T;
        const double cp = combaero::cp_R(i, T);
        EXPECT_NEAR((H_R(i, T + e) - H_R(i, T - e)) / (2 * e), cp, 1e-6 * cp) << i << " " << T;
        EXPECT_NEAR((combaero::s_R(i, T + e) - combaero::s_R(i, T - e)) / (2 * e) * T, cp,
                    1e-6 * cp)
            << i << " " << T;
      }
      for (double E : {T_lo(i), T_hi(i)}) {
        const double e = 1e-9 * E;
        const double cp = combaero::cp_R(i, E);
        // Less the probes' own separation 2e (slopes cp and cp / T).
        EXPECT_NEAR(H_R(i, E + e) - H_R(i, E - e), 2 * e * cp, 1e-12 * cp * E) << i << " " << E;
        EXPECT_NEAR(combaero::s_R(i, E + e) - combaero::s_R(i, E - e), 2 * e * cp / E, 1e-12 * cp)
            << i << " " << E;
      }
    }
  }
}

// Recovers T from h, s and u from a 300 K guess, inside the data and far
// outside it, for every pure species and a burnt mixture.
TEST(ThermoTails, TheInversesRecoverTAnywhere) {
  std::vector<std::vector<double>> Xs;
  for (std::size_t i = 0; i < combaero::num_species(); ++i) {
    Xs.push_back(pure(i));
  }
  std::vector<double> burnt(combaero::num_species(), 0.0);
  burnt[combaero::species_index_from_name("N2")] = 0.71;
  burnt[combaero::species_index_from_name("H2O")] = 0.19;
  burnt[combaero::species_index_from_name("CO2")] = 0.09;
  burnt[combaero::species_index_from_name("O2")] = 0.01;
  Xs.push_back(burnt);
  const double P = 1.0e5;
  for (ThermoModel m : kModels) {
    ModelGuard g(m);
    for (const auto& X : Xs) {
      for (double T : {20.0, 120.0, 199.0, 300.0, 1000.0, 2500.0, 6500.0, 9000.0, 30000.0, 60000.0}) {
        EXPECT_NEAR(combaero::calc_T_from_h(combaero::h(T, X), X, 300.0), T, 1e-6 * T) << T;
        EXPECT_NEAR(combaero::calc_T_from_s(combaero::s(T, X, P), P, X, 300.0), T, 1e-6 * T) << T;
        EXPECT_NEAR(combaero::calc_T_from_u(combaero::u(T, X), X, 300.0), T, 1e-6 * T) << T;
      }
    }
  }
}

// An h no positive T reaches (below h at 0 K on the constant tail) is a
// failure the caller hears about, not a silently wrong T.
TEST(ThermoTails, AnUnreachableTargetThrows) {
  const auto X = pure(combaero::species_index_from_name("N2"));
  const double h_unreachable = combaero::h(1.0, X) - 1e6;
  EXPECT_THROW(combaero::calc_T_from_h(h_unreachable, X, 300.0), std::runtime_error);
}
