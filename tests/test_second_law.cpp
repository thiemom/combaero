// Second law of a control volume (#545): entropy_generation and
// heat_entropy_flux.
#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "composition.h"
#include "solver_interface.h"
#include "state.h"
#include "thermo.h"

using namespace combaero;
namespace sv = combaero::solver;

namespace {

std::vector<double> air_Y() {
  std::vector<double> Y(num_species(), 0.0);
  for (std::size_t k = 0; k < num_species(); ++k) {
    if (species_name(k) == "N2") {
      Y[k] = 0.767;
    } else if (species_name(k) == "O2") {
      Y[k] = 0.233;
    }
  }
  return Y;
}

std::vector<double> ch4_Y() {
  std::vector<double> Y(num_species(), 0.0);
  for (std::size_t k = 0; k < num_species(); ++k) {
    if (species_name(k) == "CH4") {
      Y[k] = 1.0;
    }
  }
  return Y;
}

}  // namespace

// An adiabatic two-port of fixed composition: sigma_norm = -ln(Pt2 / Pt1),
// compressible flow included, at any temperature.
TEST(SecondLaw, AdiabaticTwoPortIsTheTotalPressureRatio) {
  const auto Y = air_Y();
  for (double T : {300.0, 900.0, 1800.0}) {
    for (double ratio : {0.5, 0.95, 1.0, 1.02}) {
      const auto r = sv::entropy_generation({{0.3, T, 2.0e5, Y}}, {{0.3, T, 2.0e5 * ratio, Y}});
      EXPECT_NEAR(r.sigma_norm, -std::log(ratio), 1e-12) << "T " << T << " ratio " << ratio;
    }
  }
}

// Mixing two streams at the same total pressure generates entropy: by
// temperature, and by composition (the ideal mixing term).
TEST(SecondLaw, MixingAtEqualPressureGeneratesEntropy) {
  const auto air = air_Y();
  const double P = 1.5e5;
  // hot + cold air, adiabatic mixing at constant cp-weighted enthalpy
  const double h_mix = 0.5 * (h_mass(300.0, mass_to_mole(air)) + h_mass(900.0, mass_to_mole(air)));
  const double T_mix = calc_T_from_h_mass(h_mix, mass_to_mole(air), 600.0);
  const auto thermal = sv::entropy_generation({{0.5, 300.0, P, air}, {0.5, 900.0, P, air}},
                                              {{1.0, T_mix, P, air}});
  EXPECT_GT(thermal.sigma, 0.0);

  // air + methane at one temperature
  const auto ch4 = ch4_Y();
  std::vector<double> Y_mix(air.size());
  for (std::size_t k = 0; k < air.size(); ++k) {
    Y_mix[k] = 0.9 * air[k] + 0.1 * ch4[k];
  }
  const double H = 0.9 * h_mass(400.0, mass_to_mole(air)) + 0.1 * h_mass(400.0, mass_to_mole(ch4));
  const double T_m = calc_T_from_h_mass(H, mass_to_mole(Y_mix), 400.0);
  const auto species = sv::entropy_generation({{0.9, 400.0, P, air}, {0.1, 400.0, P, ch4}},
                                              {{1.0, T_m, P, Y_mix}});
  EXPECT_GT(species.sigma, 0.0);
}

// Heat at M = 0 brings in exactly the entropy the stream gains at constant
// pressure, so a heated frictionless two-port generates none.
TEST(SecondLaw, HeatAtZeroMachIsTheConstantPressureEntropyRise) {
  const auto Y = air_Y();
  const auto X = mass_to_mole(Y);
  const double m = 0.02;
  for (auto [T1, T2] : {std::pair{667.0, 1351.0}, std::pair{900.0, 400.0}}) {
    const double S = sv::heat_entropy_flux(m, T1, T2, Y);
    EXPECT_NEAR(S, m * (s_mass(T2, X, 1.0e5) - s_mass(T1, X, 1.0e5)), 1e-9 * std::abs(S));
    const auto r = sv::entropy_generation({{m, T1, 2.0e5, Y}}, {{m, T2, 2.0e5, Y}}, S);
    EXPECT_NEAR(r.sigma_norm, 0.0, 1e-9);
  }
}

// At Mach M the static temperature is lower, so heating brings in MORE
// entropy than at M = 0 (Rayleigh: heating costs total pressure) and cooling
// takes out more (cooling may raise it).
TEST(SecondLaw, MachShiftsTheHeatEntropyTheRayleighWay) {
  const auto Y = air_Y();
  const double heat0 = sv::heat_entropy_flux(0.02, 600.0, 900.0, Y);
  const double heat3 = sv::heat_entropy_flux(0.02, 600.0, 900.0, Y, 0.3, 0.3);
  EXPECT_GT(heat3, heat0);
  const double cool0 = sv::heat_entropy_flux(0.02, 900.0, 600.0, Y);
  const double cool3 = sv::heat_entropy_flux(0.02, 900.0, 600.0, Y, 0.3, 0.3);
  EXPECT_LT(cool3, cool0);
  // linear Mach path lies between the two constant ones
  const double lin = sv::heat_entropy_flux(0.02, 600.0, 900.0, Y, 0.1, 0.3);
  EXPECT_GT(lin, sv::heat_entropy_flux(0.02, 600.0, 900.0, Y, 0.1, 0.1));
  EXPECT_LT(lin, heat3);
}

TEST(SecondLaw, RefusesWhatIsNotAControlVolume) {
  const auto Y = air_Y();
  EXPECT_THROW(sv::entropy_generation({}, {{0.1, 300.0, 1e5, Y}}), std::invalid_argument);
  EXPECT_THROW(sv::entropy_generation({{-0.1, 300.0, 1e5, Y}}, {}), std::invalid_argument);
  EXPECT_THROW(sv::heat_entropy_flux(0.1, 0.0, 300.0, Y), std::invalid_argument);
  EXPECT_THROW(sv::heat_entropy_flux(0.1, 300.0, 400.0, Y, -0.1, 0.0), std::invalid_argument);
}
