#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "compressible.h"
#include "thermo.h"

// critical_pressure_ratio is the PEAK of the isentropic mass flux, where the
// flow is sonic. Until #471 its golden-section search kept its probe points
// in the opposite order to its update and shrank onto a fixed point near
// 0.5636 whatever the gas or temperature: 6% high for air, choked mass flux
// 0.1-0.3% low, and choking declared early by nozzle_flow.

namespace {
std::vector<double> air() {
  std::vector<double> X(combaero::num_species(), 0.0);
  X[combaero::species_index_from_name("N2")] = 0.7808;
  X[combaero::species_index_from_name("O2")] = 0.2095;
  X[combaero::species_index_from_name("AR")] = 0.0097;
  return X;
}
}  // namespace

TEST(CriticalPressureRatio, IsWhereTheMassFluxPeaksAndTheFlowIsSonic) {
  const auto X = air();
  const double P0 = 3.0e5;
  for (double T0 : {300.0, 600.0, 1500.0}) {
    const double r = combaero::critical_pressure_ratio(T0, P0, X);
    const double G = combaero::mass_flux_isentropic(T0, P0, r * P0, X);
    EXPECT_GE(G, combaero::mass_flux_isentropic(T0, P0, (r + 0.002) * P0, X)) << T0;
    EXPECT_GE(G, combaero::mass_flux_isentropic(T0, P0, (r - 0.002) * P0, X)) << T0;
    EXPECT_NEAR(combaero::mach_from_pressure_ratio(T0, P0, r * P0, X), 1.0, 2e-3) << T0;
  }
}

TEST(CriticalPressureRatio, MovesWithTemperatureAsGammaDoes) {
  // gamma falls as air heats, so P*/P0 = (2/(gamma+1))^(gamma/(gamma-1)) rises.
  const auto X = air();
  const double cold = combaero::critical_pressure_ratio(300.0, 3.0e5, X);
  const double hot = combaero::critical_pressure_ratio(1500.0, 3.0e5, X);
  EXPECT_NEAR(cold, 0.528, 2e-3);
  EXPECT_GT(hot, cold + 0.005);
}
