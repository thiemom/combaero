#include <gtest/gtest.h>

#include "../include/compressible.h"
#include "../include/fanno_mach.h"
#include "../include/friction.h"
#include "../include/humidair.h"
#include "../include/inner_solve.h"
#include "../include/stagnation.h"
#include "../include/thermo.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace combaero;

// ---------------------------------------------------------------------------
// Inner solves report a cap exit instead of returning a non-root (#868).
//
// Each case below reaches a condition the audit found (or one a caller's
// arguments reach): on main they returned a value that is not a root. The
// converged cases next to them pin that the normal path is unchanged.
// See docs/archive/AUDIT_868_INNER_SOLVES.md.
// ---------------------------------------------------------------------------

TEST(InnerSolveCap868, CapExitIsRootRejectsNaNAndWideSteps) {
  EXPECT_TRUE(inner_solve::cap_exit_is_root(1e-12, 1.0));
  EXPECT_TRUE(inner_solve::cap_exit_is_root(-1e-12, -1.0));
  EXPECT_FALSE(inner_solve::cap_exit_is_root(1e-6, 1.0));
  EXPECT_FALSE(inner_solve::cap_exit_is_root(std::numeric_limits<double>::quiet_NaN(), 1.0));
}

// A solver probe with a garbage stagnation T (-1.4e19 K, from the cascaded tee
// manifold test) ran chamber_static_state's Newton to its cap and returned
// T = 1e11 K as the static state.
TEST(InnerSolveCap868, ChamberStaticStateThrowsOnAProbeItCannotSolve) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(chamber_static_state(-1.393e19, 7.734e4, 85.02, X), std::runtime_error);
  const ChamberStaticState cs = chamber_static_state(400.0, 1.0e5, 150.0, X);
  const double R = specific_gas_constant(X);
  const double v = 150.0 * R * cs.T / 1.0e5;
  EXPECT_NEAR(h_mass(400.0, X) - h_mass(cs.T, X), 0.5 * v * v, 1e-6);
}

// A target cp the clamped range [200, 6000] K never reaches returned 6000 K.
TEST(InnerSolveCap868, CalcTFromCpThrowsOnAnUnreachableTarget) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(calc_T_from_cp(1.0e3, X, 300.0), std::runtime_error);
  EXPECT_NEAR(calc_T_from_cp(cp(800.0, X), X, 300.0), 800.0, 1e-6);
}

// Below G ~ 1e-10 of the sonic flux the central-difference stencil crossed
// M = 0 (where the flux is held at 0), read half the slope and Newton
// 2-cycled about the root: M_in came out 36-73 % high. The inlet Mach is
// linear in G there, so the ratio below is 1.
TEST(InnerSolveCap868, FannoInletMachIsLinearAtVanishingFlux) {
  const std::vector<double> X = dry_air();
  const double G_max = fanno_sonic_mass_flux(2.0e5, 300.0, X);
  const double M_ref = fanno_inlet_mach(2.0e5, 300.0, 1e-8 * G_max, X);
  for (double r : {1e-12, 1e-20, 1e-25}) {
    const double M = fanno_inlet_mach(2.0e5, 300.0, r * G_max, X);
    EXPECT_NEAR(M / (M_ref * r / 1e-8), 1.0, 1e-9) << "G/G_max = " << r;
  }
}

// An unphysical exit stagnation T (2.1e5 K, a solver probe in the C3X vane
// example) left the choke search with no sign change down to e^-61 of the
// sonic flux; Brent then returned a point of that interval as the choked
// flux (G = 6.3e-5, choked = true). The same duct at a physical exit T solves.
TEST(InnerSolveCap868, RayleighChokeSearchWithoutASignChangeThrows) {
  const std::vector<double> X = dry_air();
  auto flow = [&](double Tt_exit) {
    return fanno_rayleigh_channel_flow(266266.38, 1603.23, Tt_exit, X, 136765.9, false, 0.0142,
                                       0.0897, 0.0, "haaland", 1.0, 0.0636, {}, false);
  };
  EXPECT_THROW(flow(211975.1), std::runtime_error);
  const FannoChannelFlow ok = flow(2000.0);
  EXPECT_GT(ok.G, 0.0);
  EXPECT_FALSE(ok.choked);
}

// A caller's iteration budget is honoured: a cap exit throws rather than
// returning the half-converged iterate.
TEST(InnerSolveCap868, CallerBudgetsAreReported) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(friction_colebrook(1.0e5, 1.0e-4, 1e-10, 1), std::runtime_error);
  EXPECT_GT(friction_colebrook(1.0e5, 1.0e-4), 0.0);
  EXPECT_THROW(Tt_from_static(300.0, 2.0, X, 1e-8, 1), std::runtime_error);
  EXPECT_THROW(T_from_stagnation(900.0, 2.0, X, 1e-8, 1), std::runtime_error);
}

// The documented T_FLOOR clamp is unchanged: an expansion whose isentropic
// static T lies below the floor returns the floor (now without running on to
// the cap), it does not throw.
TEST(InnerSolveCap868, TemperatureFloorIsStillTheAnswerBelowIt) {
  const std::vector<double> X = dry_air();
  EXPECT_EQ(T_from_stagnation(300.0, 50.0, X), thermo::T_FLOOR);
  EXPECT_TRUE(std::isfinite(mach_from_pressure_ratio(300.0, 1.0e5, 1.0e-3, X)));
  EXPECT_GT(mass_flux_isentropic(300.0, 1.0e5, 1.0e-3, X), 0.0);
}
