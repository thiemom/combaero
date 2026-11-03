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
// T = 1e11 K as the static state. Since #937 it is refused up front, as an
// invalid argument.
TEST(InnerSolveCap868, ChamberStaticStateThrowsOnAProbeItCannotSolve) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(chamber_static_state(-1.393e19, 7.734e4, 85.02, X), std::invalid_argument);
  const ChamberStaticState cs = chamber_static_state(400.0, 1.0e5, 150.0, X);
  const double R = specific_gas_constant(X);
  const double v = 150.0 * R * cs.T / 1.0e5;
  EXPECT_NEAR(h_mass(400.0, X) - h_mass(cs.T, X), 0.5 * v * v, 1e-6);
}

// A target cp the clamped range [200, 6000] K never reaches returned 6000 K.
// Since R6 the scan of that range finds no sign change and says so; invalid
// input is std::invalid_argument (#928).
TEST(InnerSolveCap868, CalcTFromCpThrowsOnAnUnreachableTarget) {
  const std::vector<double> X = dry_air();
  EXPECT_THROW(calc_T_from_cp(1.0e3, X, 300.0), std::runtime_error);
  EXPECT_THROW(calc_T_from_cp(20.0, X, 300.0), std::runtime_error);
  EXPECT_THROW(calc_T_from_cp(-1.0, X, 300.0), std::invalid_argument);
  EXPECT_THROW(calc_T_from_cp(30.0, X, -1.0), std::invalid_argument);
  EXPECT_NEAR(calc_T_from_cp(cp(800.0, X), X, 300.0), 800.0, 1e-9);
}

// cp of N2 has a minimum near 214 K, so cp(205 K) has a second root near
// 226.7 K: the one returned is the root nearest T_guess (#868 R6).
TEST(InnerSolveCap868, CalcTFromCpTakesTheRootNearestTheGuess) {
  std::vector<double> X(num_species(), 0.0);
  X[species_index_from_name("N2")] = 1.0;
  const double target = cp(205.0, X);
  EXPECT_NEAR(calc_T_from_cp(target, X, 213.0), 205.0, 1e-9);
  const double upper = calc_T_from_cp(target, X, 218.0);
  EXPECT_GT(upper, 226.0);
  EXPECT_LT(upper, 227.5);
  EXPECT_NEAR(cp(upper, X) / target, 1.0, 1e-13);
}

// T_from_stagnation's slope dropped T dgamma/dT, so Newton converged
// linearly and stopped with |F| up to 8e-10 cp T at 2500 K (#868 R3).
TEST(InnerSolveCap868, TFromStagnationMeetsTheEnergyEquation) {
  const std::vector<double> X = dry_air();
  const double Tt = 2500.0;
  for (double M : {1.7782794100389228, 3.1622776601683795, 5.623413251903491}) {
    const double T = T_from_stagnation(Tt, M, X);
    const double v = M * speed_of_sound(T, X);
    const double F = h_mass(T, X) + 0.5 * v * v - h_mass(Tt, X);
    EXPECT_LT(std::abs(F) / (cp_mass(T, X) * T), 1e-13) << "M = " << M;
  }
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

// The static-T solves go to their root below the former 10 K floor (#935):
// the gas is ideal to 0 K, so a deep expansion has a positive static T
// however small. On main each of these returned 10 K.
TEST(InnerSolveCap868, StaticTemperatureGoesToItsRootBelowTheFormerFloor) {
  const std::vector<double> X = dry_air();
  // Hypersonic stagnation inverse: T ~ Tt / (0.2 M^2) = 0.6 K.
  const double M = 50.0;
  const double T = T_from_stagnation(300.0, M, X);
  EXPECT_LT(T, 1.0);
  EXPECT_GT(T, 0.0);
  const double a = speed_of_sound(T, X);
  const double res = h_mass(T, X) + 0.5 * M * M * a * a - h_mass(300.0, X);
  EXPECT_LT(std::abs(res), 1e-6 * cp_mass(T, X) * T);
  // Static to stagnation from a 5 K static state.
  const double T0 = Tt_from_static(5.0, 0.5, X);
  EXPECT_GT(T0, 5.0);
  EXPECT_LT(T0, thermo::T_FLOOR);
  // The isentrope to P/Pt = 1e-8 (static T ~ 1.5 K), against the independent
  // safeguarded inverse of s.
  const double Tt = 300.0, Pt = 1.0e5, P = 1.0e-3;
  const double s0 = s_mass(Tt, X, Pt);
  const double T_ref = calc_T_from_s_mass(s0, P, X, 1.0, 1e-13, 200);
  EXPECT_LT(T_ref, 2.0);
  const double v_ref = std::sqrt(2.0 * (h_mass(Tt, X) - h_mass(T_ref, X)));
  const double G_ref = density(T_ref, P, X) * v_ref;
  EXPECT_NEAR(mass_flux_isentropic(Tt, Pt, P, X) / G_ref, 1.0, 1e-10);
  EXPECT_NEAR(mach_from_pressure_ratio(Tt, Pt, P, X) / (v_ref / speed_of_sound(T_ref, X)), 1.0,
              1e-10);
  // Steps in ln T reach a root many decades down within the budget.
  const double G_deep = mass_flux_isentropic(Tt, Pt, 1.0e-295, X);
  EXPECT_TRUE(std::isfinite(G_deep));
  EXPECT_GT(G_deep, 0.0);
}

// s is exact (and finite) at every T > 0 down to S_R_T_MIN: the low tail
// keeps its fade branch where t rounds to 1 (it returned NaN from 1e-14 K);
// a non-positive probe T is held, finite.
TEST(InnerSolveCap868, EntropyIsExactBelowTheFormerFloor) {
  const std::vector<double> X = dry_air();
  const double s5 = s_mass(5.0, X, 1.0e5);
  const double s_20 = s_mass(1.0e-20, X, 1.0e5);
  const double s_19 = s_mass(1.0e-19, X, 1.0e5);
  EXPECT_TRUE(std::isfinite(s_20));
  EXPECT_LT(s_20, s_19);
  EXPECT_LT(s_19, s5);
  EXPECT_LT(s5, s_mass(thermo::T_FLOOR, X, 1.0e5));
  // One decade in T at cp = 7/2 R (rigid diatomic): ds = cp ln 10.
  const double cp_low = cp_mass(1.0e-20, X);
  EXPECT_NEAR((s_19 - s_20) / (cp_low * std::log(10.0)), 1.0, 1e-9);
  EXPECT_TRUE(std::isfinite(s_mass(0.0, X, 1.0e5)));
  EXPECT_TRUE(std::isfinite(s_mass(-5.0, X, 1.0e5)));
}
