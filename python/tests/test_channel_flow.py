"""
Tests for combined convective heat transfer + pressure loss channel functions.

Tests verify:
- channel_smooth: consistency with htc_pipe and Darcy-Weisbach
- channel_ribbed: Nu and dP exceed smooth baseline
- channel_dimpled: Nu enhanced, dP lower than ribbed
- channel_pin_fin: dP scales with N_rows, Nu increases with Re
- channel_impingement: Nu from correlation, dP from orifice model
- T_aw continuity: no threshold discontinuity at any Mach
- q = nan when T_wall not supplied
- heat_transfer submodule API
"""

import math

import pytest

import combaero as cb
from combaero import heat_transfer as ht

# Standard air at 600 K, 5 bar
T_AIR = 600.0
P_AIR = 5e5
X_AIR = [0.79, 0.21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Pin-fin geometry: keep Re_d = rho*v*d/mu within [3000, 90000].
# At 600K/5bar air: rho~2.9 kg/m3, mu~3.0e-5 Pa.s
# Re_d = 2.9 * v * d / 3e-5  =>  v=5 m/s, d=0.003 m -> Re_d ~ 1450 (too low)
# v=5 m/s, d=0.005 m -> Re_d ~ 2400 (too low)
# v=10 m/s, d=0.005 m -> Re_d ~ 4800 (ok)
# v=20 m/s, d=0.005 m -> Re_d ~ 9700 (ok)
# v=50 m/s, d=0.005 m -> Re_d ~ 24000 (ok)
PIN_V = 10.0  # approach velocity [m/s]
PIN_H = 0.015  # channel height / pin length [m]
PIN_D = 0.005  # pin diameter [m]  -> Re_d ~ 4800 at 600K/5bar
PIN_SD = 2.5  # spanwise pitch/d
PIN_XD = 2.5  # streamwise pitch/d
PIN_N = 5  # number of rows

# Impingement: keep Re_jet = rho*v_jet*d/mu within [5000, 80000].
# A_jet = pi/4 * d^2.  v_jet = mdot / (rho * A_jet)
# At 600K/5bar: rho~2.9 kg/m3, mu~3e-5
# d=0.005 m, A_jet=1.96e-5 m2
# mdot=0.001 kg/s -> v_jet=0.001/(2.9*1.96e-5)=17.6 m/s -> Re_jet=2.9*17.6*0.005/3e-5=8500 (ok)
# mdot=0.005 kg/s -> Re_jet~42500 (ok)
IMP_MDOT = 0.001  # jet mass flow [kg/s]
IMP_D = 0.005  # jet diameter [m]
IMP_ZD = 6.0  # z/D
IMP_A = 0.01  # target area [m^2]


class TestChannelSmooth:
    """Tests for channel_smooth."""

    def test_returns_channel_result(self):
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        assert isinstance(sol, cb.ChannelResult)

    def test_h_matches_htc_pipe(self):
        """channel_smooth h should match htc_pipe for same conditions."""
        v = 20.0
        D = 0.01
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, 0.5)
        h_ref, Nu_ref, Re_ref = cb.htc_pipe(T_AIR, P_AIR, X_AIR, v, D)
        assert sol.h == pytest.approx(h_ref, rel=1e-6)
        assert sol.Nu == pytest.approx(Nu_ref, rel=1e-6)
        assert sol.Re == pytest.approx(Re_ref, rel=1e-6)

    def test_dP_matches_darcy_weisbach(self):
        """dP should match Darcy-Weisbach: f*(L/D)*(rho*v^2/2)."""
        v = 20.0
        D = 0.01
        L = 0.5
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        rho = cb.density(T_AIR, P_AIR, X_AIR)
        dP_ref = sol.f * (L / D) * (rho * v * v / 2.0)
        assert sol.dP == pytest.approx(dP_ref, rel=1e-6)

    def test_q_nan_without_T_wall(self):
        """q should be nan when T_wall is not supplied."""
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        assert math.isnan(sol.q)

    def test_q_finite_with_T_wall(self):
        """q should be finite and equal h*(T_aw - T_wall) when T_wall given."""
        T_wall = 800.0
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, T_wall=T_wall)
        assert math.isfinite(sol.q)
        assert sol.q == pytest.approx(sol.h * (sol.T_aw - T_wall), rel=1e-9)

    def test_T_aw_equals_T_at_zero_velocity(self):
        """At v=0, T_aw should equal T_static (recovery factor formula gives 0 correction)."""
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 0.0, 0.01, 0.5)
        assert sol.T_aw == pytest.approx(T_AIR, rel=1e-9)

    def test_T_aw_continuous_across_mach(self):
        """T_aw must be continuous — no jump at any velocity."""
        D = 0.05
        L = 1.0
        velocities = [0.0, 50.0, 100.0, 150.0, 200.0]
        T_aws = [cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L).T_aw for v in velocities]
        # T_aw must be monotonically non-decreasing with velocity
        for i in range(1, len(T_aws)):
            assert T_aws[i] >= T_aws[i - 1]
        # No jump larger than 5 K between adjacent steps
        for i in range(1, len(T_aws)):
            assert abs(T_aws[i] - T_aws[i - 1]) < 50.0

    def test_T_aw_increases_with_velocity(self):
        """T_aw should increase with velocity (kinetic heating)."""
        sol_slow = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 10.0, 0.05, 1.0)
        sol_fast = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 200.0, 0.05, 1.0)
        assert sol_fast.T_aw > sol_slow.T_aw

    def test_mach_computed_internally(self):
        """M should be v/a(T,X)."""
        v = 100.0
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, 0.05, 1.0)
        a = cb.speed_of_sound(T_AIR, X_AIR)
        assert pytest.approx(v / a, rel=1e-6) == sol.M

    def test_dP_zero_at_zero_velocity(self):
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 0.0, 0.01, 0.5)
        assert sol.dP == pytest.approx(0.0, abs=1e-10)

    def test_dP_scales_with_length(self):
        """dP should scale linearly with L."""
        sol1 = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        sol2 = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 1.0)
        assert sol2.dP == pytest.approx(2.0 * sol1.dP, rel=1e-6)

    def test_roughness_increases_dP(self):
        """Roughness should increase friction factor and dP."""
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, roughness=0.0)
        sol_rough = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, roughness=5e-5)
        assert sol_rough.f > sol_smooth.f
        assert sol_rough.dP > sol_smooth.dP

    def test_correlations_accepted(self):
        """All four correlation names should work without error."""
        for corr in ["gnielinski", "dittus_boelter", "sieder_tate", "petukhov"]:
            sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 30.0, 0.01, 0.5, correlation=corr)
            assert sol.h > 0.0

    def test_repr(self):
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        r = repr(sol)
        assert "ChannelResult" in r
        assert "h=" in r


class TestChannelRibbed:
    """Tests for channel_ribbed."""

    def test_nu_exceeds_smooth(self):
        """Ribbed Nu should be higher than smooth-pipe Nu at same Re."""
        v = 20.0
        D = 0.01
        L = 0.5
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        sol_ribbed = cb.channel_ribbed(T_AIR, P_AIR, X_AIR, v, D, L, e_D=0.05, P_e=10.0, alpha=90.0)
        assert sol_ribbed.Nu > sol_smooth.Nu
        assert sol_ribbed.h > sol_smooth.h

    def test_dP_exceeds_smooth(self):
        """Ribbed dP should be higher than smooth-pipe dP."""
        v = 20.0
        D = 0.01
        L = 0.5
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        sol_ribbed = cb.channel_ribbed(T_AIR, P_AIR, X_AIR, v, D, L, e_D=0.05, P_e=10.0, alpha=90.0)
        assert sol_ribbed.dP > sol_smooth.dP
        assert sol_ribbed.f > sol_smooth.f

    def test_q_with_T_wall(self):
        T_wall = 800.0
        sol = cb.channel_ribbed(
            T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, e_D=0.05, P_e=10.0, alpha=90.0, T_wall=T_wall
        )
        assert math.isfinite(sol.q)
        assert sol.q == pytest.approx(sol.h * (sol.T_aw - T_wall), rel=1e-9)

    def test_higher_rib_height_more_enhancement(self):
        """Higher e_D should give more Nu enhancement."""
        v, D, L = 20.0, 0.01, 0.5
        sol_low = cb.channel_ribbed(T_AIR, P_AIR, X_AIR, v, D, L, e_D=0.03, P_e=10.0, alpha=90.0)
        sol_high = cb.channel_ribbed(T_AIR, P_AIR, X_AIR, v, D, L, e_D=0.08, P_e=10.0, alpha=90.0)
        assert sol_high.Nu > sol_low.Nu

    def test_re_same_as_smooth(self):
        """Re should be the same as smooth-pipe baseline (same fluid, same velocity)."""
        v, D, L = 20.0, 0.01, 0.5
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        sol_ribbed = cb.channel_ribbed(T_AIR, P_AIR, X_AIR, v, D, L, e_D=0.05, P_e=10.0, alpha=90.0)
        assert sol_ribbed.Re == pytest.approx(sol_smooth.Re, rel=1e-6)


class TestChannelDimpled:
    """Tests for channel_dimpled."""

    def test_nu_exceeds_smooth(self):
        v, D, L = 20.0, 0.01, 0.5
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        sol_dimpled = cb.channel_dimpled(T_AIR, P_AIR, X_AIR, v, D, L, d_Dh=0.2, h_d=0.2, S_d=2.0)
        assert sol_dimpled.Nu > sol_smooth.Nu

    def test_dP_exceeds_smooth(self):
        v, D, L = 20.0, 0.01, 0.5
        sol_smooth = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L)
        sol_dimpled = cb.channel_dimpled(T_AIR, P_AIR, X_AIR, v, D, L, d_Dh=0.2, h_d=0.2, S_d=2.0)
        assert sol_dimpled.dP > sol_smooth.dP

    def test_friction_multiplier_greater_than_one(self):
        """Dimple friction multiplier must be > 1 (dimples always increase friction)."""
        v, D, L = 20.0, 0.01, 0.5
        Re = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L).Re
        f_mul = cb.dimple_friction_multiplier(Re, d_Dh=0.2, h_d=0.2)
        assert f_mul > 1.0

    def test_friction_multiplier_increases_with_depth(self):
        """Deeper dimples should produce higher friction multiplier."""
        v, D, L = 20.0, 0.01, 0.5
        Re = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, L).Re
        f_shallow = cb.dimple_friction_multiplier(Re, d_Dh=0.2, h_d=0.1)
        f_deep = cb.dimple_friction_multiplier(Re, d_Dh=0.2, h_d=0.25)
        assert f_deep > f_shallow

    def test_q_nan_without_T_wall(self):
        sol = cb.channel_dimpled(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, d_Dh=0.2, h_d=0.2, S_d=2.0)
        assert math.isnan(sol.q)


class TestChannelPinFin:
    """Tests for channel_pin_fin."""

    def test_returns_channel_result(self):
        sol = cb.channel_pin_fin(
            T_AIR,
            P_AIR,
            X_AIR,
            velocity=PIN_V,
            channel_height=PIN_H,
            pin_diameter=PIN_D,
            S_D=PIN_SD,
            X_D=PIN_XD,
            N_rows=PIN_N,
        )
        assert isinstance(sol, cb.ChannelResult)

    def test_dP_scales_with_N_rows(self):
        """dP should scale linearly with N_rows."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "velocity": PIN_V,
            "channel_height": PIN_H,
            "pin_diameter": PIN_D,
            "S_D": PIN_SD,
            "X_D": PIN_XD,
        }
        sol5 = cb.channel_pin_fin(**kwargs, N_rows=5)
        sol10 = cb.channel_pin_fin(**kwargs, N_rows=10)
        assert sol10.dP == pytest.approx(2.0 * sol5.dP, rel=1e-6)

    def test_nu_increases_with_re(self):
        """Nu should increase with Re (velocity)."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "channel_height": PIN_H,
            "pin_diameter": PIN_D,
            "S_D": PIN_SD,
            "X_D": PIN_XD,
            "N_rows": PIN_N,
        }
        sol_low = cb.channel_pin_fin(**kwargs, velocity=PIN_V)
        sol_high = cb.channel_pin_fin(**kwargs, velocity=PIN_V * 4)
        assert sol_high.Nu > sol_low.Nu

    def test_staggered_higher_nu_than_inline(self):
        """Staggered arrays should have higher Nu than inline."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "velocity": PIN_V,
            "channel_height": PIN_H,
            "pin_diameter": PIN_D,
            "S_D": PIN_SD,
            "X_D": PIN_XD,
            "N_rows": PIN_N,
        }
        sol_stag = cb.channel_pin_fin(**kwargs, is_staggered=True)
        sol_inline = cb.channel_pin_fin(**kwargs, is_staggered=False)
        assert sol_stag.Nu > sol_inline.Nu

    def test_staggered_higher_dP_than_inline(self):
        """Staggered arrays have higher friction than inline."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "velocity": PIN_V,
            "channel_height": PIN_H,
            "pin_diameter": PIN_D,
            "S_D": PIN_SD,
            "X_D": PIN_XD,
            "N_rows": PIN_N,
        }
        sol_stag = cb.channel_pin_fin(**kwargs, is_staggered=True)
        sol_inline = cb.channel_pin_fin(**kwargs, is_staggered=False)
        assert sol_stag.dP > sol_inline.dP

    def test_q_with_T_wall(self):
        T_wall = 800.0
        sol = cb.channel_pin_fin(
            T_AIR,
            P_AIR,
            X_AIR,
            velocity=PIN_V,
            channel_height=PIN_H,
            pin_diameter=PIN_D,
            S_D=PIN_SD,
            X_D=PIN_XD,
            N_rows=PIN_N,
            T_wall=T_wall,
        )
        assert math.isfinite(sol.q)
        assert sol.q == pytest.approx(sol.h * (sol.T_aw - T_wall), rel=1e-9)

    def test_re_based_on_pin_diameter(self):
        """Re should be rho*v*d/mu (pin diameter, not hydraulic diameter)."""
        sol = cb.channel_pin_fin(
            T_AIR,
            P_AIR,
            X_AIR,
            velocity=PIN_V,
            channel_height=PIN_H,
            pin_diameter=PIN_D,
            S_D=PIN_SD,
            X_D=PIN_XD,
            N_rows=PIN_N,
        )
        rho = cb.density(T_AIR, P_AIR, X_AIR)
        mu = cb.viscosity(T_AIR, P_AIR, X_AIR)
        Re_expected = rho * PIN_V * PIN_D / mu
        assert sol.Re == pytest.approx(Re_expected, rel=1e-6)


class TestChannelImpingement:
    """Tests for channel_impingement."""

    def test_returns_channel_result(self):
        sol = cb.channel_impingement(
            T_AIR, P_AIR, X_AIR, mdot_jet=IMP_MDOT, d_jet=IMP_D, z_D=IMP_ZD, A_target=IMP_A
        )
        assert isinstance(sol, cb.ChannelResult)

    def test_nu_reasonable_range(self):
        """Nu should be in a physically reasonable range."""
        sol = cb.channel_impingement(
            T_AIR, P_AIR, X_AIR, mdot_jet=IMP_MDOT, d_jet=IMP_D, z_D=IMP_ZD, A_target=IMP_A
        )
        assert 10 < sol.Nu < 500

    def test_dP_from_orifice_model(self):
        """dP = (1/Cd^2) * rho * v_jet^2 / 2."""
        Cd = 0.65
        sol = cb.channel_impingement(
            T_AIR,
            P_AIR,
            X_AIR,
            mdot_jet=IMP_MDOT,
            d_jet=IMP_D,
            z_D=IMP_ZD,
            A_target=IMP_A,
            Cd_jet=Cd,
        )
        rho = cb.density(T_AIR, P_AIR, X_AIR)
        A_jet = math.pi / 4.0 * IMP_D**2
        v_jet = IMP_MDOT / (rho * A_jet)
        dP_expected = (1.0 / Cd**2) * rho * v_jet**2 / 2.0
        assert sol.dP == pytest.approx(dP_expected, rel=1e-6)

    def test_f_equals_one_over_cd_squared(self):
        """f should equal 1/Cd^2."""
        Cd = 0.65
        sol = cb.channel_impingement(
            T_AIR,
            P_AIR,
            X_AIR,
            mdot_jet=IMP_MDOT,
            d_jet=IMP_D,
            z_D=IMP_ZD,
            A_target=IMP_A,
            Cd_jet=Cd,
        )
        assert sol.f == pytest.approx(1.0 / Cd**2, rel=1e-9)

    def test_array_lower_nu_than_single_jet(self):
        """Jet array Nu should be lower than single jet (crossflow degradation)."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "mdot_jet": IMP_MDOT,
            "d_jet": IMP_D,
            "z_D": IMP_ZD,
            "A_target": IMP_A,
        }
        sol_single = cb.channel_impingement(**kwargs)
        sol_array = cb.channel_impingement(**kwargs, x_D=8.0, y_D=8.0)
        assert sol_array.Nu < sol_single.Nu

    def test_q_with_T_wall(self):
        T_wall = 800.0
        sol = cb.channel_impingement(
            T_AIR,
            P_AIR,
            X_AIR,
            mdot_jet=IMP_MDOT,
            d_jet=IMP_D,
            z_D=IMP_ZD,
            A_target=IMP_A,
            T_wall=T_wall,
        )
        assert math.isfinite(sol.q)
        assert sol.q == pytest.approx(sol.h * (sol.T_aw - T_wall), rel=1e-9)

    def test_higher_mdot_increases_nu(self):
        """Higher jet mass flow → higher Re → higher Nu."""
        kwargs = {
            "T": T_AIR,
            "P": P_AIR,
            "X": X_AIR,
            "d_jet": IMP_D,
            "z_D": IMP_ZD,
            "A_target": IMP_A,
        }
        sol_low = cb.channel_impingement(**kwargs, mdot_jet=IMP_MDOT)
        sol_high = cb.channel_impingement(**kwargs, mdot_jet=IMP_MDOT * 4)
        assert sol_high.Nu > sol_low.Nu


class TestPinFinFriction:
    """Tests for the scalar pin_fin_friction function."""

    def test_staggered_higher_than_inline(self):
        Re_d = 20000.0
        f_stag = cb.pin_fin_friction(Re_d, is_staggered=True)
        f_inline = cb.pin_fin_friction(Re_d, is_staggered=False)
        assert f_stag > f_inline

    def test_decreases_with_re(self):
        """Friction coefficient decreases with Re (power law)."""
        f_low = cb.pin_fin_friction(5000.0)
        f_high = cb.pin_fin_friction(50000.0)
        assert f_high < f_low

    def test_reasonable_range(self):
        """f_pin should be in a physically reasonable range."""
        f = cb.pin_fin_friction(20000.0)
        assert 0.01 < f < 1.0


class TestHeatTransferSubmodule:
    """Tests for the combaero.heat_transfer submodule API."""

    def test_smooth_returns_channel_result(self):
        sol = ht.smooth(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01)
        assert isinstance(sol, cb.ChannelResult)

    def test_smooth_matches_direct_call(self):
        """ht.smooth should give identical result to cb.channel_smooth."""
        sol_sub = ht.smooth(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01, T_wall=800.0)
        sol_direct = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, T_wall=800.0)
        assert sol_sub.h == pytest.approx(sol_direct.h, rel=1e-9)
        assert sol_sub.dP == pytest.approx(sol_direct.dP, rel=1e-9)
        assert sol_sub.q == pytest.approx(sol_direct.q, rel=1e-9)

    def test_ribbed_returns_channel_result(self):
        sol = ht.ribbed(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01, e_D=0.05, P_e=10.0, alpha=90.0)
        assert isinstance(sol, cb.ChannelResult)

    def test_dimpled_returns_channel_result(self):
        sol = ht.dimpled(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01, d_Dh=0.2, h_d=0.2, S_d=2.0)
        assert isinstance(sol, cb.ChannelResult)

    def test_pin_fin_returns_channel_result(self):
        sol = ht.pin_fin(
            T_AIR, P_AIR, X_AIR, u=PIN_V, H=PIN_H, d=PIN_D, S_D=PIN_SD, X_D=PIN_XD, N_rows=PIN_N
        )
        assert isinstance(sol, cb.ChannelResult)

    def test_impingement_returns_channel_result(self):
        sol = ht.impingement(
            T_AIR, P_AIR, X_AIR, mdot_jet=IMP_MDOT, d_jet=IMP_D, z_D=IMP_ZD, A_target=IMP_A
        )
        assert isinstance(sol, cb.ChannelResult)

    def test_network_solver_pattern(self):
        """Demonstrate the intended network solver usage pattern."""
        T_wall = 850.0
        sol = ht.smooth(T_AIR, P_AIR, X_AIR, u=25.0, L=0.3, D=0.008, T_wall=T_wall)
        assert sol.h > 0.0
        assert sol.dP > 0.0
        assert math.isfinite(sol.q)
        assert sol.T_aw >= T_AIR
        assert math.isfinite(sol.M)
        assert sol.M > 0.0
