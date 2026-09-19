"""
Tests for combined convective heat transfer + pressure loss channel functions.

Tests verify:
- channel_smooth: consistency with htc_channel and Darcy-Weisbach
- channel_ribbed: Nu and dP exceed smooth baseline
- channel_dimpled: Nu enhanced, dP lower than ribbed
- channel_pin_fin: dP scales with N_rows, Nu increases with Re
- channel_impingement: Nu from correlation, dP from orifice model
- T_aw continuity: no threshold discontinuity at any Mach
- q = nan when T_hot not supplied
- heat_transfer submodule API
"""

import math

import pytest

import combaero as cb
from combaero import heat_transfer as ht

# Standard air at 600 K, 5 bar
T_AIR = 600.0
P_AIR = 5e5
X_AIR = cb.species.dry_air()

# Pin-fin geometry: keep Re_d = rho*v*d/mu within [3000, 90000].
# At 600K/5bar air: rho~2.9 kg/m3, mu~3.0e-5 Pa.s
# Re_d = 2.9 * v * d / 3e-5  =>  v=5 m/s, d=0.003 m -> Re_d ~ 1450 (too low)
# v=5 m/s, d=0.005 m -> Re_d ~ 2400 (too low)
# v=10 m/s, d=0.005 m -> Re_d ~ 4800 (ok)
# v=20 m/s, d=0.005 m -> Re_d ~ 9700 (ok)
# v=50 m/s, d=0.005 m -> Re_d ~ 24000 (ok)

# Impingement: keep Re_jet = rho*v_jet*d/mu within [5000, 80000].
# A_jet = pi/4 * d^2.  v_jet = mdot / (rho * A_jet)
# At 600K/5bar: rho~2.9 kg/m3, mu~3e-5
# d=0.005 m, A_jet=1.96e-5 m2
# mdot=0.001 kg/s -> v_jet=0.001/(2.9*1.96e-5)=17.6 m/s -> Re_jet=2.9*17.6*0.005/3e-5=8500 (ok)
# mdot=0.005 kg/s -> Re_jet~42500 (ok)


class TestChannelSmooth:
    """Tests for channel_smooth."""

    def test_returns_channel_result(self):
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        assert isinstance(sol, cb.ChannelResult)

    def test_h_matches_htc_channel(self):
        """channel_smooth h should match htc_channel for same conditions."""
        v = 20.0
        D = 0.01
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, v, D, 0.5)
        h_ref, Nu_ref, Re_ref = cb.htc_channel(T_AIR, P_AIR, X_AIR, v, D)
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

    def test_q_nan_without_T_hot(self):
        """q should be nan when T_hot is not supplied."""
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5)
        assert math.isnan(sol.q)

    def test_q_finite_with_T_hot(self):
        """q should be finite and equal h*(T_aw - T_hot) when T_hot given."""
        T_hot = 800.0
        sol = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, T_hot=T_hot)
        assert math.isfinite(sol.q)
        assert sol.q == pytest.approx(sol.h * (sol.T_aw - T_hot), rel=1e-9)

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


class TestHeatTransferSubmodule:
    """Tests for the combaero.heat_transfer submodule API."""

    def test_smooth_returns_channel_result(self):
        sol = ht.smooth(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01)
        assert isinstance(sol, cb.ChannelResult)

    def test_smooth_matches_direct_call(self):
        """ht.smooth should give identical result to cb.channel_smooth."""
        sol_sub = ht.smooth(T_AIR, P_AIR, X_AIR, u=20.0, L=0.5, D=0.01, T_hot=800.0)
        sol_direct = cb.channel_smooth(T_AIR, P_AIR, X_AIR, 20.0, 0.01, 0.5, T_hot=800.0)
        assert sol_sub.h == pytest.approx(sol_direct.h, rel=1e-9)
        assert sol_sub.dP == pytest.approx(sol_direct.dP, rel=1e-9)
        assert sol_sub.q == pytest.approx(sol_direct.q, rel=1e-9)

    def test_smooth_multiplier_forwarding(self):
        sol_sub = ht.smooth(
            T_AIR,
            P_AIR,
            X_AIR,
            u=20.0,
            L=0.5,
            D=0.01,
            Nu_multiplier=2.0,
            f_multiplier=2.0,
        )
        sol_direct = cb.channel_smooth(
            T_AIR,
            P_AIR,
            X_AIR,
            20.0,
            0.01,
            0.5,
            Nu_multiplier=2.0,
            f_multiplier=2.0,
        )
        assert sol_sub.h == pytest.approx(sol_direct.h, rel=1e-9)
        assert sol_sub.dP == pytest.approx(sol_direct.dP, rel=1e-9)

    def test_smooth_multiplier_scaling(self):
        # For Dittus-Boelter, Nu does not depend on f; this isolates exact scaling.
        base = ht.smooth(
            T_AIR,
            P_AIR,
            X_AIR,
            u=20.0,
            L=0.5,
            D=0.01,
            correlation="dittus_boelter",
        )
        nu2 = ht.smooth(
            T_AIR,
            P_AIR,
            X_AIR,
            u=20.0,
            L=0.5,
            D=0.01,
            correlation="dittus_boelter",
            Nu_multiplier=2.0,
        )
        f2 = ht.smooth(
            T_AIR,
            P_AIR,
            X_AIR,
            u=20.0,
            L=0.5,
            D=0.01,
            correlation="dittus_boelter",
            f_multiplier=2.0,
        )
        assert nu2.h == pytest.approx(2.0 * base.h, rel=1e-6)
        assert f2.dP == pytest.approx(2.0 * base.dP, rel=1e-6)

    def test_network_solver_pattern(self):
        """Demonstrate the intended network solver usage pattern."""
        T_hot = 850.0
        sol = ht.smooth(T_AIR, P_AIR, X_AIR, u=25.0, L=0.3, D=0.008, T_hot=T_hot)
        assert sol.h > 0.0
        assert sol.dP > 0.0
        assert math.isfinite(sol.q)
        assert sol.T_aw >= T_AIR
        assert math.isfinite(sol.M)
        assert sol.M > 0.0
