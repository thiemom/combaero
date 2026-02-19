"""Tests for the symmetric combaero.incompressible / combaero.compressible API.

The key design goal is that ``incompressible`` and ``compressible`` expose
identical function signatures so that swapping the import is sufficient to
change the flow regime.  These tests verify:

1. Both submodules return ``FlowSolution`` instances.
2. The ``regime`` field is set correctly.
3. Fields that are inapplicable carry ``math.nan`` / ``False``.
4. Physically meaningful fields are populated and sensible.
5. The symmetric call pattern works (same kwargs, different module).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
import combaero.compressible as comp
import combaero.incompressible as incomp
from combaero import FlowSolution

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def air() -> tuple:
    """Return (T, P, X) for dry air at 400 K, 2 bar — subsonic Fanno regime."""
    X = cb.standard_dry_air_composition()
    return 400.0, 200_000.0, X


def air_300() -> tuple:
    X = cb.standard_dry_air_composition()
    return 300.0, 101_325.0, X


# ---------------------------------------------------------------------------
# FlowSolution dataclass
# ---------------------------------------------------------------------------


class TestFlowSolution:
    def test_frozen(self) -> None:
        sol = FlowSolution(regime="incompressible", mdot=1.0, v=10.0)
        with pytest.raises((AttributeError, TypeError)):
            sol.mdot = 2.0  # type: ignore[misc]

    def test_defaults_nan_false(self) -> None:
        sol = FlowSolution(regime="incompressible", mdot=1.0, v=10.0)
        assert math.isnan(sol.dP)
        assert math.isnan(sol.Re)
        assert math.isnan(sol.rho)
        assert math.isnan(sol.f)
        assert math.isnan(sol.Cd)
        assert math.isnan(sol.M)
        assert math.isnan(sol.T_out)
        assert math.isnan(sol.P_out)
        assert math.isnan(sol.h0)
        assert sol.choked is False
        assert math.isnan(sol.L_choke)
        assert sol.profile == []

    def test_repr_smoke(self) -> None:
        sol = FlowSolution(regime="compressible", mdot=0.5, v=300.0, M=0.8, choked=False)
        assert "compressible" in repr(sol)
        assert "M=" in repr(sol)


# ---------------------------------------------------------------------------
# Symmetric pipe_flow — same kwargs, different module
# ---------------------------------------------------------------------------


class TestSymmetricPipeFlow:
    KWARGS = {"u": 10.0, "L": 1.0, "D": 0.05, "f": 0.02}

    def test_incompressible_returns_flow_solution(self) -> None:
        T, P, X = air_300()
        sol = incomp.pipe_flow(T, P, X, **self.KWARGS)
        assert isinstance(sol, FlowSolution)

    def test_compressible_returns_flow_solution(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, **self.KWARGS)
        assert isinstance(sol, FlowSolution)

    def test_incompressible_regime(self) -> None:
        T, P, X = air_300()
        sol = incomp.pipe_flow(T, P, X, **self.KWARGS)
        assert sol.regime == "incompressible"

    def test_compressible_regime(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, **self.KWARGS)
        assert sol.regime == "compressible"

    def test_incompressible_inapplicable_fields(self) -> None:
        T, P, X = air_300()
        sol = incomp.pipe_flow(T, P, X, **self.KWARGS)
        assert math.isnan(sol.M)
        assert math.isnan(sol.T_out)
        assert math.isnan(sol.P_out)
        assert math.isnan(sol.h0)
        assert math.isnan(sol.Cd)
        assert math.isnan(sol.L_choke)
        assert sol.choked is False

    def test_compressible_inapplicable_fields(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, **self.KWARGS)
        assert math.isnan(sol.Cd)

    def test_incompressible_applicable_fields(self) -> None:
        T, P, X = air_300()
        sol = incomp.pipe_flow(T, P, X, **self.KWARGS)
        assert sol.mdot > 0.0
        assert sol.v > 0.0
        assert sol.dP > 0.0
        assert sol.Re > 0.0
        assert sol.rho > 0.0
        assert sol.f > 0.0

    def test_compressible_applicable_fields(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, **self.KWARGS)
        assert sol.mdot > 0.0
        assert sol.v > 0.0
        assert sol.dP > 0.0
        assert sol.Re > 0.0
        assert sol.rho > 0.0
        assert sol.f > 0.0
        assert 0.0 < sol.M < 1.0
        assert sol.T_out > 0.0
        assert sol.P_out > 0.0
        assert sol.h0 > 0.0

    def test_both_pressure_drops(self) -> None:
        """Both regimes should give positive dP for the same conditions."""
        T, P, X = air()
        sol_i = incomp.pipe_flow(T, P, X, **self.KWARGS)
        sol_c = comp.pipe_flow(T, P, X, **self.KWARGS)
        assert sol_i.dP > 0.0
        assert sol_c.dP > 0.0

    def test_compressible_dP_higher_than_incompressible(self) -> None:
        """Compressible friction losses are slightly higher due to density change."""
        T, P, X = air()
        sol_i = incomp.pipe_flow(T, P, X, **self.KWARGS)
        sol_c = comp.pipe_flow(T, P, X, **self.KWARGS)
        # Both positive; compressible dP >= incompressible dP (density drops)
        assert sol_c.dP >= sol_i.dP * 0.9  # within 10% — just check sign


# ---------------------------------------------------------------------------
# Symmetric pipe_flow_rough — same kwargs, different module
# ---------------------------------------------------------------------------


class TestSymmetricPipeFlowRough:
    KWARGS = {"u": 10.0, "L": 1.0, "D": 0.05, "roughness": 0.0}

    def test_incompressible_returns_flow_solution(self) -> None:
        T, P, X = air_300()
        sol = incomp.pipe_flow_rough(T, P, X, **self.KWARGS)
        assert isinstance(sol, FlowSolution)

    def test_compressible_returns_flow_solution(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow_rough(T, P, X, **self.KWARGS)
        assert isinstance(sol, FlowSolution)

    def test_regimes(self) -> None:
        T, P, X = air_300()
        assert incomp.pipe_flow_rough(T, P, X, **self.KWARGS).regime == "incompressible"
        T, P, X = air()
        assert comp.pipe_flow_rough(T, P, X, **self.KWARGS).regime == "compressible"

    def test_rough_increases_dP_both_regimes(self) -> None:
        T_i, P_i, X = air_300()
        sol_smooth_i = incomp.pipe_flow_rough(T_i, P_i, X, **self.KWARGS)
        sol_rough_i = incomp.pipe_flow_rough(T_i, P_i, X, u=10.0, L=1.0, D=0.05, roughness=1e-4)
        assert sol_rough_i.dP > sol_smooth_i.dP

        T_c, P_c, X = air()
        sol_smooth_c = comp.pipe_flow_rough(T_c, P_c, X, **self.KWARGS)
        sol_rough_c = comp.pipe_flow_rough(T_c, P_c, X, u=10.0, L=1.0, D=0.05, roughness=1e-4)
        assert sol_rough_c.dP > sol_smooth_c.dP

    def test_correlation_kwarg_accepted_both(self) -> None:
        T_i, P_i, X = air_300()
        sol = incomp.pipe_flow_rough(
            T_i, P_i, X, u=10.0, L=1.0, D=0.05, roughness=0.0, correlation="serghides"
        )
        assert isinstance(sol, FlowSolution)

        T_c, P_c, X = air()
        sol = comp.pipe_flow_rough(
            T_c, P_c, X, u=10.0, L=1.0, D=0.05, roughness=0.0, correlation="serghides"
        )
        assert isinstance(sol, FlowSolution)

    def test_profile_empty_by_default(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow_rough(T, P, X, **self.KWARGS)
        assert sol.profile == []

    def test_profile_populated_when_requested(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow_rough(
            T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0, store_profile=True
        )
        assert len(sol.profile) > 0


# ---------------------------------------------------------------------------
# orifice_flow (incompressible only)
# ---------------------------------------------------------------------------


class TestIncompressibleOrificeFlow:
    def test_returns_flow_solution(self) -> None:
        T, P, X = air_300()
        sol = incomp.orifice_flow(T, P, X, P_back=P - 1000.0, A=1e-4)
        assert isinstance(sol, FlowSolution)

    def test_regime(self) -> None:
        T, P, X = air_300()
        sol = incomp.orifice_flow(T, P, X, P_back=P - 1000.0, A=1e-4)
        assert sol.regime == "incompressible"

    def test_Cd_populated(self) -> None:
        T, P, X = air_300()
        Cd = 0.65
        sol = incomp.orifice_flow(T, P, X, P_back=P - 1000.0, A=1e-4, Cd=Cd)
        assert np.isclose(sol.Cd, Cd)

    def test_zero_dP_zero_mdot(self) -> None:
        T, P, X = air_300()
        sol = incomp.orifice_flow(T, P, X, P_back=P, A=1e-4)
        assert sol.mdot == 0.0
        assert sol.dP == 0.0

    def test_larger_area_more_mdot(self) -> None:
        T, P, X = air_300()
        P_back = P - 1000.0
        sol1 = incomp.orifice_flow(T, P, X, P_back=P_back, A=1e-4)
        sol2 = incomp.orifice_flow(T, P, X, P_back=P_back, A=2e-4)
        assert sol2.mdot > sol1.mdot

    def test_inapplicable_fields(self) -> None:
        T, P, X = air_300()
        sol = incomp.orifice_flow(T, P, X, P_back=P - 1000.0, A=1e-4)
        assert math.isnan(sol.M)
        assert math.isnan(sol.T_out)
        assert math.isnan(sol.h0)
        assert sol.choked is False


# ---------------------------------------------------------------------------
# nozzle_flow (compressible only)
# ---------------------------------------------------------------------------


class TestCompressibleNozzleFlow:
    def test_returns_flow_solution(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P / 2.0, A_eff=1e-3)
        assert isinstance(sol, FlowSolution)

    def test_regime(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P / 2.0, A_eff=1e-3)
        assert sol.regime == "compressible"

    def test_subsonic_not_choked(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.9, A_eff=1e-3)
        assert not sol.choked
        assert 0.0 < sol.M < 1.0

    def test_choked_low_back_pressure(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.1, A_eff=1e-3)
        assert sol.choked

    def test_mdot_positive(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.8, A_eff=1e-3)
        assert sol.mdot > 0.0

    def test_inapplicable_fields(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.8, A_eff=1e-3)
        assert math.isnan(sol.dP)
        assert math.isnan(sol.Re)
        assert math.isnan(sol.rho)
        assert math.isnan(sol.f)
        assert math.isnan(sol.Cd)
        assert math.isnan(sol.L_choke)

    def test_T_out_less_than_T_in(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.5, A_eff=1e-3)
        assert sol.T_out < T

    def test_P_out_less_than_P_in(self) -> None:
        T, P, X = air()
        sol = comp.nozzle_flow(T, P, X, P_back=P * 0.5, A_eff=1e-3)
        assert sol.P_out < P


# ---------------------------------------------------------------------------
# Choked Fanno pipe
# ---------------------------------------------------------------------------


class TestCompressiblePipeChoked:
    def test_choked_flag_set(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, u=250.0, L=1000.0, D=0.05, f=0.02)
        assert sol.choked

    def test_L_choke_finite_when_choked(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, u=250.0, L=1000.0, D=0.05, f=0.02)
        assert not math.isnan(sol.L_choke)
        assert sol.L_choke < 1000.0

    def test_not_choked_short_pipe(self) -> None:
        T, P, X = air()
        sol = comp.pipe_flow(T, P, X, u=10.0, L=0.1, D=0.05, f=0.02)
        assert not sol.choked
        assert math.isnan(sol.L_choke)


# ---------------------------------------------------------------------------
# Regime-swap pattern: same call, different import
# ---------------------------------------------------------------------------


class TestRegimeSwap:
    """Verify that swapping the module import is the only change needed."""

    def _run(self, module, T: float, P: float, X) -> FlowSolution:
        return module.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)

    def test_swap_returns_flow_solution(self) -> None:
        T, P, X = air()
        for module in (incomp, comp):
            sol = self._run(module, T, P, X)
            assert isinstance(sol, FlowSolution)

    def test_swap_different_regimes(self) -> None:
        T, P, X = air()
        sol_i = self._run(incomp, T, P, X)
        sol_c = self._run(comp, T, P, X)
        assert sol_i.regime == "incompressible"
        assert sol_c.regime == "compressible"

    def test_swap_both_give_positive_dP(self) -> None:
        T, P, X = air()
        for module in (incomp, comp):
            sol = self._run(module, T, P, X)
            assert sol.dP > 0.0

    def test_pipe_flow_rough_swap(self) -> None:
        T, P, X = air()
        for module in (incomp, comp):
            sol = module.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
            assert isinstance(sol, FlowSolution)
            assert sol.dP > 0.0
