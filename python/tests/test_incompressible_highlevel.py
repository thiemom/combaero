"""Tests for the thermo-aware incompressible high-level API.

Covers IncompressibleFlowSolution, orifice_flow_thermo, pipe_flow,
pipe_flow_rough, and pressure_drop_pipe (moved from pipe_flow.h).
"""

import math

import numpy as np
import pytest

import combaero as cb

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def air_at_300() -> tuple:
    """Return (T, P, X) for dry air at 300 K, 1 atm."""
    X = cb.standard_dry_air_composition()
    return 300.0, 101325.0, X


# ---------------------------------------------------------------------------
# IncompressibleFlowSolution struct
# ---------------------------------------------------------------------------


class TestIncompressibleFlowSolution:
    def test_pipe_flow_returns_solution(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)
        assert isinstance(sol, cb.IncompressibleFlowSolution)

    def test_pipe_flow_rough_returns_solution(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
        assert isinstance(sol, cb.IncompressibleFlowSolution)

    def test_orifice_flow_thermo_returns_solution(self) -> None:
        T, P, X = air_at_300()
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P - 1000.0, A=1e-4)
        assert isinstance(sol, cb.IncompressibleFlowSolution)


# ---------------------------------------------------------------------------
# pipe_flow (constant f)
# ---------------------------------------------------------------------------


class TestPipeFlow:
    def test_basic_pressure_drop(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)
        assert sol.dP > 0.0

    def test_zero_velocity_zero_dP(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow(T, P, X, u=0.0, L=1.0, D=0.05, f=0.02)
        assert sol.dP == 0.0
        assert sol.mdot == 0.0

    def test_re_populated(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)
        assert sol.Re > 0.0

    def test_rho_matches_density(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)
        rho_ref = cb.density(T, P, X)
        assert np.isclose(sol.rho, rho_ref, rtol=1e-6)

    def test_f_matches_input(self) -> None:
        T, P, X = air_at_300()
        f_in = 0.025
        sol = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=f_in)
        assert np.isclose(sol.f, f_in)

    def test_mdot_consistent_with_velocity(self) -> None:
        T, P, X = air_at_300()
        u = 10.0
        D = 0.05
        sol = cb.pipe_flow(T, P, X, u=u, L=1.0, D=D, f=0.02)
        A = math.pi * D**2 / 4.0
        mdot_ref = sol.rho * u * A
        assert np.isclose(sol.mdot, mdot_ref, rtol=1e-6)

    def test_dP_matches_darcy_weisbach(self) -> None:
        T, P, X = air_at_300()
        u, L, D, f = 10.0, 1.0, 0.05, 0.02
        sol = cb.pipe_flow(T, P, X, u=u, L=L, D=D, f=f)
        dP_ref = cb.pipe_dP(u, L, D, f, sol.rho)
        assert np.isclose(sol.dP, dP_ref, rtol=1e-6)

    def test_longer_pipe_higher_dP(self) -> None:
        T, P, X = air_at_300()
        sol1 = cb.pipe_flow(T, P, X, u=10.0, L=1.0, D=0.05, f=0.02)
        sol2 = cb.pipe_flow(T, P, X, u=10.0, L=2.0, D=0.05, f=0.02)
        assert sol2.dP > sol1.dP

    def test_invalid_T(self) -> None:
        _, P, X = air_at_300()
        with pytest.raises((ValueError, RuntimeError)):
            cb.pipe_flow(-1.0, P, X, u=10.0, L=1.0, D=0.05, f=0.02)


# ---------------------------------------------------------------------------
# pipe_flow_rough (roughness-based f)
# ---------------------------------------------------------------------------


class TestPipeFlowRough:
    def test_basic_pressure_drop(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
        assert sol.dP > 0.0

    def test_re_populated(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
        assert sol.Re > 0.0

    def test_f_positive(self) -> None:
        T, P, X = air_at_300()
        sol = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
        assert sol.f > 0.0

    def test_rough_higher_dP_than_smooth(self) -> None:
        T, P, X = air_at_300()
        sol_s = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0)
        sol_r = cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=1e-4)
        assert sol_r.dP > sol_s.dP

    def test_smooth_matches_pipe_flow_with_haaland_f(self) -> None:
        """pipe_flow_rough(roughness=0) should match pipe_flow(f=haaland(Re))."""
        T, P, X = air_at_300()
        u, L, D = 10.0, 1.0, 0.05
        Re = cb.reynolds(T, P, X, u, D)
        f_ref = cb.friction_haaland(Re, 0.0)

        sol_f = cb.pipe_flow(T, P, X, u=u, L=L, D=D, f=f_ref)
        sol_r = cb.pipe_flow_rough(T, P, X, u=u, L=L, D=D, roughness=0.0)

        assert np.isclose(sol_f.dP, sol_r.dP, rtol=1e-6)

    def test_correlations_consistent(self) -> None:
        T, P, X = air_at_300()
        kw = {"u": 10.0, "L": 1.0, "D": 0.05, "roughness": 0.0}
        sol_h = cb.pipe_flow_rough(T, P, X, **kw, correlation="haaland")
        sol_s = cb.pipe_flow_rough(T, P, X, **kw, correlation="serghides")
        sol_c = cb.pipe_flow_rough(T, P, X, **kw, correlation="colebrook")
        assert np.isclose(sol_h.dP, sol_s.dP, rtol=0.01)
        assert np.isclose(sol_h.dP, sol_c.dP, rtol=0.01)

    def test_invalid_correlation(self) -> None:
        T, P, X = air_at_300()
        with pytest.raises(ValueError, match="unknown friction correlation"):
            cb.pipe_flow_rough(T, P, X, u=10.0, L=1.0, D=0.05, roughness=0.0, correlation="bad")


# ---------------------------------------------------------------------------
# orifice_flow_thermo
# ---------------------------------------------------------------------------


class TestOrificeFlowThermo:
    def test_basic_mdot(self) -> None:
        T, P, X = air_at_300()
        P_back = P - 1000.0
        A = 1e-4
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P_back, A=A)
        assert sol.mdot > 0.0

    def test_zero_dP_zero_mdot(self) -> None:
        T, P, X = air_at_300()
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P, A=1e-4)
        assert sol.mdot == 0.0
        assert sol.dP == 0.0

    def test_rho_matches_density(self) -> None:
        T, P, X = air_at_300()
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P - 500.0, A=1e-4)
        rho_ref = cb.density(T, P, X)
        assert np.isclose(sol.rho, rho_ref, rtol=1e-6)

    def test_f_field_stores_Cd(self) -> None:
        T, P, X = air_at_300()
        Cd = 0.65
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P - 1000.0, A=1e-4, Cd=Cd)
        assert np.isclose(sol.f, Cd)

    def test_mdot_matches_scalar_primitive(self) -> None:
        T, P, X = air_at_300()
        P_back = P - 1000.0
        A = 1e-4
        Cd = 0.65
        rho = cb.density(T, P, X)
        mdot_ref = cb.orifice_mdot(P, P_back, A, Cd, rho)
        sol = cb.orifice_flow_thermo(T, P, X, P_back=P_back, A=A, Cd=Cd)
        assert np.isclose(sol.mdot, mdot_ref, rtol=1e-6)

    def test_larger_area_more_mdot(self) -> None:
        T, P, X = air_at_300()
        P_back = P - 1000.0
        sol1 = cb.orifice_flow_thermo(T, P, X, P_back=P_back, A=1e-4)
        sol2 = cb.orifice_flow_thermo(T, P, X, P_back=P_back, A=2e-4)
        assert sol2.mdot > sol1.mdot


# ---------------------------------------------------------------------------
# pressure_drop_pipe (legacy tuple API, now in incompressible.cpp)
# ---------------------------------------------------------------------------


class TestPressureDropPipeRelocated:
    def test_returns_tuple(self) -> None:
        T, P, X = air_at_300()
        result = cb.pressure_drop_pipe(T, P, X, v=10.0, D=0.05, L=1.0)
        dP, Re, f = result
        assert dP > 0.0
        assert Re > 0.0
        assert f > 0.0

    def test_matches_pipe_flow_rough(self) -> None:
        T, P, X = air_at_300()
        u, D, L = 10.0, 0.05, 1.0
        dP_tuple, Re_tuple, f_tuple = cb.pressure_drop_pipe(T, P, X, v=u, D=D, L=L)
        sol = cb.pipe_flow_rough(T, P, X, u=u, L=L, D=D, roughness=0.0)
        assert np.isclose(dP_tuple, sol.dP, rtol=1e-6)
        assert np.isclose(Re_tuple, sol.Re, rtol=1e-6)
        assert np.isclose(f_tuple, sol.f, rtol=1e-6)


# ---------------------------------------------------------------------------
# Scalar primitives (zeta/Cd, pressure_loss)
# ---------------------------------------------------------------------------


class TestScalarPrimitives:
    def test_zeta_from_Cd_roundtrip(self) -> None:
        Cd = 0.65
        zeta = cb.zeta_from_Cd(Cd)
        Cd_back = cb.Cd_from_zeta(zeta)
        assert np.isclose(Cd_back, Cd, rtol=1e-10)

    def test_zeta_from_Cd_value(self) -> None:
        assert np.isclose(cb.zeta_from_Cd(1.0), 1.0)
        assert np.isclose(cb.zeta_from_Cd(0.5), 4.0)

    def test_pressure_loss_consistency(self) -> None:
        v, rho, zeta = 10.0, 1.2, 2.5
        dP = cb.pressure_loss(v, rho, zeta)
        v_back = cb.velocity_from_pressure_loss(dP, rho, zeta)
        assert np.isclose(v_back, v, rtol=1e-10)

    def test_dynamic_pressure(self) -> None:
        v, rho = 10.0, 1.2
        q = cb.dynamic_pressure(v, rho)
        assert np.isclose(q, 0.5 * rho * v**2)
        v_back = cb.velocity_from_q(q, rho)
        assert np.isclose(v_back, v)
