"""Tests for incompressible flow functions."""

import numpy as np
import pytest
import combaero as cb


# Water at 20°C
RHO_WATER = 998.0  # kg/m³


class TestBernoulli:
    """Tests for Bernoulli equation."""

    def test_bernoulli_P2_basic(self) -> None:
        """Pressure drops when velocity increases."""
        P1 = 200000.0  # Pa
        v1 = 2.0  # m/s
        v2 = 5.0  # m/s
        rho = RHO_WATER

        P2 = cb.bernoulli_P2(P1, v1, v2, rho)

        # P2 < P1 because v2 > v1
        assert P2 < P1
        # Check formula: P2 = P1 + 0.5*rho*(v1² - v2²)
        expected = P1 + 0.5 * rho * (v1**2 - v2**2)
        assert np.isclose(P2, expected)

    def test_bernoulli_P2_with_elevation(self) -> None:
        """Pressure drops with elevation increase."""
        P1 = 200000.0
        v1 = v2 = 2.0  # Same velocity
        rho = RHO_WATER
        dz = 10.0  # 10 m up

        P2 = cb.bernoulli_P2(P1, v1, v2, rho, dz=dz)

        # P2 < P1 due to elevation
        assert P2 < P1
        # ΔP ≈ ρgh = 998 * 9.81 * 10 ≈ 97.9 kPa
        assert np.isclose(P1 - P2, rho * 9.80665 * dz)

    def test_bernoulli_v2_basic(self) -> None:
        """Velocity increases when pressure drops."""
        P1 = 200000.0
        P2 = 150000.0
        v1 = 2.0
        rho = RHO_WATER

        v2 = cb.bernoulli_v2(P1, P2, v1, rho)

        assert v2 > v1

    def test_bernoulli_roundtrip(self) -> None:
        """P2 -> v2 -> P2 roundtrip."""
        P1 = 200000.0
        v1 = 2.0
        v2 = 5.0
        rho = RHO_WATER

        P2 = cb.bernoulli_P2(P1, v1, v2, rho)
        v2_back = cb.bernoulli_v2(P1, P2, v1, rho)

        assert np.isclose(v2_back, v2)

    def test_bernoulli_v2_insufficient_energy(self) -> None:
        """Raises when energy insufficient for flow."""
        with pytest.raises(ValueError):
            # P2 > P1 and v1 = 0 -> impossible
            cb.bernoulli_v2(P1=100000, P2=200000, v1=0.0, rho=RHO_WATER)


class TestOrifice:
    """Tests for orifice flow."""

    def test_orifice_mdot_basic(self) -> None:
        """Basic orifice mass flow."""
        P1 = 200000.0
        P2 = 100000.0
        A = 0.001  # m²
        Cd = 0.62
        rho = RHO_WATER

        mdot = cb.orifice_mdot(P1, P2, A, Cd, rho)

        assert mdot > 0
        # Check formula: mdot = Cd * A * sqrt(2 * rho * dP)
        expected = Cd * A * np.sqrt(2 * rho * (P1 - P2))
        assert np.isclose(mdot, expected)

    def test_orifice_Q_vs_mdot(self) -> None:
        """Q = mdot / rho."""
        P1, P2 = 200000.0, 100000.0
        A, Cd, rho = 0.001, 0.62, RHO_WATER

        mdot = cb.orifice_mdot(P1, P2, A, Cd, rho)
        Q = cb.orifice_Q(P1, P2, A, Cd, rho)

        assert np.isclose(Q, mdot / rho)

    def test_orifice_area_inverse(self) -> None:
        """Solve for area, then verify mdot."""
        P1, P2 = 200000.0, 100000.0
        Cd, rho = 0.62, RHO_WATER
        mdot_target = 5.0  # kg/s

        A = cb.orifice_area(mdot_target, P1, P2, Cd, rho)
        mdot_check = cb.orifice_mdot(P1, P2, A, Cd, rho)

        assert np.isclose(mdot_check, mdot_target)

    def test_orifice_dP_inverse(self) -> None:
        """Solve for dP, then verify mdot."""
        A, Cd, rho = 0.001, 0.62, RHO_WATER
        mdot = 5.0

        dP = cb.orifice_dP(mdot, A, Cd, rho)
        P1, P2 = 200000.0, 200000.0 - dP
        mdot_check = cb.orifice_mdot(P1, P2, A, Cd, rho)

        assert np.isclose(mdot_check, mdot)

    def test_orifice_velocity(self) -> None:
        """Ideal orifice velocity."""
        P1, P2 = 200000.0, 100000.0
        rho = RHO_WATER

        v = cb.orifice_velocity(P1, P2, rho)

        expected = np.sqrt(2 * (P1 - P2) / rho)
        assert np.isclose(v, expected)

    def test_orifice_invalid_pressure(self) -> None:
        """Raises when P1 < P2."""
        with pytest.raises(ValueError):
            cb.orifice_mdot(P1=100000, P2=200000, A=0.001, Cd=0.62, rho=RHO_WATER)


class TestPipe:
    """Tests for pipe flow."""

    def test_pipe_dP_basic(self) -> None:
        """Basic pipe pressure drop."""
        v = 2.0  # m/s
        L = 10.0  # m
        D = 0.05  # m
        f = 0.02
        rho = RHO_WATER

        dP = cb.pipe_dP(v, L, D, f, rho)

        # Check formula: dP = f * (L/D) * (rho * v² / 2)
        expected = f * (L / D) * 0.5 * rho * v**2
        assert np.isclose(dP, expected)

    def test_pipe_dP_mdot_consistency(self) -> None:
        """pipe_dP and pipe_dP_mdot should agree."""
        L, D, f, rho = 10.0, 0.05, 0.02, RHO_WATER
        v = 2.0

        mdot = cb.pipe_mdot(v, D, rho)
        dP_from_v = cb.pipe_dP(v, L, D, f, rho)
        dP_from_mdot = cb.pipe_dP_mdot(mdot, L, D, f, rho)

        assert np.isclose(dP_from_v, dP_from_mdot)

    def test_pipe_velocity_mdot_roundtrip(self) -> None:
        """velocity -> mdot -> velocity roundtrip."""
        D, rho = 0.05, RHO_WATER
        v = 3.0

        mdot = cb.pipe_mdot(v, D, rho)
        v_back = cb.pipe_velocity(mdot, D, rho)

        assert np.isclose(v_back, v)


class TestHydraulic:
    """Tests for hydraulic utilities."""

    def test_dynamic_pressure(self) -> None:
        """q = 0.5 * rho * v²."""
        v, rho = 10.0, RHO_WATER

        q = cb.dynamic_pressure(v, rho)

        assert np.isclose(q, 0.5 * rho * v**2)

    def test_velocity_from_q_roundtrip(self) -> None:
        """v -> q -> v roundtrip."""
        v, rho = 10.0, RHO_WATER

        q = cb.dynamic_pressure(v, rho)
        v_back = cb.velocity_from_q(q, rho)

        assert np.isclose(v_back, v)

    def test_hydraulic_diameter_circle(self) -> None:
        """For circle, Dh = D."""
        D = 0.1
        A = np.pi * D**2 / 4
        P_wetted = np.pi * D

        Dh = cb.hydraulic_diameter(A, P_wetted)

        assert np.isclose(Dh, D)

    def test_hydraulic_diameter_rect(self) -> None:
        """Rectangular duct Dh = 2ab/(a+b)."""
        a, b = 0.1, 0.05

        Dh = cb.hydraulic_diameter_rect(a, b)

        expected = 2 * a * b / (a + b)
        assert np.isclose(Dh, expected)

    def test_hydraulic_diameter_square(self) -> None:
        """Square duct: Dh = a."""
        a = 0.1

        Dh = cb.hydraulic_diameter_rect(a, a)

        assert np.isclose(Dh, a)

    def test_hydraulic_diameter_annulus(self) -> None:
        """Annulus Dh = D_outer - D_inner."""
        D_outer, D_inner = 0.1, 0.05

        Dh = cb.hydraulic_diameter_annulus(D_outer, D_inner)

        assert np.isclose(Dh, D_outer - D_inner)
