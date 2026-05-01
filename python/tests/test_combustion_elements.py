"""Tests for combustion elements and functions."""

import pytest

import combaero as cb
from combaero.network import (
    CombustionResult,
    combustion_from_phi,
    mix_streams,
)
from combaero.network.components import NetworkMixtureState


def test_combustion_result_dataclass():
    """Test CombustionResult dataclass creation."""
    result = CombustionResult(
        X=[0.1] * 14,
        Y=[0.1] * 14,
        mw=29.0,
        T=1000.0,
        P=101325.0,
        m_dot=0.1,
        h=50000.0,
        cp=1000.0,
        rho=1.0,
        gamma=1.4,
        a=300.0,
        phi=1.0,
        T_adiabatic=2000.0,
        eta=0.95,
        Q_released=1000.0,
    )

    assert result.T == 1000.0
    assert result.phi == 1.0
    assert result.Q_released == 1000.0


def test_mix_streams_basic():
    """Test basic stream mixing functionality."""
    state1 = NetworkMixtureState(
        P=101325.0,
        Pt=101325.0,
        T=300.0,
        Tt=300.0,
        m_dot=0.05,
        Y=cb.species.empty(),
    )
    state1.Y[0] = 1.0

    state2 = NetworkMixtureState(
        P=101325.0,
        Pt=101325.0,
        T=400.0,
        Tt=400.0,
        m_dot=0.05,
        Y=cb.species.empty(),
    )
    state2.Y[1] = 1.0

    result = mix_streams(state1, state2)

    assert result.m_dot == 0.1
    assert pytest.approx(101325.0, rel=1e-6) == result.P
    assert result.Y[0] == pytest.approx(0.5, abs=1e-5)
    assert result.Y[1] == pytest.approx(0.5, abs=1e-5)
    assert sum(result.X) == pytest.approx(1.0, abs=1e-6)
    assert 300.0 < result.T < 400.0


def test_mix_streams_energy_conservation():
    """Test that mixing identical streams gives the same enthalpy."""
    state1 = NetworkMixtureState(
        P=101325.0,
        Pt=101325.0,
        T=300.0,
        Tt=300.0,
        m_dot=0.05,
        Y=cb.species.dry_air_mass(),
    )

    state2 = NetworkMixtureState(
        P=101325.0,
        Pt=101325.0,
        T=300.0,
        Tt=300.0,
        m_dot=0.05,
        Y=cb.species.dry_air_mass(),
    )

    result = mix_streams(state1, state2)

    expected_h = cb.h_mass(300.0, list(cb.mass_to_mole(list(state1.Y))))
    assert result.h == pytest.approx(expected_h, rel=1e-6)


def test_combustion_from_phi_basic():
    """Test combustion_from_phi with basic parameters."""
    air_state = NetworkMixtureState(
        P=500000.0,
        Pt=500000.0,
        T=700.0,
        Tt=700.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )

    X_ch4 = cb.species.pure_species("CH4")

    result = combustion_from_phi(air_state, X_ch4, phi=1.0)

    assert result.m_dot > air_state.m_dot
    assert result.T > air_state.T
    assert result.phi > 0.0
    assert result.T_adiabatic == pytest.approx(result.T, rel=1e-3)
    assert result.Q_released > 0


def test_combustion_from_phi_lean_vs_rich():
    """Test lean vs rich combustion."""
    air_state = NetworkMixtureState(
        P=500000.0,
        Pt=500000.0,
        T=700.0,
        Tt=700.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )

    X_ch4 = cb.species.pure_species("CH4")

    result_lean = combustion_from_phi(air_state, X_ch4, phi=0.8)
    result_rich = combustion_from_phi(air_state, X_ch4, phi=1.2)

    assert result_lean.phi < 1.0
    assert result_rich.phi > 1.0
    assert result_rich.m_dot > result_lean.m_dot


def test_combustion_efficiency():
    """Test combustion efficiency parameter."""
    air_state = NetworkMixtureState(
        P=500000.0,
        Pt=500000.0,
        T=700.0,
        Tt=700.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )

    X_ch4 = cb.species.pure_species("CH4")

    result_full = combustion_from_phi(air_state, X_ch4, phi=1.0, eta=1.0)
    result_partial = combustion_from_phi(air_state, X_ch4, phi=1.0, eta=0.8)

    assert result_partial.T < result_full.T


def test_pressure_drop():
    """Test pressure drop parameter."""
    air_state = NetworkMixtureState(
        P=500000.0,
        Pt=500000.0,
        T=700.0,
        Tt=700.0,
        m_dot=0.1,
        Y=cb.species.dry_air_mass(),
    )

    X_ch4 = cb.species.pure_species("CH4")

    result_no_drop = combustion_from_phi(air_state, X_ch4, phi=1.0, delta_P_frac=0.0)
    result_drop = combustion_from_phi(air_state, X_ch4, phi=1.0, delta_P_frac=0.04)

    assert result_drop.P < result_no_drop.P
    assert pytest.approx(result_no_drop.P * 0.96, rel=1e-3) == result_drop.P
