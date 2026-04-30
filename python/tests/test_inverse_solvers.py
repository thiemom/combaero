"""Round-trip tests for thermodynamic inverse solvers.

Each test computes a thermodynamic property forward (T -> property),
then inverts it (property -> T) and checks the recovered temperature.
"""

import pytest

import combaero as cb

AIR = cb.species.dry_air()
TOL = 1e-3  # K — tighter than the solver default (1e-6 K) to catch regressions


@pytest.fixture
def air():
    return AIR


# ---------------------------------------------------------------------------
# Molar-basis inverses
# ---------------------------------------------------------------------------


def test_calc_T_from_h(air):
    T_ref = 850.0
    h_val = cb.h(T_ref, air)
    T_rec = cb.calc_T_from_h(h_target=h_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_s(air):
    T_ref = 600.0
    P = 300_000.0
    s_val = cb.s(T_ref, air, P)
    T_rec = cb.calc_T_from_s(s_target=s_val, P=P, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_cp(air):
    T_ref = 700.0
    cp_val = cb.cp(T_ref, air)
    T_rec = cb.calc_T_from_cp(cp_target=cp_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_u(air):
    T_ref = 900.0
    u_val = cb.u(T_ref, air)
    T_rec = cb.calc_T_from_u(u_target=u_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


# ---------------------------------------------------------------------------
# Mass-specific inverses
# ---------------------------------------------------------------------------


def test_calc_T_from_h_mass(air):
    T_ref = 850.0
    h_val = cb.h_mass(T_ref, air)
    T_rec = cb.calc_T_from_h_mass(h_mass_target=h_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_s_mass(air):
    T_ref = 600.0
    P = 300_000.0
    s_val = cb.s_mass(T_ref, air, P)
    T_rec = cb.calc_T_from_s_mass(s_mass_target=s_val, P=P, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_u_mass(air):
    T_ref = 900.0
    u_val = cb.u_mass(T_ref, air)
    T_rec = cb.calc_T_from_u_mass(u_mass_target=u_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


# ---------------------------------------------------------------------------
# Flash calculations
# ---------------------------------------------------------------------------


def test_calc_T_from_sv_mass(air):
    T_ref = 700.0
    P_ref = 200_000.0
    rho = cb.density(T_ref, P_ref, air)
    v_val = 1.0 / rho
    s_val = cb.s_mass(T_ref, air, P_ref)
    T_rec = cb.calc_T_from_sv_mass(s_mass_target=s_val, v_mass_target=v_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


def test_calc_T_from_sh_mass(air):
    T_ref = 800.0
    P_ref = 101_325.0
    s_val = cb.s_mass(T_ref, air, P_ref)
    h_val = cb.h_mass(T_ref, air)
    T_rec = cb.calc_T_from_sh_mass(s_mass_target=s_val, h_mass_target=h_val, X=air, T_guess=300.0)
    assert abs(T_rec - T_ref) < TOL


# ---------------------------------------------------------------------------
# API hygiene: all inverse solvers must be importable and callable
# ---------------------------------------------------------------------------


def test_all_inverse_solvers_are_exported():
    expected = [
        "calc_T_from_h",
        "calc_T_from_h_mass",
        "calc_T_from_s",
        "calc_T_from_s_mass",
        "calc_T_from_cp",
        "calc_T_from_u",
        "calc_T_from_u_mass",
        "calc_T_from_sv_mass",
        "calc_T_from_sh_mass",
    ]
    for name in expected:
        assert hasattr(cb, name), f"cb.{name} not found"
        assert callable(getattr(cb, name)), f"cb.{name} is not callable"
