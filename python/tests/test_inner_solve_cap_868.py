"""Inner solves report a cap exit instead of returning a non-root (#868).

The C++ side is pinned in tests/test_inner_solve_cap_868.cpp; these are the
same conditions through the Python API. Each failed on main: the solve
returned a value that was not a root. See docs/archive/AUDIT_868_INNER_SOLVES.md.
"""

import pytest

import combaero as cb

X_AIR = list(cb.dry_air())


@pytest.mark.parametrize("ratio", [1e-12, 1e-20, 1e-25])
def test_fanno_inlet_mach_is_linear_at_vanishing_flux(ratio: float) -> None:
    """Below ~1e-10 of the sonic flux the slope stencil crossed M = 0 and
    Newton 2-cycled: M_in was 36-73 % high on main."""
    g_max = cb.fanno_sonic_mass_flux(2.0e5, 300.0, X_AIR)
    m_ref = cb.fanno_inlet_mach(2.0e5, 300.0, 1e-8 * g_max, X_AIR)
    m = cb.fanno_inlet_mach(2.0e5, 300.0, ratio * g_max, X_AIR)
    assert m / (m_ref * ratio / 1e-8) == pytest.approx(1.0, abs=1e-9)


def test_calc_t_from_cp_refuses_an_unreachable_target() -> None:
    """The 6000 K clamp was returned for a cp no T in [200, 6000] K has
    (#924 made it a cap exit; since R6 the scan reports no root in range)."""
    with pytest.raises(RuntimeError, match=r"calc_T_from_cp: no root in \[200, 6000\] K"):
        cb.calc_T_from_cp(1.0e3, X_AIR, 300.0)
    with pytest.raises(RuntimeError, match="no root"):
        cb.calc_T_from_cp(20.0, X_AIR, 300.0)
    assert cb.calc_T_from_cp(cb.cp(800.0, X_AIR), X_AIR, 300.0) == pytest.approx(800.0, abs=1e-9)


def _pure(name: str) -> list[float]:
    x = [0.0] * cb.num_species()
    x[cb.species_index_from_name(name)] = 1.0
    return x


def test_calc_t_from_cp_takes_the_root_nearest_the_guess() -> None:
    """cp of N2 has a minimum at ~214 K, so cp(205 K) is also reached at
    ~226.7 K: the root is the one nearest T_guess (#868 R6)."""
    x = _pure("N2")
    target = cb.cp(205.0, x)
    lower = cb.calc_T_from_cp(target, x, 213.0)
    upper = cb.calc_T_from_cp(target, x, 218.0)
    assert lower == pytest.approx(205.0, abs=1e-9)
    assert 226.0 < upper < 227.5
    assert cb.cp(upper, x) == pytest.approx(target, rel=1e-13)
    # A guess outside the range starts the scan at its end.
    assert cb.calc_T_from_cp(target, x, 50.0) == pytest.approx(205.0, abs=1e-9)
    assert cb.calc_T_from_cp(target, x, 9000.0) == pytest.approx(upper, abs=1e-9)


@pytest.mark.parametrize(
    ("cp_target", "t_guess"), [(float("nan"), 300.0), (-1.0, 300.0), (30.0, -1.0), (30.0, 0.0)]
)
def test_calc_t_from_cp_refuses_invalid_input(cp_target: float, t_guess: float) -> None:
    """A negative T_guess returned 530 K on main; invalid input is ValueError (#928)."""
    with pytest.raises(ValueError, match="calc_T_from_cp"):
        cb.calc_T_from_cp(cp_target, X_AIR, t_guess)


@pytest.mark.parametrize("mach", [1.7782794100389228, 3.1622776601683795, 5.623413251903491])
def test_t_from_stagnation_meets_the_energy_equation(mach: float) -> None:
    """The exact slope (T dgamma/dT kept) converges quadratically: main's
    approximate slope stopped with |F| up to 8e-10 cp T here (#868 R3)."""
    tt = 2500.0
    t = cb.T_from_stagnation(tt, mach, X_AIR)
    v = mach * cb.speed_of_sound(t, X_AIR)
    f = cb.h_mass(t, X_AIR) + 0.5 * v * v - cb.h_mass(tt, X_AIR)
    assert abs(f) / (cb.cp_mass(t, X_AIR) * t) < 1e-13


def test_set_vh_mass_is_the_enthalpy_inverse() -> None:
    """u + R T = h: set_VH_mass is calc_T_from_h_mass from the current T
    (its hand Newton stopped at 1e-4 K, #868 R5)."""
    st = cb.State()
    st.X = X_AIR
    st.TP = (1234.5, 3.0e5)
    h, v = st.h, 1.0 / st.rho
    st.TP = (400.0, 1.0e5)
    st.set_VH_mass(v, h)
    assert cb.calc_T_from_h_mass(h, X_AIR, 400.0) == st.T
    assert pytest.approx(3.0e5, rel=1e-9) == st.P


def test_rayleigh_choke_search_without_a_sign_change_raises() -> None:
    """An unphysical exit Tt (a solver probe in the C3X vane example) left the
    choke search unbracketed; main returned G = 6.3e-5, choked, as its root."""
    args = (266266.38, 1603.23)
    tail = (X_AIR, 136765.9, False, 0.0142, 0.0897, 0.0, "haaland", 1.0, 0.0636, [], False)
    with pytest.raises(RuntimeError, match="no sign change"):
        cb.fanno_rayleigh_channel_flow(*args, 211975.1, *tail)
    ok = cb.fanno_rayleigh_channel_flow(*args, 2000.0, *tail)
    assert ok.G > 0.0 and not ok.choked
