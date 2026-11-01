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
    """main returned the 6000 K clamp for a cp no T in [200, 6000] K has."""
    with pytest.raises(RuntimeError, match="calc_T_from_cp: no convergence"):
        cb.calc_T_from_cp(1.0e3, X_AIR, 300.0)
    assert cb.calc_T_from_cp(cb.cp(800.0, X_AIR), X_AIR, 300.0) == pytest.approx(800.0, abs=1e-6)


def test_rayleigh_choke_search_without_a_sign_change_raises() -> None:
    """An unphysical exit Tt (a solver probe in the C3X vane example) left the
    choke search unbracketed; main returned G = 6.3e-5, choked, as its root."""
    args = (266266.38, 1603.23)
    tail = (X_AIR, 136765.9, False, 0.0142, 0.0897, 0.0, "haaland", 1.0, 0.0636, [], False)
    with pytest.raises(RuntimeError, match="no sign change"):
        cb.fanno_rayleigh_channel_flow(*args, 211975.1, *tail)
    ok = cb.fanno_rayleigh_channel_flow(*args, 2000.0, *tail)
    assert ok.G > 0.0 and not ok.choked
