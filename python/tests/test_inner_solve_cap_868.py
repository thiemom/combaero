"""Inner solves report a cap exit instead of returning a non-root (#868).

The C++ side is pinned in tests/test_inner_solve_cap_868.cpp; these are the
same conditions through the Python API. Each failed on main: the solve
returned a value that was not a root. See docs/archive/AUDIT_868_INNER_SOLVES.md.
"""

import math

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
    (#924 made it a cap exit; since R6 the scan reports no root in range,
    now [1, 6000] K)."""
    with pytest.raises(RuntimeError, match=r"calc_T_from_cp: no root in \[1, 6000\] K"):
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
    # A guess above the range starts the scan at its end.
    assert cb.calc_T_from_cp(target, x, 9000.0) == pytest.approx(upper, abs=1e-9)
    # The tail below the data has a maximum near 127 K, so cp(205 K) is also
    # reached at ~37.4 K, the root nearest a 50 K guess.
    low = cb.calc_T_from_cp(target, x, 50.0)
    assert 37.0 < low < 38.0
    assert cb.cp(low, x) == pytest.approx(target, rel=1e-13)


@pytest.mark.parametrize(("species", "t_true"), [("H2", 50.0), ("C3H8", 1.5), ("NC7H16", 10.0)])
def test_calc_t_from_cp_reaches_the_low_tail(species: str, t_true: float) -> None:
    """Targets only the tail below the 200 K data reaches (the [200, 6000] K
    range threw): the search goes down to CP_INVERSE_T_MIN = 1 K."""
    x = _pure(species)
    assert cb.calc_T_from_cp(cb.cp(t_true, x), x, 300.0) == pytest.approx(t_true, rel=1e-8)


def test_calc_t_from_cp_resolves_roots_below_the_log_scan_switch() -> None:
    """CO2's cp has a minimum at 6.2 K, so cp(4 K) is also reached at ~7.6 K:
    below 100 K the scan steps 10 % of T, which separates the two roots that
    one 10 K step would span."""
    x = _pure("CO2")
    target = cb.cp(4.0, x)
    assert cb.calc_T_from_cp(target, x, 3.0) == pytest.approx(4.0, abs=1e-8)
    assert 7.5 < cb.calc_T_from_cp(target, x, 9.0) < 7.7


def test_calc_t_from_cp_nearest_is_in_kelvin() -> None:
    """O2's cp(1.5 K) is also reached at ~189.9 K. From a 20 K guess the
    1.5 K root is nearer: the scan advances the side whose next node is
    nearer the guess (one step per side per round returned 189.9 K)."""
    x = _pure("O2")
    assert cb.calc_T_from_cp(cb.cp(1.5, x), x, 20.0) == pytest.approx(1.5, abs=1e-8)


@pytest.mark.parametrize(
    ("cp_target", "t_guess"), [(float("nan"), 300.0), (-1.0, 300.0), (30.0, -1.0), (30.0, 0.0)]
)
def test_calc_t_from_cp_refuses_invalid_input(cp_target: float, t_guess: float) -> None:
    """A negative T_guess returned 530 K on main; invalid input is ValueError (#928)."""
    with pytest.raises(ValueError, match="calc_T_from_cp"):
        cb.calc_T_from_cp(cp_target, X_AIR, t_guess)


def _cp_extremum(x: list[float], lo: float, hi: float) -> float:
    """dcp/dT's sign change in [lo, hi], bisected to adjacent doubles."""
    s_lo = cb.dcp_dT(lo, x) > 0.0
    while True:
        mid = 0.5 * (lo + hi)
        if not lo < mid < hi:
            return lo
        if (cb.dcp_dT(mid, x) > 0.0) == s_lo:
            lo = mid
        else:
            hi = mid


def _check_roots(roots: list, target: float, x: list[float]) -> None:
    ts = [r.T for r in roots]
    assert ts == sorted(ts)
    for r in roots:
        if not r.double_root:
            assert abs(cb.cp(r.T, x) - target) <= 1e-14 * target
            assert r.dT_dcp == 1.0 / cb.dcp_dT(r.T, x)


def test_calc_t_from_cp_roots_returns_every_root_sorted() -> None:
    """N2's cp has a tail maximum at 127 K and a minimum at 215 K, so cp(205 K)
    is reached three times; calc_T_from_cp returns one of them (#963)."""
    x = _pure("N2")
    target = cb.cp(205.0, x)
    roots = cb.calc_T_from_cp_roots(target, x)
    assert [r.T for r in roots] == pytest.approx([37.41129, 205.0, 226.73870], abs=1e-5)
    assert not any(r.double_root for r in roots)
    _check_roots(roots, target, x)
    assert cb.calc_T_from_cp(target, x, 213.0) == pytest.approx(roots[1].T, abs=1e-9)
    upper = cb.calc_T_from_cp_roots(target, x, T_min=200.0)
    assert [r.T for r in upper] == pytest.approx([205.0, 226.73870], abs=1e-5)


def test_calc_t_from_cp_roots_finds_a_pair_within_one_step() -> None:
    """O2's cp(127 K) is reached again at 120.8 K, inside one 10 K step:
    calc_T_from_cp sees no sign change and throws; splitting at the 124 K
    minimum finds both."""
    x = _pure("O2")
    target = cb.cp(127.0, x)
    with pytest.raises(RuntimeError, match="no root"):
        cb.calc_T_from_cp(target, x, 300.0)
    roots = cb.calc_T_from_cp_roots(target, x)
    assert [r.T for r in roots] == pytest.approx([120.83142, 127.0], abs=1e-5)
    _check_roots(roots, target, x)


def test_calc_t_from_cp_roots_finds_an_extremum_pair_within_one_step() -> None:
    """N2/O2 at 77.28/22.72 has a maximum at 156.1 K and a minimum at 159.3 K,
    both inside the 150-160 K step, so dcp/dT has one sign on every node: the
    dip of |dcp/dT| is searched, and cp(157.7 K) has three roots."""
    x = [0.0] * cb.num_species()
    x[cb.species_index_from_name("N2")] = 0.7728
    x[cb.species_index_from_name("O2")] = 0.2272
    target = cb.cp(157.7, x)
    roots = cb.calc_T_from_cp_roots(target, x)
    assert [r.T for r in roots] == pytest.approx([154.94551, 157.7, 160.51748], abs=1e-5)
    _check_roots(roots, target, x)


def test_calc_t_from_cp_roots_reports_a_double_root_once() -> None:
    """A target equal to cp at the 127 K maximum of N2, to rounding (here 5e-14
    above it, where strictly there is no root near 127 K), is one double root
    (dT_dcp undefined); a target 1e-9 below it is two roots 0.1 K apart."""
    x = _pure("N2")
    t_max = _cp_extremum(x, 120.0, 135.0)
    for rel in (0.0, 5e-14, -5e-14):
        roots = cb.calc_T_from_cp_roots(cb.cp(t_max, x) * (1.0 + rel), x)
        assert len(roots) == 2
        assert roots[0].double_root
        assert pytest.approx(t_max, abs=1e-9) == roots[0].T
        assert math.isnan(roots[0].dT_dcp)
        assert not roots[1].double_root
        assert pytest.approx(334.16406, abs=1e-5) == roots[1].T
    below = cb.calc_T_from_cp_roots(cb.cp(t_max, x) * (1.0 - 1e-9), x)
    assert [r.T for r in below] == pytest.approx([126.63728, 126.74086, 334.16402], abs=1e-5)
    assert not any(r.double_root for r in below)


def test_calc_t_from_cp_roots_is_empty_without_a_root() -> None:
    """A query: no root in the range is an empty list, not an exception."""
    assert cb.calc_T_from_cp_roots(1.0e3, X_AIR) == []
    assert cb.calc_T_from_cp_roots(20.0, X_AIR) == []
    assert cb.calc_T_from_cp_roots(cb.cp(205.0, _pure("N2")), _pure("N2"), 300.0, 6000.0) == []


@pytest.mark.parametrize(
    ("cp_target", "t_min", "t_max"),
    [
        (float("nan"), 1.0, 6000.0),
        (-1.0, 1.0, 6000.0),
        (30.0, 0.5, 6000.0),
        (30.0, 1.0, 7000.0),
        (30.0, 500.0, 500.0),
        (30.0, float("nan"), 6000.0),
    ],
)
def test_calc_t_from_cp_roots_refuses_invalid_input(
    cp_target: float, t_min: float, t_max: float
) -> None:
    with pytest.raises(ValueError, match="calc_T_from_cp_roots"):
        cb.calc_T_from_cp_roots(cp_target, X_AIR, t_min, t_max)


def test_calc_t_from_cp_roots_refuses_a_constant_cp() -> None:
    """Argon's cp is exactly 5/2 R up to ~950 K: every T there is a root."""
    x = _pure("AR")
    with pytest.raises(ValueError, match="not isolated"):
        cb.calc_T_from_cp_roots(cb.cp(300.0, x), x)


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
