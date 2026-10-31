"""Reader-Harris/Gallagher below its lower Re_D bound (#873).

ISO 5167-2 Eq. (4) was evaluated far below its Re_D range: the Re terms grow
without bound as Re_D falls, so a 4 mm, beta 0.5 plate read Cd 1.16 at Re_D
100 and hit the 1.5 cap below Re_D ~ 40, and a Re_D < 1 floor sat under it.
Both were derivative kinks. Below the bound Re_b = metering_cd_re_min the Cd
is now HELD, entered C1 in ln Re_D over w = 0.1 just below Re_b, and
metering_cd_in_range keeps flagging it.

Stolz has NO hold: ISO 5167:1980's lower bound is not held as a source
(validation/sources.yaml, wanted; follow-up of #873). It is pinned here to
its pre-#873 values, floor and cap included.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver, OrificeElement, PressureBoundary

M = cb.MeteringCdCorrelation
T = cb.OrificeTap
W = 0.1  # orifice::metering_hold::blend_ln_re


def _geom(beta: float, D: float, tap: object) -> cb.OrificeGeometry:
    g = cb.OrificeGeometry()
    g.D = D
    g.d = beta * D
    g.tap = tap
    return g


def _cd(corr: object, g: cb.OrificeGeometry, Re: float) -> tuple[float, float]:
    c, d = cb.metering_cd_and_dRe(corr, g, Re)
    return float(c), float(d)


CASES = [
    (M.ReaderHarrisGallagher, tap, beta, D)
    for tap in (T.Corner, T.D_D2, T.Flange)
    for beta in (0.2, 0.5, 0.6, 0.75)
    for D in (0.004, 0.1, 1.0)
]


def test_the_lower_bounds_are_the_standards() -> None:
    rhg = M.ReaderHarrisGallagher
    # ISO 5167-2:2003 5.3.1
    assert cb.metering_cd_re_min(rhg, _geom(0.5, 0.1, T.Corner)) == 5000.0
    assert cb.metering_cd_re_min(rhg, _geom(0.7, 0.1, T.D_D2)) == pytest.approx(16000 * 0.49)
    assert cb.metering_cd_re_min(rhg, _geom(0.2, 0.1, T.Flange)) == 5000.0
    assert cb.metering_cd_re_min(rhg, _geom(0.7, 0.5, T.Flange)) == pytest.approx(170 * 0.49 * 500)
    with pytest.raises(ValueError, match="Reynolds range"):
        cb.metering_cd_re_min(M.Constant, _geom(0.5, 0.1, T.Corner))


@pytest.mark.parametrize("tap", [T.Corner, T.D_D2, T.Flange])
def test_stolz_has_no_encoded_bound(tap: object) -> None:
    """ISO 5167:1980 is not held, so no bound is encoded: asking for one is
    an error, not a number from memory."""
    with pytest.raises(ValueError, match="Stolz has no encoded Reynolds range"):
        cb.metering_cd_re_min(M.Stolz, _geom(0.5, 0.1, tap))


@pytest.mark.parametrize(("corr", "tap", "beta", "D"), CASES)
def test_cd_is_held_below_the_bound(corr: object, tap: object, beta: float, D: float) -> None:
    """Below Re_b e^-w the value is constant, the slope exactly 0, and the
    held value is C_b - s w / 2 (s = Re_b dC/dRe_D at the bound). On the old
    code Cd kept rising to the 1.5 cap."""
    g = _geom(beta, D, tap)
    Rb = cb.metering_cd_re_min(corr, g)
    C_b, dC_b = _cd(corr, g, Rb)
    held = C_b - 0.5 * dC_b * Rb * W
    for Re in (Rb * math.exp(-W) * 0.999, Rb * 0.5, Rb * 1e-3, 1.0, 0.0, -5.0):
        c, d = _cd(corr, g, Re)
        assert c == pytest.approx(held, abs=1e-13), Re
        assert d == 0.0, Re
    # monotone: the hold sits above the bound value, never below it
    assert held >= C_b


@pytest.mark.parametrize(("corr", "tap", "beta", "D"), CASES)
def test_the_blend_is_c1_and_its_derivative_exact(
    corr: object, tap: object, beta: float, D: float
) -> None:
    g = _geom(beta, D, tap)
    Rb = cb.metering_cd_re_min(corr, g)
    s = abs(_cd(corr, g, Rb)[1]) * Rb
    # Cd and Re dCd/dRe continuous at both ends of the blend
    for edge in (Rb, Rb * math.exp(-W)):
        lo, hi = _cd(corr, g, edge * (1 - 1e-11)), _cd(corr, g, edge * (1 + 1e-11))
        assert lo[0] == pytest.approx(hi[0], abs=3e-11 * s + 1e-14)
        assert lo[1] * edge == pytest.approx(hi[1] * edge, abs=1e-8 * s)
    # analytic vs central FD inside the blend and around the bound
    for u in (-0.09, -0.05, -0.01, -1e-4, 1e-4, 0.05):
        Re = Rb * math.exp(u)
        h = 1e-6 * Re
        fd = (_cd(corr, g, Re + h)[0] - _cd(corr, g, Re - h)[0]) / (2 * h)
        assert _cd(corr, g, Re)[1] * Re == pytest.approx(fd * Re, abs=1e-7 * s), u
        # the slope ramps from s at the bound to 0 at its lower edge
        if u < 0:
            assert _cd(corr, g, Re)[1] * Re == pytest.approx(
                _cd(corr, g, Rb)[1] * Rb * (1 + u / W), abs=1e-9 * s
            )


@pytest.mark.parametrize(("corr", "tap", "beta", "D"), CASES)
def test_at_and_above_the_bound_nothing_moved(
    corr: object, tap: object, beta: float, D: float
) -> None:
    """At Re_D >= Re_b the correlation is evaluated as printed: the
    independent ISO Eq. (4) form here, to the last bits."""
    g = _geom(beta, D, tap)
    Rb = cb.metering_cd_re_min(corr, g)
    for Re in (Rb, Rb * 1.01, Rb * 10.0, 1e7):
        c = _cd(corr, g, Re)[0]
        L1, L2 = _tap(tap, D)
        assert c == pytest.approx(_printed(corr, beta, Re, D, L1, L2), abs=1e-13), Re


def _tap(tap: object, D: float) -> tuple[float, float]:
    if tap == T.Corner:
        return 0.0, 0.0
    if tap == T.D_D2:
        return 1.0, 0.47
    L = 25.4 / (D * 1000.0)
    return L, L


def _printed(corr: object, b: float, Re: float, D: float, L1: float, L2: float) -> float:
    """ISO 5167-2 Eq. (4) as printed, typed independently of the C++."""
    assert corr == M.ReaderHarrisGallagher
    A = (19000.0 * b / Re) ** 0.8
    M2 = 2.0 * L2 / (1.0 - b)
    C = (
        0.5961
        + 0.0261 * b**2
        - 0.216 * b**8
        + 0.000521 * (1e6 * b / Re) ** 0.7
        + (0.0188 + 0.0063 * A) * b**3.5 * (1e6 / Re) ** 0.3
        + (0.043 + 0.080 * math.exp(-10.0 * L1) - 0.123 * math.exp(-7.0 * L1))
        * (1.0 - 0.11 * A)
        * b**4
        / (1.0 - b**4)
        - 0.031 * (M2 - 0.8 * M2**1.1) * b**1.3
    )
    if D * 1000.0 < 71.12:
        C += 0.011 * (0.75 - b) * (2.8 - D * 1000.0 / 25.4)
    return C


def test_the_issue_example() -> None:
    """D 4 mm, beta 0.5, D-D/2 taps: 0.687 at Re_D 1000 and the 1.5 cap at
    Re_D 1-40 before; the held value 0.6345 now, everywhere below 0.905 Re_b."""
    g = _geom(0.5, 0.004, T.D_D2)
    for Re in (1.0, 10.0, 40.0, 100.0, 300.0, 1000.0):
        assert _cd(M.ReaderHarrisGallagher, g, Re)[0] == pytest.approx(0.6345, abs=5e-5)
    assert _cd(M.ReaderHarrisGallagher, g, 1e4)[0] == pytest.approx(0.625, abs=5e-4)  # unchanged


def test_the_hold_stays_inside_half_the_iso_uncertainty() -> None:
    """The width w = 0.1 was picked once from this bound (orifice.h,
    metering_hold): over every in-range RHG geometry the held value is within
    half of ISO 5167-2 5.3.3.1's stated uncertainty of C at the bound."""
    worst = 0.0
    for tap in (T.Corner, T.D_D2, T.Flange):
        for D in (0.05, 0.1, 0.3, 1.0):
            for i in range(66):
                beta = 0.1 + 0.01 * i
                g = _geom(beta, D, tap)
                Rb = cb.metering_cd_re_min(M.ReaderHarrisGallagher, g)
                C_b = _cd(M.ReaderHarrisGallagher, g, Rb)[0]
                held = _cd(M.ReaderHarrisGallagher, g, 0.5 * Rb)[0]
                u = (0.7 - beta) if beta < 0.2 else (0.5 if beta <= 0.6 else 1.667 * beta - 0.5)
                worst = max(worst, (held - C_b) / C_b * 100.0 / u)
    assert 0.0 < worst <= 0.5


def test_the_flag_still_says_out_of_range_below_the_bound() -> None:
    corr = M.ReaderHarrisGallagher
    for tap in (T.Corner, T.D_D2, T.Flange):
        g = _geom(0.5, 0.1, tap)
        Rb = cb.metering_cd_re_min(corr, g)
        assert cb.metering_cd_in_range(corr, g, Rb)
        assert not cb.metering_cd_in_range(corr, g, Rb * 0.999)
        assert not cb.metering_cd_in_range(corr, g, 0.0)
        # Stolz's 1980 limits are not encoded: the flag stays True.
        assert cb.metering_cd_in_range(M.Stolz, g, 0.0)


# metering_cd_and_dRe(Stolz, ...) on main before #873 (af22aedc), repr'd:
# the Re_D < 1 floor (value of Re_D = 1, slope 0), the 1.5 cap, and the
# equation itself, unchanged. (tap, beta, D, Re_D, Cd, dCd/dRe_D)
STOLZ_PRE_873 = [
    (T.Corner, 0.1, 0.1, 0.0, 0.886160312783276, 0.0),
    (T.Corner, 0.1, 0.1, 0.5, 0.886160312783276, 0.0),
    (T.Corner, 0.1, 0.1, 1.0, 0.886160312783276, -0.21750936316053165),
    (T.Corner, 0.1, 0.1, 10.0, 0.6477201515024408, -0.0038679242199905197),
    (T.Corner, 0.1, 0.1, 1000.0, 0.5977786886161806, -1.2231450352100224e-06),
    (T.Corner, 0.5, 0.1, 10.0, 1.5, 0.0),
    (T.Corner, 0.5, 0.1, 100.0, 1.1151333928752514, -0.0038450586415594823),
    (T.Corner, 0.5, 0.1, 1000.0, 0.6936267555030199, -6.837588612677457e-05),
    (T.Corner, 0.5, 0.1, 100000.0, 0.6053418878298926, -2.1622353719291627e-08),
    (T.D_D2, 0.5, 0.004, 1000.0, 0.69424688050302, -6.837588612677457e-05),
    (T.D_D2, 0.5, 0.004, 100000.0, 0.6059620128298927, -2.1622353719291627e-08),
    (T.Flange, 0.6, 0.05, 1000.0, 0.7494032055346991, -0.00010785892383659225),
    (T.Flange, 0.6, 0.05, 100000.0, 0.6101390386232195, -3.410798652982583e-08),
]


@pytest.mark.parametrize(("tap", "beta", "D", "Re", "cd", "dcd"), STOLZ_PRE_873)
def test_stolz_is_unchanged(
    tap: object, beta: float, D: float, Re: float, cd: float, dcd: float
) -> None:
    assert _cd(M.Stolz, _geom(beta, D, tap), Re) == (cd, dcd)


def test_a_low_flow_orifice_reads_the_held_cd() -> None:
    """The network effect (#873): a plenum-fed 0.1 m RHG bore (D = 1 m, beta
    0.1) at 1 Pa runs at Re_D ~ 420, below its bound of 5000. It read Cd
    0.6173 from the extrapolated Re terms (and a small piped plate climbed to
    the 1.5 cap); it now reads the held value, 0.6005, still flagged."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.0e5 + 1.0, Tt=300.0))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(
        OrificeElement("o", "hi", "ho", diameter=0.1, correlation="ReaderHarrisGallagher")
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r["__message__"]
    geom = g.elements["o"]._orifice_geom
    Re_D = 4.0 * r["o.m_dot"] / (math.pi * geom.D * r["o.mu_ref"])
    Rb = cb.metering_cd_re_min(M.ReaderHarrisGallagher, geom)
    assert 1.0 < Re_D < Rb * math.exp(-W)
    assert r["o.Cd"] == pytest.approx(_cd(M.ReaderHarrisGallagher, geom, 0.5 * Rb)[0], abs=1e-12)
    assert r["o.Cd"] == pytest.approx(0.6005, abs=5e-5)
    assert r["o.Cd_in_range"] == 0.0
