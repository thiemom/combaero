"""#893: compression and expansion along a path of given efficiency --
cb.compression / cb.expansion, their closed forms, their limits and their
Jacobian."""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb

AIR = list(cb.species.dry_air())
AR = list(cb.species.pure_species("AR"))
CASES = [
    (cb.compression, 288.15, 10.0),
    (cb.compression, 400.0, 2.0),
    (cb.expansion, 1400.0, 10.0),
    (cb.expansion, 900.0, 2.5),
]


def _products() -> list[float]:
    """CH4 burnt lean in air: a turbine's gas."""
    Y = list(cb.mole_to_mass(AIR))
    Y_f = list(cb.species.to_mass(cb.species.pure_species("CH4")))
    st = cb.mix(
        [_stream(700.0, Y, 1.0), _stream(300.0, Y_f, 0.02)],
    )
    return list(cb.complete_combustion(st.T, list(st.X), st.P).X)


def _stream(T: float, Y: list[float], m: float) -> cb.Stream:
    s = cb.Stream()
    s.set_T(T)
    s.set_P(1.0e6)
    s.set_X(list(cb.mass_to_mole(Y)))
    s.set_mdot(m)
    return s


@pytest.mark.parametrize("eta", [1.0, 0.9, 0.6])
@pytest.mark.parametrize("PR", [1.5, 10.0, 40.0])
def test_argon_has_the_constant_gamma_closed_forms(PR: float, eta: float) -> None:
    """Argon's cp is 5/2 R: (gamma - 1)/gamma = 0.4, and the textbook forms
    are exact. Its data hold 2.5 R to 6e-9 below 3000 K (measured), hence 1e-8."""
    T1 = 300.0
    c = cb.compression(T1, AR, PR, eta_poly=eta)
    assert pytest.approx(T1 * PR ** (0.4 / eta), rel=1e-8) == c.T2
    c = cb.compression(T1, AR, PR, eta_is=eta)
    assert pytest.approx(T1 * (1.0 + (PR**0.4 - 1.0) / eta), rel=1e-8) == c.T2
    T1 = 1500.0
    e = cb.expansion(T1, AR, PR, eta_poly=eta)
    assert pytest.approx(T1 * PR ** (-0.4 * eta), rel=1e-8) == e.T2
    e = cb.expansion(T1, AR, PR, eta_is=eta)
    assert pytest.approx(T1 * (1.0 - eta * (1.0 - PR**-0.4)), rel=1e-8) == e.T2


@pytest.mark.parametrize(("fn", "T1", "PR"), CASES)
def test_the_polytropic_path_is_the_limit_of_small_isentropic_stages(fn, T1, PR) -> None:
    """The definition of the polytropic efficiency: n small stages, each of
    isentropic efficiency eta_p, independent of the closed form. They approach
    it at first order (the gap falls 4.00x per 4x stages, measured); the
    Richardson limit of n = 1000 and 4000 removes that term."""
    eta = 0.88

    def staged(n: int) -> float:
        T = T1
        for _ in range(n):
            T = fn(T, AIR, PR ** (1.0 / n), eta_is=eta).T2
        return T

    exact = fn(T1, AIR, PR, eta_poly=eta).T2
    T1k, T4k = staged(1000), staged(4000)
    assert (T1k - exact) / (T4k - exact) == pytest.approx(4.0, rel=1e-3)
    assert (4.0 * T4k - T1k) / 3.0 == pytest.approx(exact, rel=1e-8)


@pytest.mark.parametrize(("fn", "T1", "PR"), CASES)
def test_each_efficiency_round_trips_through_the_other(fn, T1, PR) -> None:
    p = fn(T1, AIR, PR, eta_poly=0.9)
    i = fn(T1, AIR, PR, eta_is=p.eta_is)
    assert pytest.approx(p.T2, rel=1e-10) == i.T2
    assert i.eta_poly == pytest.approx(0.9, rel=1e-9)
    # compression: eta_is < eta_poly; expansion: eta_is > eta_poly
    assert (p.eta_is < 0.9) if fn is cb.compression else (p.eta_is > 0.9)


@pytest.mark.parametrize(("fn", "T1", "PR"), CASES)
def test_ideal_paths_coincide_and_lose_no_entropy(fn, T1, PR) -> None:
    a = fn(T1, AIR, PR, eta_poly=1.0)
    b = fn(T1, AIR, PR, eta_is=1.0)
    assert a.T2 == pytest.approx(b.T2, rel=1e-10) == pytest.approx(a.T2s, rel=1e-10)
    P1 = 1.0e5
    P2 = P1 * PR if fn is cb.compression else P1 / PR
    for eta in (1.0, 0.95, 0.7):
        r = fn(T1, AIR, PR, eta_poly=eta)
        ds = cb.s_mass(r.T2, AIR, P2) - cb.s_mass(T1, AIR, P1)
        assert ds >= -1e-9 * abs(cb.s_mass(T1, AIR, P1)), (eta, ds)
        assert (ds > 0.0) if eta < 1.0 else True


def test_unit_pressure_ratio_is_no_change() -> None:
    for fn in (cb.compression, cb.expansion):
        r = fn(500.0, AIR, 1.0, eta_poly=0.8)
        assert r.T2 == 500.0
        assert r.dh == 0.0
        assert r.eta_poly == r.eta_is == 0.8


def test_work_signs_and_magnitudes() -> None:
    c = cb.compression(288.15, AIR, 10.0, eta_poly=0.9)
    e = cb.expansion(1400.0, AIR, 10.0, eta_poly=0.9)
    assert c.dh > 0.0 > e.dh
    assert c.dh == pytest.approx(cb.h_mass(c.T2, AIR) - cb.h_mass(288.15, AIR), rel=1e-12)


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"eta_poly": 1.05}, "outside"),
        ({"eta_is": 0.0}, "outside"),
        ({"eta_poly": 0.9, "eta_is": 0.9}, "exactly one"),
        ({}, "exactly one"),
    ],
)
def test_invalid_efficiencies_are_refused(kw: dict, match: str) -> None:
    for fn in (cb.compression, cb.expansion):
        with pytest.raises(ValueError, match=match):
            fn(300.0, AIR, 2.0, **kw)


def test_a_pressure_ratio_below_one_is_refused() -> None:
    with pytest.raises(ValueError, match="< 1"):
        cb.compression(300.0, AIR, 0.5, eta_poly=0.9)


# -- Jacobian -----------------------------------------------------------------


def _fd(f, x: float, h: float) -> float:
    return (f(x + h) - f(x - h)) / (2.0 * h)


def _check(analytic: float, numeric: float, scale: float, what: str) -> None:
    assert analytic == pytest.approx(numeric, rel=2e-6, abs=1e-9 * scale), what


@pytest.mark.parametrize("gas", ["air", "products"])
@pytest.mark.parametrize("mode", ["eta_poly", "eta_is"])
@pytest.mark.parametrize(("fn", "T1", "PR"), CASES)
def test_the_jacobian_is_exact(fn, T1, PR, mode: str, gas: str) -> None:
    X = AIR if gas == "air" else _products()
    Y = np.asarray(cb.mole_to_mass(X))
    eta = 0.87
    r = fn(T1, X, PR, **{mode: eta})

    def at(T1_: float = T1, PR_: float = PR, eta_: float = eta, Y_: np.ndarray = Y):
        return fn(T1_, list(cb.mass_to_mole(list(Y_))), PR_, **{mode: eta_})

    # Mass-fraction directions along the simplex: zero sum, present species only.
    rng = np.random.default_rng(7)
    dirs = []
    for _ in range(3):
        v = rng.standard_normal(len(Y)) * (Y > 1e-6)
        dirs.append(v - (Y > 1e-6) * v.sum() / max((Y > 1e-6).sum(), 1))

    for attr, d in (("T2", "dT2"), ("dh", "ddh")):
        scale = max(abs(getattr(r, attr)), 1.0)
        _check(
            getattr(r, f"{d}_dT1"),
            _fd(lambda v, a=attr: getattr(at(T1_=v), a), T1, 1e-3),
            scale,
            f"{attr}/T1",
        )
        _check(
            getattr(r, f"{d}_dPR"),
            _fd(lambda v, a=attr: getattr(at(PR_=v), a), PR, 1e-5 * PR),
            scale,
            f"{attr}/PR",
        )
        _check(
            getattr(r, f"{d}_deta"),
            _fd(lambda v, a=attr: getattr(at(eta_=v), a), eta, 1e-6),
            scale,
            f"{attr}/eta",
        )
        grad = np.asarray(getattr(r, f"{d}_dY"))
        for v in dirs:
            _check(
                float(grad @ v),
                _fd(lambda t, a=attr, v=v: getattr(at(Y_=Y + t * v), a), 0.0, 1e-7),
                scale,
                f"{attr}/Y",
            )


def test_the_temperature_does_not_depend_on_the_pressure_level() -> None:
    """At a fixed PR the outlet T is the same at any inlet pressure: the
    kernel takes none, and a direct isentrope at two levels agrees."""
    r = cb.compression(300.0, AIR, 8.0, eta_is=1.0)
    for P1 in (1.0e4, 1.0e5, 3.0e6):
        s1 = cb.s(300.0, AIR, P1)
        assert cb.calc_T_from_s(s1, 8.0 * P1, AIR, 500.0) == pytest.approx(r.T2, rel=1e-9)
    assert math.isfinite(r.T2)
