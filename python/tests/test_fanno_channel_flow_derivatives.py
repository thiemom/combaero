"""fanno_channel_flow derivatives: low drive, the drive floor, composition.

At low drive the exit and inlet Mach nearly coincide, (Me - M_in)/Me ~
drive/P, and the implicit derivatives combine terms ~P/drive times larger than
their result. The slopes of the isentropic inlet flux and of the exit pressure
per unit flux were differenced, and that rounding, amplified ~1e5, put up to
8e-3 into dG/dPt0 and dG/dP at a few Pa of drive and 1e-5 noise into G below
the 1 Pa floor (through the blend exponent). Both are exact now.

dG/dY_k (with_dG_dY) holds the other mass fractions fixed, the convention of an
element's "{node}.Y[k]" Jacobian column.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb

N = cb.num_species()
PT, TT, L, D, ROUGH = 2.0e5, 700.0, 1.0, 0.02, 1e-5


def _X(d: dict[str, float]) -> list[float]:
    X = [0.0] * N
    for k, v in d.items():
        X[cb.species_index_from_name(k)] = v
    s = sum(X)
    return [x / s for x in X]


MIXES = {
    "air": list(cb.species.dry_air()),
    "products": _X({"N2": 0.71, "H2O": 0.18, "CO2": 0.09, "O2": 0.02}),
    "co2_rich": _X({"N2": 0.5, "O2": 0.13, "CO2": 0.3, "H2O": 0.07}),
}


def _flow(
    X: list[float],
    P_target: float,
    exit_total: bool = False,
    Pt0: float = PT,
    Tt0: float = TT,
    with_dY: bool = False,
) -> cb.FannoChannelFlow:
    return cb.fanno_channel_flow(
        Pt0, Tt0, X, P_target, exit_total, L, D, ROUGH, "haaland", 1.0, -1.0, with_dY
    )


@pytest.mark.parametrize("drive", [1.5, 3.0, 10.0, 100.0, 1.0e4])
@pytest.mark.parametrize("exit_total", [False, True])
def test_pressure_and_temperature_derivatives_hold_at_low_drive(
    drive: float, exit_total: bool
) -> None:
    """<= 1e-5 down to 1.5 Pa on 2 bar (was up to 8e-3)."""
    X = MIXES["products"]
    P = PT - drive
    f = _flow(X, P, exit_total)
    hp = 2e-3 * drive
    ht = 1e-4 * TT
    fd_pt = (_flow(X, P, exit_total, Pt0=PT + hp).G - _flow(X, P, exit_total, Pt0=PT - hp).G) / (
        2 * hp
    )
    fd_p = (_flow(X, P + hp, exit_total).G - _flow(X, P - hp, exit_total).G) / (2 * hp)
    fd_t = (_flow(X, P, exit_total, Tt0=TT + ht).G - _flow(X, P, exit_total, Tt0=TT - ht).G) / (
        2 * ht
    )
    assert f.dG_dPt0 == pytest.approx(fd_pt, rel=1e-5)
    assert f.dG_dP_target == pytest.approx(fd_p, rel=1e-5)
    assert f.dG_dTt0 == pytest.approx(fd_t, rel=1e-5)


@pytest.mark.parametrize("drive", [0.4, 0.8])
def test_the_flux_is_smooth_below_the_drive_floor(drive: float) -> None:
    """Second differences over 1 mK steps in Tt: rounding level (~1e-9 of G).
    The blend exponent sigma carried the differenced slopes' noise, 1e-5 of G."""
    X = MIXES["air"]
    T = TT + np.linspace(0.0, 0.05, 51)
    G = np.array([_flow(X, PT - drive, Tt0=t).G for t in T])
    assert np.max(np.abs(np.diff(G, 2))) / G.mean() < 1e-7


def _dG_dY_fd(X: list[float], P: float, exit_total: bool, k: int) -> float:
    """Richardson-extrapolated central difference in Y_k (others fixed)."""
    Y = list(cb.mole_to_mass(X))

    def D(h: float) -> float:
        Yp, Ym = list(Y), list(Y)
        Yp[k] += h
        Ym[k] -= h
        Gp = _flow(list(cb.mass_to_mole(Yp)), P, exit_total).G
        Gm = _flow(list(cb.mass_to_mole(Ym)), P, exit_total).G
        return (Gp - Gm) / (2 * h)

    h = min(2e-4, 0.4 * Y[k])
    return (4 * D(h / 2) - D(h)) / 3


@pytest.mark.parametrize("mix", sorted(MIXES))
@pytest.mark.parametrize("P_target", [1.6e5, 1.0e5, 0.5e5], ids=["unchoked", "near", "choked"])
@pytest.mark.parametrize("exit_total", [False, True])
def test_the_composition_derivative_matches_differences(
    mix: str, P_target: float, exit_total: bool
) -> None:
    X = MIXES[mix]
    f = _flow(X, P_target, exit_total, with_dY=True)
    Y = cb.mole_to_mass(X)
    assert len(f.dG_dY) == N
    for k in range(N):
        if Y[k] <= 0.0:
            assert f.dG_dY[k] == 0.0
            continue
        fd = _dG_dY_fd(X, P_target, exit_total, k)
        assert f.dG_dY[k] == pytest.approx(fd, rel=1e-5, abs=1e-6 * f.G)


@pytest.mark.parametrize("drive", [3.0, 30.0])
def test_the_composition_derivative_holds_at_low_drive(drive: float) -> None:
    """The cancellation is ~P/drive there: <= 3e-3 of the column at 3 Pa on
    2 bar, <= 1e-4 at 30 Pa, measured against the column's own scale."""
    X = MIXES["products"]
    P = PT - drive
    f = _flow(X, P, with_dY=True)
    Y = cb.mole_to_mass(X)
    for k in range(N):
        if Y[k] > 0.0:
            fd = _dG_dY_fd(X, P, False, k)
            scale = max(abs(fd), 0.05 * f.G)
            assert abs(f.dG_dY[k] - fd) / scale < (3e-3 if drive < 10 else 1e-4)


def test_the_composition_derivative_is_opt_in() -> None:
    X = MIXES["products"]
    assert len(_flow(X, 1.6e5).dG_dY) == 0
    zero = _flow(X, PT + 10.0, with_dY=True)
    assert zero.G == 0.0 and list(zero.dG_dY) == [0.0] * N
    floor = _flow(X, PT - 0.5, with_dY=True)
    assert len(floor.dG_dY) == N and any(d != 0.0 for d in floor.dG_dY)
