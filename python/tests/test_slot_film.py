"""Slot film effectiveness, Marek & Tacina (1975) NASA TN D-7958 (#714).

The kernel is Eq. (4), eta = 1/(1 + C_m x/(M S)), with C_m = 0.25 Tu/100
(Eq. 5). Fidelity is judged on the report's own Table I (the transcription
is CONFIRMED, validation/cooling/extractions/marek_tacina_1975_d7958.md):

- with the authors' per-blockage C_m (p. 8) the kernel must reproduce the
  +-30% deviation the report states for that case;
- with C_m from Eq. (5) at the measured Tu the score is pinned and judged
  against the authors' +-0.04 eta measurement band (p. 6).

The C_m_multiplier tuner is the user's; nothing here scores it.
"""

from __future__ import annotations

import csv
import math

import pytest

import combaero as cb
from validation.cooling.schema import DATA_ROOT

DATA = DATA_ROOT / "marek_tacina1975"
H_DUCT_CM, S_SLOT_CM = 10.5, 0.63  # report p. 3
# p. 8: the average C_m per blockage percent.
C_M_AUTHORS = {"0": 0.025, "52": 0.04, "72": 0.05, "90": 0.1}


def _points() -> list[tuple[str, float, float, float, float]]:
    """(blockage, x/S, M, Tu, measured eta) for every cooled, printed cell."""
    with (DATA / "blockage_turbulence.csv").open() as f:
        tu = {r["blockage_pct"]: float(r["Tu_axial_pct"]) for r in csv.DictReader(f)}
    out = []
    with (DATA / "table1_film_temperatures.csv").open() as f:
        for r in csv.DictReader(f):
            m_s = float(r["m_cool_kg_s"])
            if m_s <= 0.0 or not r["T_wall_K"]:
                continue
            T_H, T_s, T_w = float(r["T_hot_K"]), float(r["T_cool_K"]), float(r["T_wall_K"])
            M = m_s / float(r["m_hot_kg_s"]) * H_DUCT_CM / S_SLOT_CM
            b = r["blockage_pct"]
            out.append((b, float(r["x_cm"]) / S_SLOT_CM, M, tu[b], (T_H - T_w) / (T_H - T_s)))
    return out


def _quiet(*args: float, **kw: float) -> float:
    with cb.suppress_warnings():
        return cb.slot_film_effectiveness_marek_tacina_1975(*args, **kw)


def _warnings_of(*args: float) -> list[str]:
    messages: list[str] = []
    cb.set_warning_handler(messages.append)
    try:
        cb.slot_film_effectiveness_marek_tacina_1975(*args)
    finally:
        cb.set_warning_handler(None)
    return messages


def test_eq4_and_eq5_closed_form() -> None:
    """eta = 1/(1 + 0.25 Tu/100 x/(M S)); 1 at the lip."""
    x_S, M, Tu = 15.0, 1.2, 14.0
    expected = 1.0 / (1.0 + 0.25 * Tu / 100.0 * x_S / M)
    assert cb.slot_film_effectiveness_marek_tacina_1975(x_S, M, Tu) == pytest.approx(
        expected, rel=1e-14
    )
    assert cb.slot_film_effectiveness_marek_tacina_1975(0.0, M, Tu) == 1.0


def test_the_multiplier_scales_c_m_only() -> None:
    """C_m_multiplier k is the same as Tu scaled by k: it acts on C_m alone."""
    a = cb.slot_film_effectiveness_marek_tacina_1975(20.0, 0.8, 14.0, C_m_multiplier=1.5)
    b = cb.slot_film_effectiveness_marek_tacina_1975(20.0, 0.8, 21.0)
    assert a == pytest.approx(b, rel=1e-14)


def test_the_derivative_is_the_central_difference() -> None:
    for x_S, M, Tu in ((2.0, 0.5, 7.0), (15.0, 1.0, 14.0), (30.0, 1.8, 35.0)):
        eta, d = cb.slot_film_effectiveness_marek_tacina_1975_and_derivatives(x_S, M, Tu)
        assert eta == cb.slot_film_effectiveness_marek_tacina_1975(x_S, M, Tu)
        h = 1e-6 * M
        fd = (
            cb.slot_film_effectiveness_marek_tacina_1975(x_S, M + h, Tu)
            - cb.slot_film_effectiveness_marek_tacina_1975(x_S, M - h, Tu)
        ) / (2 * h)
        assert d == pytest.approx(fd, rel=1e-7)
        assert d > 0.0  # more coolant, more film


@pytest.mark.parametrize(
    ("args", "kw"),
    [
        ((-1.0, 1.0, 14.0), {}),
        ((5.0, 0.0, 14.0), {}),
        ((5.0, 1.0, -1.0), {}),
        ((5.0, 1.0, 14.0), {"C_m_multiplier": 0.0}),
    ],
)
def test_refuses_what_the_form_cannot_represent(
    args: tuple[float, ...], kw: dict[str, float]
) -> None:
    with pytest.raises(ValueError):
        cb.slot_film_effectiveness_marek_tacina_1975(*args, **kw)


def test_outside_the_run_it_warns_and_extrapolates() -> None:
    for args, name in (
        ((10.0, 1.0, 3.0), "Tu"),
        ((10.0, 3.0, 14.0), "M"),
        ((50.0, 1.0, 14.0), "x_over_S"),
    ):
        (msg,) = _warnings_of(*args)
        assert f"{name} = " in msg and "Extrapolating" in msg


def test_inside_the_run_it_is_silent() -> None:
    """Every Table I point is inside the declared span: no warning."""
    for _, x_S, M, Tu, _ in _points():
        assert _warnings_of(x_S, M, Tu) == []


def test_reproduces_the_reports_30_percent_with_its_own_c_m() -> None:
    """FIDELITY. p. 8: "If a constant value of C_m is taken for each blockage
    plate, the deviation between the correlation and the film-cooling
    effectiveness is +-30 percent." With the printed per-blockage C_m (as a
    multiplier on Eq. 5's line at that Tu) the kernel gives 136 of 140 points
    inside 30%, the worst 33% -- the statement, to its rounding.
    """
    rel = []
    for b, x_S, M, Tu, eta in _points():
        k = C_M_AUTHORS[b] / (0.25 * Tu / 100.0)
        rel.append(_quiet(x_S, M, Tu, C_m_multiplier=k) / eta - 1.0)
    assert len(rel) == 140
    assert sum(abs(r) <= 0.30 for r in rel) == 136
    assert max(abs(r) for r in rel) == pytest.approx(0.33, abs=0.005)


def test_eq5_at_the_measured_turbulence_scores_on_table_one() -> None:
    """FIDELITY of Eq. (5) on its own data, against the +-0.04 eta band.

    Eq. (5) is a line through four per-blockage averages, and C_m varies
    +-50% along the plate (p. 8), so most points sit outside a +-0.04 band:
    57 of 140, eta RMS 0.0696, bias +0.006. Per blockage the bias changes
    sign (+0.068 at 0%, -0.038 at 72%): the line's slope, not an offset.
    Pinned as measured; a change to the kernel or the data moves it.
    """
    pts = _points()
    err = [_quiet(x_S, M, Tu) - eta for _, x_S, M, Tu, eta in pts]
    rms = math.sqrt(sum(e * e for e in err) / len(err))
    assert rms == pytest.approx(0.0696, abs=2e-4)
    assert sum(err) / len(err) == pytest.approx(0.006, abs=1e-3)
    assert sum(abs(e) <= 0.04 for e in err) == 57
    by = {b: [e for (bb, *_), e in zip(pts, err, strict=True) if bb == b] for b in C_M_AUTHORS}
    assert sum(by["0"]) / len(by["0"]) == pytest.approx(0.068, abs=1e-3)
    assert sum(by["72"]) / len(by["72"]) == pytest.approx(-0.038, abs=1e-3)
