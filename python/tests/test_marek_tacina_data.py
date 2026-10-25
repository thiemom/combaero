"""Marek & Tacina (1975), NASA TN D-7958, Table I (slot film cooling at four
free-stream turbulence levels), transcribed from the page image (#714). The
extraction record is
validation/cooling/extractions/marek_tacina_1975_d7958.md (UNCONFIRMED).

These tests check the transcription against the report's own identities;
nothing scores against it until the record is signed off.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict

from validation.cooling.schema import DATA_ROOT

DATA = DATA_ROOT / "marek_tacina1975"
H_DUCT_CM, S_SLOT_CM = 10.5, 0.63  # report p. 3


def _rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def test_the_table_shape() -> None:
    """19 runs (blockage 0/52/72/90: 5/4/5/5) by 9 stations; the 72% plate's
    16.20 cm station is printed '---' in every run."""
    rows = _rows("table1_film_temperatures.csv")
    assert len(rows) == 171
    runs = defaultdict(set)
    for r in rows:
        runs[r["blockage_pct"]].add(r["run"])
    assert {b: len(v) for b, v in runs.items()} == {"0": 5, "52": 4, "72": 5, "90": 5}
    missing = [(r["blockage_pct"], r["x_cm"]) for r in rows if not r["T_wall_K"]]
    assert missing == [("72", "16.20")] * 5


def test_the_cooling_percent_is_of_the_total_flow() -> None:
    """Printed percent = 100 m_s / (m_s + m_H) to the printed rounding
    (measured: worst 0.043)."""
    for r in _rows("table1_film_temperatures.csv"):
        mc, mh = float(r["m_cool_kg_s"]), float(r["m_hot_kg_s"])
        assert abs(100 * mc / (mc + mh) - float(r["cool_pct"])) < 0.05, r


def test_more_cooling_gives_a_cooler_wall() -> None:
    """At every blockage and station the wall cools as the slot flow rises."""
    by = defaultdict(list)
    for r in _rows("table1_film_temperatures.csv"):
        if r["T_wall_K"]:
            by[(r["blockage_pct"], r["x_cm"])].append(
                (float(r["m_cool_kg_s"]), float(r["T_wall_K"]))
            )
    for key, v in by.items():
        t = [tw for _, tw in sorted(v)]
        assert all(a > b for a, b in zip(t, t[1:], strict=False)), key


def test_equation_4_on_table_i_returns_the_reports_mixing_coefficients() -> None:
    """eta = (T_H - T_w)/(T_H - T_s) (Eq. 1, wall = film temperature) and
    Eq. 4, eta = 1/(1 + C_m x/(M S)), M = (m_s/m_H)(H/S), give each cooled
    point's C_m. Their mean per blockage lands inside the report's own
    average and band (p. 8): measured 2026-10-25 0.026, 0.037, 0.046, 0.079
    against 0.025 +-0.01, 0.04 +-0.02, 0.05 +-0.03, 0.1 +-0.05."""
    cm = defaultdict(list)
    for r in _rows("table1_film_temperatures.csv"):
        mc = float(r["m_cool_kg_s"])
        if mc == 0.0 or not r["T_wall_K"]:
            continue
        Th, Tc, Tw = float(r["T_hot_K"]), float(r["T_cool_K"]), float(r["T_wall_K"])
        eta = (Th - Tw) / (Th - Tc)
        M = mc / float(r["m_hot_kg_s"]) * H_DUCT_CM / S_SLOT_CM
        cm[r["blockage_pct"]].append((1.0 / eta - 1.0) * M * S_SLOT_CM / float(r["x_cm"]))
    for b in _rows("blockage_turbulence.csv"):
        mean = statistics.mean(cm[b["blockage_pct"]])
        assert abs(mean - float(b["C_m_avg"])) <= float(b["C_m_band"]), (b, mean)
