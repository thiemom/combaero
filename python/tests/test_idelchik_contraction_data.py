"""Idelchik (1966), Handbook of Hydraulic Resistance, diagrams 3-6 (conical
converging bellmouth with end wall: the zeta' of diagram 3-7's in-line conical
contraction) and 3-10 (sudden contraction, Re < 1e4, Karev), transcribed from
the page images (#731). The extraction record is
validation/area_change/extractions/idelchik_1966_contraction.md (CONFIRMED).

These tests check the transcription against the diagrams' own identities;
the closures score against them in test_area_change.py.
"""

from __future__ import annotations

import csv
from collections import defaultdict

from validation.area_change import DATA_ROOT

DATA = DATA_ROOT / "idelchik1966"


def _rows(name: str) -> list[dict[str, float]]:
    with (DATA / name).open() as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def test_the_conical_bellmouth_table() -> None:
    """3-6: 6 lengths by 9 full angles. Both ends are the sharp edge (0.50:
    no cone at 0 deg, a flat wall at 180 deg), and at every angle between a
    longer cone loses no more."""
    rows = _rows("diag3-6_conical_bellmouth_end_wall.csv")
    assert len(rows) == 54
    tab: dict[float, dict[float, float]] = defaultdict(dict)
    for r in rows:
        tab[r["alpha_deg"]][r["l_Dh"]] = r["zeta_prime"]
    assert sorted(tab) == [0, 10, 20, 30, 40, 60, 100, 140, 180]
    for a in (0.0, 180.0):
        assert set(tab[a].values()) == {0.5}
    for a, col in tab.items():
        z = [col[length] for length in sorted(col)]
        assert all(x >= y for x, y in zip(z, z[1:], strict=False)), a
        assert all(0.0 < x <= 0.5 for x in z), a


def test_the_low_re_contraction_table() -> None:
    """3-10: 6 area ratios by 14 printed columns. Both of the last two are
    headed 10^4; the last is 3-9's 0.5 (1 - n) exactly, the high-Re limit.
    At Re 10 every ratio prints 5.00, and at every Re a larger ratio loses
    no more."""
    rows = _rows("diag3-10_sudden_contraction_low_re.csv")
    assert len(rows) == 84
    by_c: dict[int, dict[float, float]] = defaultdict(dict)
    for r in rows:
        by_c[int(r["column"])][r["F0_F1"]] = r["zeta"]
    assert len(by_c) == 14
    for n, z in by_c[14].items():
        assert z == round(0.5 * (1.0 - n), 2), n
    assert set(by_c[1].values()) == {5.0}
    for c, col in by_c.items():
        z = [col[n] for n in sorted(col)]
        assert all(x >= y for x, y in zip(z, z[1:], strict=False)), c
    # The two "10^4" columns differ at n 0.1 only (0.50 against 0.45).
    assert {n for n in by_c[13] if by_c[13][n] != by_c[14][n]} == {0.1}
