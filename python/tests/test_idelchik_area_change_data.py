"""Idelchik (1966), Handbook of Hydraulic Resistance, diagrams 4-1 (sudden
expansion) and 5-2 (conical diffuser), transcribed from the page images
(#731). The extraction record is
validation/area_change/extractions/idelchik_1966_area_change.md (UNCONFIRMED).

These tests check the transcription against the diagrams' own identities;
scoring the area change element against them comes once it is signed off.
"""

from __future__ import annotations

import csv
from collections import defaultdict

import pytest

from validation.area_change import DATA_ROOT

DATA = DATA_ROOT / "idelchik1966"


def _rows(name: str) -> list[dict[str, float]]:
    with (DATA / name).open() as f:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def test_the_high_re_expansion_is_borda_carnot() -> None:
    """Part 1 prints zeta = (1 - F0/F2)^2 (without baffles) and ~0.6 of it
    (with baffles), both to two decimals."""
    rows = _rows("diag4-1_sudden_expansion_high_re.csv")
    assert len(rows) == 10
    for r in rows:
        borda = (1.0 - r["F0_F2"]) ** 2
        # The printed 0.50 at F0/F2 0.3 is (0.7)^2 = 0.49 rounded up.
        assert r["zeta_no_baffles"] == pytest.approx(borda, abs=0.011), r
        assert r["zeta_baffles"] == pytest.approx(0.6 * borda, abs=0.011), r


def test_the_low_re_table_closes_on_borda_carnot() -> None:
    """Part 2's last column, Re 3.5e3, is where part 1 takes over: it prints
    part 1's own values. 6 area ratios by 13 Reynolds numbers."""
    rows = _rows("diag4-1_sudden_expansion_low_re.csv")
    assert len(rows) == 78
    by_n = defaultdict(dict)
    for r in rows:
        by_n[r["F0_F2"]][r["Re"]] = r["zeta"]
    high = {r["F0_F2"]: r["zeta_no_baffles"] for r in _rows("diag4-1_sudden_expansion_high_re.csv")}
    assert sorted(by_n) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    for n, z in by_n.items():
        assert len(z) == 13
        assert z[3500.0] == high[n], n
        # Re 10: Re alone sets it, every ratio prints 3.10.
        assert z[10.0] == 3.10, n
    # At every Re above 10 a larger area ratio loses less (the curves on
    # graph b never cross).
    for re in sorted({r["Re"] for r in rows} - {10.0}):
        col = [by_n[n][re] for n in sorted(by_n)]
        assert all(a >= b for a, b in zip(col, col[1:], strict=False)), re


def test_the_conical_diffuser_table() -> None:
    """Diagram 5-2: 11 area ratios by 14 full angles. At fixed angle a larger
    ratio loses less; the 180 deg column is close to Borda-Carnot (the
    printed values sit 0.000-0.024 above it, measured); every row peaks at 60 deg, above
    its own 180 deg value (the wide-angle overshoot)."""
    rows = _rows("diag5-2_conical_diffuser.csv")
    assert len(rows) == 154
    tab = defaultdict(dict)
    for r in rows:
        tab[r["F0_F1"]][r["alpha_deg"]] = r["zeta_exp"]
    ratios = sorted(tab)
    assert len(ratios) == 11
    for a in sorted(tab[0.0]):
        col = [tab[n][a] for n in ratios]
        assert all(x >= y for x, y in zip(col, col[1:], strict=False)), a
    for n in ratios:
        assert -1e-9 <= tab[n][180.0] - (1.0 - n) ** 2 <= 0.025, n
        peak = max(tab[n], key=tab[n].get)
        assert peak == 60.0, n
        # Rising with angle up to the peak.
        up = [tab[n][a] for a in sorted(tab[n]) if a <= 60.0]
        assert all(x <= y for x, y in zip(up, up[1:], strict=False)), n
