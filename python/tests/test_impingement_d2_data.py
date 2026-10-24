"""CR-3630 Table D.2 (Florschuetz, Metzger and Su, 1982): heat transfer with
an initial crossflow, Nu and eta per segment (eta referenced to the entrance
crossflow temperature, Eq. 5.3), transcribed from the text layer (#467).
The extraction record is
validation/cooling/extractions/florschuetz_metzger_su_1982_cr3630_table_d2.md.
These tests check the transcription against the table's own identities;
the score against it is test_impingement_d2_validation.py.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from functools import cache

import pytest

from validation.cooling.schema import DATA_ROOT

TABLE = DATA_ROOT / "florschuetz_metzger_su1982" / "table_d2_heat_transfer.csv"
X_L = [-0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05]

Test = tuple[int, int, int, str, int, float, float, float]


@cache
def _tests() -> dict[Test, list[dict[str, float | None]]]:
    out: dict[Test, list[dict[str, float | None]]] = defaultdict(list)
    with TABLE.open() as f:
        for r in csv.DictReader(f):
            key = (
                int(r["xn_d"]),
                int(r["yn_d"]),
                int(r["z_d"]),
                r["pattern"],
                int(r["pdf_page"]),
                float(r["Re_mean_k"]),
                float(r["Nu_mean"]),
                float(r["mc_mj"]),
            )
            out[key].append(
                {
                    "x_L": float(r["x_L"]),
                    "Nu": float(r["Nu"]) if r["Nu"] else None,
                    "eta": float(r["eta"]) if r["eta"] else None,
                }
            )
    return dict(out)


def test_every_page_and_test_is_there() -> None:
    """12 pages (11 inline geometries and B(5,4,3)S), four tests each,
    every test across the 14 segments x/L -0.25 ... 1.05."""
    tests = _tests()
    assert len(tests) == 48
    pages = defaultdict(int)
    for key, rows in tests.items():
        pages[key[4]] += 1
        assert [r["x_L"] for r in rows] == X_L, key
    assert sorted(pages) == list(range(162, 174))
    assert set(pages.values()) == {4}


def test_the_mean_nusselt_number_closes() -> None:
    """NU-bar is the mean over the heat transfer surface opposite the array
    (0 < x/L < 1): the ten in-array segments. A misread tens or units digit
    in any of them breaks this."""
    for key, rows in _tests().items():
        inside = [r["Nu"] for r in rows if 0.0 < r["x_L"] < 1.0]
        assert len(inside) == 10 and None not in inside, key
        assert sum(inside) / 10 == pytest.approx(key[6], rel=0.006), key


def test_eta_is_one_upstream_of_the_array() -> None:
    """Upstream of row 1 only the initial crossflow sweeps the surface, so
    T_aw is its temperature: eta = 1. The zero-crossflow test has none."""
    for key, rows in _tests().items():
        mc = key[7]
        etas = [r["eta"] for r in rows]
        if mc == 0.0:
            assert all(e is None for e in etas), key
            continue
        assert [r["eta"] for r in rows if r["x_L"] < 0.0] == [1.0, 1.0, 1.0], key
        inside = [r["eta"] for r in rows if 0.0 < r["x_L"] < 1.0]
        # In-array eta spans 0.05-1.02: one printed overshoot past 1 (p163,
        # mc/mj 0.56, row 1), as CR-3936 reports where the crossflow
        # controls the surface.
        assert None not in inside and all(0.0 < e <= 1.02 for e in inside), key
