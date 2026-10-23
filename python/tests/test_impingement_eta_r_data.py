"""CR-3936 Appendix A (Florschuetz and Su, 1985): the crossflow temperature
influence factor eta_r per spanwise row, transcribed from the report's text
layer with each garbled glyph read off the page image (#467). The
extraction record is validation/cooling/extractions/florschuetz_su_1985_cr3936.md.
These tests check the transcription against the table's own identities; the
closure fitted to it comes once the extraction is signed off.
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from functools import cache

from validation.cooling.schema import DATA_ROOT

TABLE = DATA_ROOT / "florschuetz_su1985" / "appendix_a_row.csv"
N_ROWS = 10  # the array; runs after the row-9 heater burned out stop at row 7

Test = tuple[int, int, int, str, int, float, float]


@cache
def _tests() -> dict[Test, list[dict[str, float]]]:
    out: dict[Test, list[dict[str, float]]] = defaultdict(list)
    with TABLE.open() as f:
        for r in csv.DictReader(f):
            key = (
                int(r["xn_d"]),
                int(r["yn_d"]),
                int(r["z_d"]),
                r["pattern"],
                int(r["pdf_page"]),
                float(r["Re_mean_k"]),
                float(r["mc_mj"]),
            )
            out[key].append(
                {k: float(r[k]) for k in ("row", "closure", "Nu_r", "eta_r", "Re_j_k", "Gc_Gj")}
                | {"repeat_of_page": float(r["repeat_of_page"] or 0)}
            )
    return dict(out)


def test_every_page_and_block_is_there() -> None:
    """16 pages, three tests each. The segment heater opposite row 9 burned
    out before the last runs (Sec. 6.1): B(10,8,3)S on p156 and the
    Re_j-variation pages 157-160 report rows 1-7 only."""
    tests = _tests()
    assert len(tests) == 48
    pages = defaultdict(int)
    for (*_, page, _re, _mc), rows in tests.items():
        pages[page] += 1
        assert len(rows) == (7 if page >= 156 else 10), page
        assert [r["row"] for r in rows] == list(range(1, len(rows) + 1))
    assert sorted(pages) == list(range(145, 161))
    assert set(pages.values()) == {3}


def test_the_repeated_runs_are_marked() -> None:
    """Three Re_j-variation pages reprint a base test digit for digit; a
    fit must count each once (one test in two folds leaks)."""
    tests = _tests()
    repeats = {
        (k[4], int(v[0]["repeat_of_page"])) for k, v in tests.items() if v[0]["repeat_of_page"]
    }
    assert repeats == {(158, 148), (159, 149), (160, 153)}
    for (*_, re_mean, mc), rows in tests.items():
        src = int(rows[0]["repeat_of_page"])
        if not src:
            continue
        base = next(v for k, v in tests.items() if k[4] == src and k[6] == mc and k[5] == re_mean)
        assert [r["eta_r"] for r in rows] == [r["eta_r"] for r in base[: len(rows)]]


def test_the_transcription_closes_on_itself() -> None:
    """Each Gc/Gj follows from the initial crossflow and the jets upstream:

        Gc_i/Gj_i = (pi/4)/((yn/d)(z/d)) (S N (mc/mj) Re_mean + sum_{k<i} Re_k) / Re_i

    S is the test's closure (metered initial crossflow against the inferred
    jets). CR-3630 prints it as MC, 0.91-1.02; CR-3936 does not, so it is
    the median over the test's rows, and a misread digit shows as a row
    breaking from it. This check found the lost B(5,8,2) header (a factor
    of 2: z/d 2 filed as 1)."""
    for (_xn, yn, z, _pat, page, re_mean, mc), rows in _tests().items():
        r = [row["Re_j_k"] for row in rows]
        scale = []
        for i, row in enumerate(rows):
            if row["Gc_Gj"] >= 0.3:
                c = row["Gc_Gj"] * r[i] * yn * z / (math.pi / 4) - sum(r[:i])
                scale.append(c / (N_ROWS * mc * re_mean))
        if not scale:  # crossflow below 0.3 on every row: nothing to check
            assert rows[0]["closure"] == 1.0, page
            continue
        s = statistics.median_high(scale)  # as the parser stored it
        assert abs(s - rows[0]["closure"]) < 6e-4, page  # stored to 3 decimals
        # p150's B(10,4,2)I at mc/mj 0.20 closes at 1.25 on every row: as
        # printed (page image checked). Every other test is within 0.8-1.2.
        assert 0.8 < s < (1.3 if (page, mc) == (150, 0.20) else 1.2), (page, mc, s)
        for i, row in enumerate(rows):
            if row["Gc_Gj"] < 0.3:
                continue
            want = (math.pi / 4) / (yn * z) * (s * N_ROWS * mc * re_mean + sum(r[:i])) / r[i]
            # Kept as printed (page images checked): p155 B(5,4,3)S at mc/mj
            # 1.00, row 8, 0.86 against 0.92; p150's row 1, 0.32 against
            # 0.34. Everything else within 5%.
            tol = 0.07 if (page, mc, i + 1) in {(155, 1.0, 8), (150, 0.2, 1)} else 0.05
            assert abs(row["Gc_Gj"] / want - 1) <= tol, (page, mc, i + 1)


def test_eta_r_spans_what_the_authors_describe() -> None:
    """eta_r rises from the jet-dominated end toward 1 and past it once the
    crossflow controls the surface (Sec. 6.2: four geometries exceed 1,
    only at Gc/Gj of order one or more)."""
    rows = [r for v in _tests().values() for r in v]
    eta = [r["eta_r"] for r in rows]
    assert min(eta) > 0.3 and max(eta) < 1.5
    over = [r for r in rows if r["eta_r"] > 1.0]
    assert over and min(r["Gc_Gj"] for r in over) > 0.5
