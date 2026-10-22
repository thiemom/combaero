"""An impingement array with an initial crossflow, against the authors' own
flow data: Florschuetz, Metzger and Su (1982), NASA CR-3630, Table D.1 (#467).

42 tests: 11 inline geometries (xn/d 5, 10; yn/d 4, 8; z/d 1, 2, 3), each at
initial-to-jet flow ratios mc/mj of 0 to about 1, ten rows each. Per row the
table gives P/P0 (channel static over jet plenum), Gj/Gj_mean and Gc/Gj. The
jet flows are the authors' own, inferred from the measured P with a constant
plate Cd (their Eq. 4.1): that is the quantity this network computes, so the
comparison is like for like. Fidelity: the source's own data
(docs/VALIDATION_POLICY.md).

The rig is reproduced as built: jet plenum at the measured P0 and T0; the
total jet flow from Re_j; the initial crossflow a fixed flow delivered along
the gap's upstream end (CR-3630 ran it 16 to 95 hydraulic diameters before
the array, so it arrives with its momentum, which the lead-in carries); the
plate's own mean Cd from Table 3.1; channel width 18.29 cm (18 holes at
yn/d 4, d 0.254 cm).

What the network leaves out, and where it shows: the authors' own analysis
says the discharge coefficient falls with crossflow and wall shear matters
at z/d = 1 with large crossflow, and that no adequate friction model exists
for these flows. Those cases are reported, with a looser bound, not dropped.
"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from functools import cache

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    ImpingementArray,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)
from validation.cooling.schema import DATA_ROOT

TABLE = DATA_ROOT / "florschuetz_metzger_su1982" / "table_d1_flow.csv"
N_ROWS = 10
SPAN = 0.1829  # channel width [m]: 18 holes at yn/d 4, d 0.254 cm
# Table 3.1: hole diameter by streamwise pitch, mean plate Cd by (xn/d, yn/d).
D_JET = {5: 0.00254, 10: 0.00127}
PLATE_CD = {(5, 4): 0.85, (5, 8): 0.80, (10, 4): 0.76, (10, 8): 0.76}
Y_AIR = cb.mole_to_mass(cb.species.dry_air())
X_AIR = cb.species.dry_air()

# Over 20 s on one worker: runs in CI's parallel slow job.
pytestmark = pytest.mark.slow

Key = tuple[int, int, int, float]


@cache
def _table() -> dict[Key, list[dict[str, float]]]:
    out: dict[Key, list[dict[str, float]]] = defaultdict(list)
    with TABLE.open() as fh:
        for r in csv.DictReader(fh):
            key = (int(r["xn_d"]), int(r["yn_d"]), int(r["z_d"]), float(r["mc_mj"]))
            out[key].append({k: (float(v) if v else math.nan) for k, v in r.items()})
    return dict(out)


@cache
def _solve(key: Key, with_bypass: bool = True) -> tuple[list[float], list[float], list[float]]:
    """Per row (Gj/Gj_mean, Gc/Gj, P/P0). ``with_bypass=False`` drops the
    initial crossflow, keeping the jets: the falsification."""
    xn, yn, z, mc = key
    rows = _table()[key]
    d = D_JET[xn]
    T0, P0, Re_j = rows[0]["T0_K"], rows[0]["P0_kPa"] * 1.0e3, rows[0]["Re_j"]
    mu = cb.complete_state(T0, P0, X_AIR).transport.mu
    n_holes = round(SPAN / (yn * d))
    m_j = Re_j * mu * N_ROWS * n_holes * math.pi * d / 4.0
    m_c = mc * m_j if with_bypass else 0.0

    net = FlowNetwork()
    supply = PressureBoundary("supply")
    supply.Pt, supply.Tt, supply.Y = P0, T0, Y_AIR
    net.add_node(supply)
    net.add_node(MassFlowBoundary("exit", m_dot=m_j + m_c, Tt=T0, Y=Y_AIR))
    arr = ImpingementArray(
        "a",
        n_rows=N_ROWS,
        d_jet=d,
        xn_d=xn,
        yn_d=yn,
        z_d=z,
        span=SPAN,
        plate_thickness=d,  # Table 3.1: b = d
        Cd=PLATE_CD[(xn, yn)],
    )
    bypass = None
    if m_c > 0.0:
        # The metered initial crossflow, straight into the lead-in. The rig's
        # 24 cm development run only sets the feed pressure, which floats
        # here: at a fixed bypass flow it cannot move the split (see
        # test_impingement_bypass.py::test_the_feed_node_does_not_matter).
        net.add_node(MassFlowBoundary("feed", m_dot=m_c, Tt=T0, Y=Y_AIR))
        bypass = "feed"
    arr.add_to(net, "supply", "exit", bypass=bypass)
    # The default start is warm for impingement networks (#706).
    res = NetworkSolver(net).solve()
    assert res["__success__"], f"{key}: {res.get('__message__')}"
    diag = res["__element_diag__"]
    m = [diag[arr.plate_id(i)]["m_dot"] for i in range(1, N_ROWS + 1)]
    mean = sum(m) / N_ROWS
    return (
        [x / mean for x in m],
        [diag[arr.plate_id(i)]["Gc_Gj"] for i in range(1, N_ROWS + 1)],
        [res[f"{arr.node_id(i)}.P"] / P0 for i in range(1, N_ROWS + 1)],
    )


def _gj_err(key: Key, with_bypass: bool = True) -> float:
    """Worst row's |Gj/Gj_mean - measured| for one test."""
    gj, _, _ = _solve(key, with_bypass)
    return max(abs(g - r["Gj_Gjbar"]) for g, r in zip(gj, _table()[key], strict=True))


def _closed_form(xn: int, yn: int, z: int, mc: float) -> list[float]:
    """The authors' own model, CR-3630 Eq. 4.10 (constant Cd, no friction)."""
    B = N_ROWS * math.sqrt(2.0) * PLATE_CD[(xn, yn)] * (math.pi / 4.0) / (yn * z)
    return [
        B * ((1 + mc) * math.cosh(B * x) - mc * math.cosh(B * (1 - x))) / math.sinh(B)
        for x in ((i + 0.5) / N_ROWS for i in range(N_ROWS))
    ]


def test_the_transcription_closes_on_itself() -> None:
    """Table D.1 is parsed from the report's text layer. Two of its own
    identities catch a misread digit: the jet flows average to 1, and each
    Gc/Gj follows from the initial crossflow and the jets upstream of it,

        Gc_i/Gj_i = (pi/4)/((yn/d)(z/d)) (N (mc/mj)/MC + sum_{k<i} g_k) / g_i

    MC is the table's mass closure (inferred jet flows over the metered
    total, 0.91 to 1.02): the authors' Gc counts the metered initial
    crossflow against the inferred jets. Without it B(5,8,1) at mc/mj 0.52,
    MC 0.91, misses by 11% -- the printed 1.18, confirmed on the page."""
    table = _table()
    assert len(table) == 42
    for (xn, yn, z, mc), rows in table.items():
        g = [r["Gj_Gjbar"] for r in rows]
        assert len(g) == N_ROWS
        assert sum(g) / N_ROWS == pytest.approx(1.0, abs=0.003), (xn, yn, z, mc)
        closure = rows[0]["mass_closure"]
        for i, r in enumerate(rows):
            # 0.30: below that the printed two decimals dominate the ratio;
            # a few rows have Gj 0.0-0.02, where nothing can be checked.
            if r["Gc_Gj"] < 0.30 or g[i] < 0.10:
                continue
            want = (math.pi / 4) / (yn * z) * (N_ROWS * mc / closure + sum(g[:i])) / g[i]
            # Measured worst 4.4%, two-decimal rounding.
            assert r["Gc_Gj"] == pytest.approx(want, rel=0.06), (xn, yn, z, mc, i + 1)


def test_the_flow_distribution_matches_the_authors_data() -> None:
    """Row-by-row Gj/Gj_mean over all 42 tests.

    Measured 2026-10-22, worst row per test: z/d >= 2 median 0.021, max
    0.100 (B(5,4,2) at mc/mj 1.01); z/d = 1 up to 0.265 (B(5,4,1) at mc/mj
    0.99), where the authors' own constant-Cd closed form misses by 1.22
    and the network's compressibility and wall shear do better. P/P0 mean
    absolute error 0.007.
    """
    wide, narrow, p_err = [], [], []
    for key, rows in sorted(_table().items()):
        err = _gj_err(key)
        (narrow if key[2] == 1 else wide).append(err)
        _, _, p = _solve(key)
        p_err += [abs(a - r["P_P0"]) for a, r in zip(p, rows, strict=True)]
    print(
        f"\nGj/Gj_mean worst row: z/d>=2 median {statistics.median(wide):.3f} max {max(wide):.3f};"
        f" z/d=1 max {max(narrow):.3f}; P/P0 mean |err| {statistics.mean(p_err):.4f}"
    )
    assert statistics.median(wide) < 0.03
    assert max(wide) < 0.11
    assert max(narrow) < 0.30
    assert statistics.mean(p_err) < 0.010


def test_row_ones_crossflow_is_the_bypass() -> None:
    """Gc/Gj at row 1 is the initial crossflow alone: where it is large the
    table pins the hand-off directly. Rows with z/d >= 2 and measured Gc/Gj
    above 0.3; measured worst 2026-10-22: 11%, which includes the table's own
    mass closure (down to 0.91) between the metered bypass and the jets."""
    errs = []
    for key, rows in _table().items():
        if key[2] == 1 or key[3] == 0.0 or rows[0]["Gc_Gj"] < 0.3:
            continue
        _, gc, _ = _solve(key)
        errs.append(abs(gc[0] / rows[0]["Gc_Gj"] - 1.0))
    print(f"\nrow-1 Gc/Gj: n {len(errs)}, worst {max(errs):.1%}")
    assert len(errs) >= 8
    assert max(errs) < 0.15


def test_the_network_discretises_the_authors_model() -> None:
    """At xn/d 5, z/d >= 2 -- low pressure ratio, little wall shear -- the
    network is the authors' Eq. 4.10 discretised: within 0.045 of it on
    every row. At xn/d 10 (P0 up to 134 kPa through 1.27 mm holes) the
    compressible jets make the supply measurably more uniform than the
    incompressible closed form; that is physics the closed form omits."""
    worst = 0.0
    for (xn, yn, z, mc), _rows in _table().items():
        if xn != 5 or z == 1:
            continue
        gj, _, _ = _solve((xn, yn, z, mc))
        cf = _closed_form(xn, yn, z, mc)
        worst = max(worst, max(abs(a - b) for a, b in zip(gj, cf, strict=True)))
    print(f"\nnetwork vs Eq. 4.10, xn/d 5, z/d >= 2: worst row {worst:.3f}")
    assert worst < 0.045


def test_dropping_the_bypass_scores_worse() -> None:
    """Falsify the metric: the same jets with no initial crossflow must miss
    the tests that had one. Every z/d >= 2 test at mc/mj near 1."""
    for key in _table():
        if key[2] == 1 or key[3] < 0.95:
            continue
        assert _gj_err(key, with_bypass=False) > 2.0 * _gj_err(key), key
