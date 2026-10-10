"""Tier 2: the junction closure at finite Mach (#272).

Two sources were planned (the local docs/junction/tier2_reference_data.md;
the tracked record is validation/junction/README.md, "Tier 2"). Neither is the
closure's own paper (Mynard 2015), so both can only answer the ACCURACY question,
and the results are labelled CROSS-SOURCE (docs/VALIDATION_POLICY.md).

**Wang 2014** measures a Miller-type loss coefficient directly, at 45 deg, three
area ratios and M_3 from 0.09 to 0.60. It set no constant in the closure
(network_runner.SOURCE_ROLES), so it is the one independent accuracy check the
junction has. Wang states NO measurement uncertainty -- the full text has no
uncertainty, repeatability or error figure -- so there is no band to be "within";
the assertions are regression floors on the measured error, with margin, and a
falsification showing the score sees the closure.

**Perez-Garcia 2010** is a declared negative. Its K_hat (Eq 41) is a function of
the two branch Mach numbers only, and a loss reaches it only through the branch
density: across most of its envelope its own 4% U95 is worth tens to
thousands of units of Miller's K. Measured, a junction with NO loss at all
reproduces Table 1 within U95 in 57 of 72 cells and in every cell at
M3* <= 0.3. A "within U95" test -- the plan in the reference doc -- would pass
for any junction model, so it is not wired. The tests below pin why, so nobody
wires it later.
"""

from __future__ import annotations

import statistics

import pytest

from validation.junction.models import perez_garcia2010 as pg
from validation.junction.models.mpce_network import MPCENetwork
from validation.junction.network_runner import SOURCE_ROLES, _mach_bin, run_network

# ---------------------------------------------------------------------------
# Perez-Garcia 2010: why it is not a validation source for the closure
# ---------------------------------------------------------------------------


def test_perez_garcia_table_1_is_the_printed_table():
    """Read from the table image at 700 dpi (accepted manuscript). Seven values
    had been transcribed wrong; pinned so a revert is a visible diff."""
    printed = {
        ("C1", "K_hat_1"): (0.6871, 2.0283, 0.1296, 5.75),
        ("C1", "K_hat_2"): (0.7559, 2.0378, -0.0392, 4.35),
        ("C2", "K_hat_2"): (0.7504, 2.0223, -0.0957, 3.54),
        ("D1", "K_hat_1"): (0.7312, 2.0374, 0.0398, 4.64),
        ("D1", "K_hat_2"): (0.6718, 1.9543, 0.0242, 6.25),
        ("D2", "K_hat_2"): (0.7307, 1.9927, -0.1538, 3.99),
    }
    for key, (s, m, n_m1, u95) in printed.items():
        row = pg.TABLE_1[key]
        assert (row["s"], row["m"], row["n_m1"], row["U95pct"]) == (s, m, n_m1, u95), key


def test_k_hat_has_no_minus_one_in_its_denominator():
    """Eq 41: K_hat = (f(M3) - 1) / f(Mj). With the "- 1" the reference doc
    used to carry, the kinematic value would be ~4x the printed correlation;
    without it the two agree to 0.4% at C1, M3 = 0.5, q = 0.5 (each branch
    near M = 0.25)."""
    with_definition = pg.K_hat_from_mach(0.5, 0.25)
    printed = pg.K_hat("C1", "K_hat_1", 0.5, 0.5)
    assert with_definition == pytest.approx(printed, rel=0.01)


def _branch_fraction(flow_type: str, K_id: str, q: float) -> float:
    """Fig 2: q = G2/G3 for every flow type, so branch 1 carries 1 - q."""
    return 1.0 - q if K_id == "K_hat_1" else q


def _lossless_k_hat(M3: float, qj: float) -> float:
    """K_hat of an equal-area branch with Miller K = 0: p0_j = p0_3, Mach from
    continuity at the branch's share of the mass flux."""
    gamma = pg.GAMMA
    exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0))

    def flux(mach: float) -> float:
        return mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** exponent

    target = qj * flux(M3)
    lo, hi = 1e-12, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if flux(mid) < target else (lo, mid)
    return pg.K_hat_from_mach(M3, 0.5 * (lo + hi))


def test_a_junction_with_no_loss_passes_the_planned_criterion():
    """The falsification of the metric. Zero loss on every branch, scored the
    way the reference doc proposed (|K_hat / K_hat_corr - 1| < U95)."""
    inside = total = 0
    low_mach_inside = low_mach_total = 0
    for (flow_type, K_id), row in pg.TABLE_1.items():
        for M3 in (0.15, 0.3, 0.5, 0.7):
            for q in (0.25, 0.5, 0.75):
                qj = _branch_fraction(flow_type, K_id, q)
                err = _lossless_k_hat(M3, qj) / pg.K_hat(flow_type, K_id, q, M3) - 1.0
                ok = abs(err) * 100.0 <= row["U95pct"]
                inside += ok
                total += 1
                if M3 <= 0.3:
                    low_mach_inside += ok
                    low_mach_total += 1

    assert (inside, total) == (57, 72)
    assert low_mach_inside == low_mach_total


@pytest.mark.parametrize(
    ("M3", "qj", "at_least"),
    [(0.15, 0.5, 300.0), (0.3, 0.5, 20.0), (0.5, 0.5, 3.0), (0.7, 0.75, 0.4)],
)
def test_its_band_in_units_of_miller_k(M3, qj, at_least):
    """What 4% of K_hat is worth in Miller's K. Only the corner M3* >= 0.5
    with a branch carrying >= 3/4 of the flow gets near the model's own
    error against measurements (0.1-0.4); q = 0/1 would go further but is
    a dead-branch endpoint the harness cannot score."""
    assert 0.04 * pg.miller_K_per_unit_K_hat(M3, qj) > at_least


# ---------------------------------------------------------------------------
# Wang 2014: CROSS-SOURCE accuracy
# ---------------------------------------------------------------------------


def _wang_errors(**model_kwargs: float) -> list:
    return [
        r
        for r in run_network(MPCENetwork(**model_kwargs), topologies=("imposed_q",))
        if r.paper == "wang2014" and r.converged and r.K_extracted is not None
    ]


def _mae(records: list) -> float:
    return statistics.fmean(abs(r.K_extracted - r.K_measured) for r in records)


@pytest.fixture(scope="module")
def wang():
    return _wang_errors()


def test_wang_is_declared_cross_source():
    assert SOURCE_ROLES["wang2014"][0] == "x-source"


def test_wang_cross_source_accuracy(wang):
    """CROSS-SOURCE. Measured 2026-10-10: 116 of 200 points scored, MAE 0.108;
    by area ratio 0.072 / 0.089 / 0.181 at a = 1 / 1.56 / 2.44 (biases +0.06,
    0.00, -0.12). The 84 unscored are the dead-branch q = 0 and q = 1 curves
    (80, rejected by the direction verifier) and four a = 2.44, q = 0.8 points
    above M 0.49 -- coverage, tracked on #272, not error."""
    assert len(wang) >= 116
    assert _mae(wang) < 0.12
    by_area = {a: _mae([r for r in wang if r.psi == a]) for a in (1.0, 1.56, 2.44)}
    assert by_area[1.0] < 0.09
    assert by_area[1.56] < 0.11
    assert by_area[2.44] < 0.21


def test_wang_error_does_not_grow_with_mach(wang):
    """The compressible part of the claim. MAE by Mach band, measured:
    0.097 below 0.15, 0.117, 0.101, 0.113 above 0.45 -- flat. An
    incompressible reference head would add ~9% of K at M 0.6."""
    lowest = _mae([r for r in wang if _mach_bin(r.mach) == "M<0.15"])
    highest = _mae([r for r in wang if _mach_bin(r.mach) == "M>0.45"])
    assert highest < lowest + 0.04


def test_the_wang_score_sees_the_closure(wang):
    """Falsify the metric: switching off the joining energy transfer (alpha,
    which acts only when psi != 1) must move the a = 2.44 score, and must
    leave a = 1 untouched. Measured 0.181 -> 0.285 and 0.072 -> 0.072."""
    off = _wang_errors(joining_etransfer_alpha=0.0)
    on_244 = _mae([r for r in wang if r.psi == 2.44])
    off_244 = _mae([r for r in off if r.psi == 2.44])
    assert off_244 > on_244 + 0.05
    assert _mae([r for r in off if r.psi == 1.0]) == pytest.approx(
        _mae([r for r in wang if r.psi == 1.0]), abs=1e-6
    )
