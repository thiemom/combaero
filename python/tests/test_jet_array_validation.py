"""Score florschuetz_1981_inline() against Figure 6 of its own source.

jet_array_runner.py is deliberately separate from runner.py: Nu/Nu1 needs
no Reynolds-number bisection (see the module docstring), so it has none of
runner.py's machinery and should not share its dispatch. See
validation/cooling/extractions/han_impingement.md and issue #337.
"""

from __future__ import annotations

import math

import pytest

from validation.cooling.jet_array_runner import Record, run_all, run_series
from validation.cooling.schema import load_dataset


@pytest.fixture(scope="module")
def dataset():
    return [s for s in load_dataset() if s.source.name == "florschuetz1981"]


@pytest.fixture(scope="module")
def records(dataset) -> list[Record]:
    return run_all(dataset)


def _series(dataset, stem: str):
    for s in dataset:
        if s.path.stem == stem:
            return s
    raise AssertionError(f"no series {stem}")


def test_all_27_series_are_measured_and_scored(dataset, records) -> None:
    """Every Figure 6 panel/shape combination loads and scores.

    27 = 9 panels x 3 z/d symbol classes. A count drift here means a panel
    was dropped from metadata.yaml or a file went missing, not a scoring
    change.
    """
    assert len(dataset) == 27
    assert all(s.kind == "measured" for s in dataset)
    unscored = [r for r in records if r.predicted is None]
    assert not unscored, [r.reason for r in unscored]


def test_figure_6_reproduces_its_own_source(records) -> None:
    """Pooled across all 242 points: bias and RMS against the real Nu ratio.

    Nu_j/Nu1_j = jet_array_impingement_nu(..., Gc_Gj) / jet_array_impingement_nu(
    ..., Gc_Gj=0), both calls hitting the real implementation -- the formula
    is never re-derived here. Bands are set from what this scoring pass
    actually measured, not fitted to make the test pass: 9.6% RMS sits
    close to the fit's own stated 5.6% standard error, which is expected
    since standard_error is the FIT's training residual and this is raw
    digitised scatter including the one outlier panel below.
    """
    errs = [r.rel_error for r in records if r.rel_error is not None]
    assert len(errs) == 242
    bias = sum(errs) / len(errs)
    rms = math.sqrt(sum(e * e for e in errs) / len(errs))
    assert abs(bias) < 0.06, f"bias {bias:.1%}"
    assert rms < 0.12, f"RMS {rms:.1%}"


def test_b_5_4_circle_is_the_documented_validity_corner_outlier(dataset) -> None:
    """B(5,4)I, z/d=1 sits at the correlation's own validity corner.

    xn/d=5, yn/d=4, z/d=1 are simultaneously at florschuetz_1981_inline()'s
    lower bound on all three. The fit underpredicts measured crossflow
    degradation there badly (documented in this series' own cross_check).
    Pinned as a known, understood discrepancy rather than a silent
    exception -- if this ever tightens up, the cross_check text and this
    bound both need revisiting together, not just the number here.
    """
    recs = run_series(_series(dataset, "rc_00_circles_zd_1"))
    errs = [r.rel_error for r in recs if r.rel_error is not None]
    assert len(errs) == 9
    bias = sum(errs) / len(errs)
    assert bias > 0.25, f"outlier bias {bias:.1%} no longer stands out as documented"


def test_wrong_geometry_changes_the_score(dataset) -> None:
    """Falsification: swapping z/d must move the RMS, or nothing is measured.

    Confirms jet_array_runner actually reads geometry off the series rather
    than defaulting silently -- a runner that ignored xn_d/yn_d/z_d would
    pass every test above by accident.
    """
    import dataclasses

    real = _series(dataset, "rc_00_circles_zd_1")
    wrong = dataclasses.replace(real, geometry={"xn_d": 15.0, "yn_d": 8.0, "z_d": 3.0})

    def rmse(recs):
        errs = [r.rel_error for r in recs if r.rel_error is not None]
        return math.sqrt(sum(e * e for e in errs) / len(errs))

    assert rmse(run_series(wrong)) != pytest.approx(rmse(run_series(real)), rel=1e-6)


def test_gc_gj_zero_scores_a_perfect_ratio(dataset) -> None:
    """At Gc/Gj = 0 the ratio must be exactly 1 by the formula's own construction.

    Not a property of the digitised data (no real point sits at exactly
    zero) -- a direct check on _nu_ratio's own algebra, independent of any
    figure.
    """

    import combaero as cb
    from validation.cooling.jet_array_runner import _nu_ratio

    real = _series(dataset, "rc_00_circles_zd_1")
    jet_set = cb.florschuetz_1981_inline()
    ratio, _ = _nu_ratio(
        jet_set, 0.0, real.geometry["xn_d"], real.geometry["yn_d"], real.geometry["z_d"]
    )
    assert ratio == pytest.approx(1.0, abs=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
