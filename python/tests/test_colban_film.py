"""Colban, Thole and Bogard (2011) shaped-hole film effectiveness (#386).

The physics is pinned in tests/test_film_effectiveness.cpp. These cover the
binding, the transcription check against the Eq. (19) line the paper draws
(committed under validation/cooling/data/colban2011), and the runner that
scores the correlation on its own Figs. 6-7 data.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb  # noqa: E402
from validation.cooling import film_runner, scorecard  # noqa: E402
from validation.cooling import shaped_film_runner as sf  # noqa: E402
from validation.cooling.schema import load_dataset, load_points  # noqa: E402

# Gritsch configuration B, Figs. 6-7.
B = {"area_ratio": 3.5, "P_over_D": 6.0, "t_over_P": 0.49}


def _series(stem: str):
    return next(s for s in load_dataset() if s.label == f"colban2011/{stem}")


def test_breakout_value_is_the_coverage_ratio() -> None:
    assert cb.film_effectiveness_colban_2011(0.0, 1.0, **B) == pytest.approx(0.49, rel=1e-15)


def test_derivative_tuple_matches_the_value_and_finite_differences() -> None:
    eta, d_dM = cb.film_effectiveness_colban_2011_and_derivatives(15.0, 1.2, **B)
    assert eta == cb.film_effectiveness_colban_2011(15.0, 1.2, **B)
    h = 1e-6
    fd = (
        cb.film_effectiveness_colban_2011(15.0, 1.2 + h, **B)
        - cb.film_effectiveness_colban_2011(15.0, 1.2 - h, **B)
    ) / (2 * h)
    assert d_dM == pytest.approx(fd, rel=1e-6)
    assert d_dM > 0.0  # shaped holes do not lift off in the fitted range


def test_refuses_a_contracting_hole() -> None:
    with pytest.raises(ValueError, match="area_ratio"):
        cb.film_effectiveness_colban_2011(5.0, 1.0, 0.8, 6.0, 0.49)


@pytest.mark.parametrize(
    ("stem", "M"), [("fig6_eq19_line_M0p5", 0.5), ("fig7_eq19_line_M2p5", 2.5)]
)
def test_reproduces_the_eq19_line_the_paper_draws(stem: str, M: float) -> None:
    """Our Eq. (19) against the authors' own evaluation of it.

    The line was read raw off the figure, independent of this code. One
    pixel is 0.0008 (Fig. 6) and 0.0013 (Fig. 7); the bound is ~2-3 px.
    """
    pts = load_points(_series(stem))
    assert len(pts) >= 12
    for p in pts:
        ours = cb.film_effectiveness_colban_2011(p.x, M, **B)
        assert ours == pytest.approx(p.y, abs=0.003), f"x/D = {p.x}"


def test_fitting_data_scores_as_recorded() -> None:
    """Eq. (19) on the Gritsch points it was fitted on.

    Pinned so a change to the correlation or the data moves a number a
    reader can see. At M = 2.5 the correlation is LOW near the hole and HIGH
    downstream -- the same shape the paper's Fig. 7 shows -- so the RMS is
    an honest 12%, not a digitisation artefact.
    """
    out = {}
    for stem in ("fig6_eta_gritsch_M0p5", "fig7_eta_gritsch_M2p5"):
        recs = sf.run_series(_series(stem))
        errs = [r.rel_error for r in recs]
        assert all(e is not None for e in errs)
        assert not any(r.extrapolated for r in recs)
        out[stem] = (
            sum(errs) / len(errs),
            math.sqrt(sum(e * e for e in errs) / len(errs)),
        )
    bias6, rms6 = out["fig6_eta_gritsch_M0p5"]
    bias7, rms7 = out["fig7_eta_gritsch_M2p5"]
    assert rms6 == pytest.approx(0.036, abs=0.003)
    assert rms7 == pytest.approx(0.121, abs=0.005)
    assert abs(bias6) < 0.03


def test_the_line_series_are_not_scored() -> None:
    recs = sf.run_series(_series("fig6_eq19_line_M0p5"))
    assert all(r.predicted is None for r in recs)


def test_only_the_shaped_runner_owns_colban_series() -> None:
    for s in load_dataset():
        if s.source.name != "colban2011":
            continue
        assert sf.owns(s)
        assert not film_runner.owns(s), s.label


def test_basis_is_fidelity_on_the_correlations_own_data() -> None:
    s = _series("fig6_eta_gritsch_M0p5")
    assert scorecard.basis_of(s, "colban_2011_shaped") == "fidelity"
