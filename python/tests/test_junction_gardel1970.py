"""Gardel & Rechsteiner 1970 cross-source score: plumbing and regression floors.

The score itself is ACCURACY on another lab's rig (validation/junction/
gardel1970_score.py); the source states no band, so these tests pin how the
models are evaluated and where the numbers sit, not a pass mark. A miss in the
plumbing tests is our bug; a move in the floors needs an explanation.
"""

from __future__ import annotations

import math

import pytest

from validation.junction import gardel1970_score as g
from validation.junction.models import bassett2001 as bassett


@pytest.mark.parametrize("delta,phi", [(90.0, 1.0), (90.0, 0.69), (135.0, 1.0)])
@pytest.mark.parametrize("q", [0.2, 0.5, 0.8])
def test_dividing_mode_is_bassett_type_3(delta: float, phi: float, q: float) -> None:
    # Gardel's delta is to the UPSTREAM main, so the branch deflection is pi - delta.
    k = g.model_losses("bassett", delta, phi, "1", q)
    th = math.radians(180.0 - delta)
    assert k[g.BETA] == pytest.approx(bassett.K6(q, 1.0 / phi, th), abs=1e-12)
    assert k[g.GAMMA] == pytest.approx(bassett.K5(1.0 - q), abs=1e-12)


@pytest.mark.parametrize("q", [0.2, 0.5, 0.8])
def test_ideal_kernel_is_the_vector_borda_loss(q: float) -> None:
    # Mode 1 at 135 deg, phi 0.69: |v_j - v_alpha|^2 on V = 1.
    phi, d = 0.69, math.radians(45.0)
    # 1e-6: the kernel blends supplier and receiver roles over 1e-4 of the
    # throughput, which moves a row by O(1e-8) (include/ideal_junction.h).
    k = g.model_losses("ideal", 135.0, phi, "1", q)
    r = q / phi
    assert k[g.BETA] == pytest.approx(1.0 + r * r - 2.0 * r * math.cos(d), abs=1e-6)
    assert k[g.GAMMA] == pytest.approx(q * q, abs=1e-6)


def test_measured_losses_at_the_panel_ends_match_the_tables() -> None:
    # Table F/G/E, GR 90/1/0: the sign mapping must turn each into a loss.
    assert g.loss_from_h(-1.03, "h_ba", g.MODES["1"][1](1.0)) == pytest.approx(1.03)
    assert g.loss_from_h(-1.01, "h_yb", g.MODES["2"][1](-1.0)) == pytest.approx(1.01)
    assert g.loss_from_h(-0.02, "h_ya", g.MODES["2"][1](0.0)) == pytest.approx(0.02)
    # Modes 3 and 4 run through the branch the other way.
    assert g.loss_from_h(1.0, "h_ba", g.MODES["3"][1](0.0)) == pytest.approx(1.0)
    assert g.loss_from_h(-1.0, "h_ba", g.MODES["4"][1](0.0)) == pytest.approx(1.0)


def test_pooled_scores_sit_where_they_were_measured() -> None:
    cells, taps = g.score()
    stats = g.pooled(cells)
    assert stats["ideal"][0] == 1426
    assert taps == 18
    # 2026-10-12: ideal 0.550 (+0.535), Bassett 0.214 (+0.115), shipped 0.219 (+0.113).
    assert stats["ideal"][1] == pytest.approx(0.550, abs=0.01)
    assert stats["bassett"][1] == pytest.approx(0.214, abs=0.01)
    assert stats["shipped"][1] == pytest.approx(0.219, abs=0.01)


def test_the_score_moves_when_a_model_does(monkeypatch: pytest.MonkeyPatch) -> None:
    # Falsify the metric: a constant 0.3 on Bassett's losses must show up in his bias.
    base = g.pooled(g.score()[0])["bassett"]
    real = g.model_losses

    def shifted(model: str, *a: object) -> dict[int, float]:
        out = real(model, *a)
        return {j: v + 0.3 for j, v in out.items()} if model == "bassett" else out

    monkeypatch.setattr(g, "model_losses", shifted)
    moved = g.pooled(g.score()[0])["bassett"]
    assert moved[2] - base[2] == pytest.approx(0.3, abs=1e-9)
