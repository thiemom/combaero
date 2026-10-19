"""A junction with several inlets AND several outlets feeds each outlet its
own flow, not the whole throughput.

The solver handed every sink of a multi-source element each source's own
flow (flow_at_node at the source): right for a merge, which has one outlet,
but with two or more outlets each outlet chamber took its velocity head at the
junction's TOTAL flow. Converged roots were wrong (port flows off by up to 2x on stress seeds 11
and 13), and as a port's flow crossed zero its chamber
closure jumped by 10 kPa, which stalled 13 of 40 stress ideal_nport draws.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from combaero.network import NetworkSolver

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402


def _solve(seed: int) -> tuple:
    net, tags = FAMILIES["ideal_nport"](seed)
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=60.0)
    return net, s, r


@pytest.mark.parametrize("seed", [11, 13])
def test_each_outlet_chamber_carries_its_own_flow(seed: int) -> None:
    """Seeds 11 (2+4) and 13 (2+3), fully forward at the root: every outlet
    chamber's flow is its own outlet channel's, and they sum to the inlets'.

    Judged on the converged root: seed 13 draws loss_factor = 0, a lossless
    split at port Mach ~0.3 that creates flow work, and the second-law check
    refuses it (#598). Its port flows are this test's subject all the same."""
    net, s, r = _solve(seed)
    assert r["__converged__"], r["__message__"]
    if seed == 13:
        assert not r["__success__"] and "second law" in r["__message__"]
    s._residuals_and_jacobian(np.array(r["__x_solution__"]), compute_jacobian=False)
    outs = sorted(n for n in net.nodes if n.startswith("po"))
    assert len(outs) >= 2
    for n in outs:
        own = r[f"co{n[2:]}.m_dot"]
        assert own > 0.0
        assert net.nodes[n]._total_m_dot == pytest.approx(own, rel=1e-12)
    m_in = sum(r[f"{e}.m_dot"] for e in net.elements if e.startswith("ci"))
    assert sum(r[f"co{n[2:]}.m_dot"] for n in outs) == pytest.approx(m_in, rel=1e-9)


def test_the_stalled_draw_converges() -> None:
    """Seed 7 (2+2) stalled with outlet co0 at -2e-8 kg/s, on the edge of
    its chamber's 10 kPa jump."""
    _, _, r = _solve(7)
    assert r["__success__"], r["__message__"]


def test_an_outlet_chamber_flow_is_continuous_through_zero() -> None:
    """At seed 7's stalled iterate, the outlet chamber's flow on either side
    of its port flow's zero: about 0 both ways (it was the 0.069 kg/s
    throughput on one side)."""
    net, s, _ = _solve(7)
    i = list(s.unknown_names).index("co0.m_dot")
    x0 = np.array(s._last_x, dtype=float)
    totals = []
    for m in (-1e-8, 0.0, 1e-8):
        xx = x0.copy()
        xx[i] = m
        s._residuals_and_jacobian(xx, compute_jacobian=False)
        totals.append(net.nodes["po0"]._total_m_dot)
    assert max(abs(t) for t in totals) < 1e-6, totals
