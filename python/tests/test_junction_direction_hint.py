"""A failed solve says when a declared-direction junction is why (stress tee3 triage).

MultiPortChamberElement and ConstantKTeeElement take a flow_direction and do
not reverse a port. The 8 stress tee3 draws that stalled all solve once the
junction is an IdealMomentumTeeElement, with one port reversed: they have no
root for the declared junction. The message said only "not making good
progress"; it now names the junction and the port.

Where a stall ends depends on wall-clock budgets, so the hint is tested on a
state built to hold a port reversed, and the message wiring separately.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from combaero.network import NetworkSolver

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import tee3  # noqa: E402


def _solve(seed: int) -> tuple[NetworkSolver, dict]:
    net, _ = tee3(seed)
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve(timeout=30.0)


@pytest.mark.parametrize(
    ("seed", "kind"), [(31, "MultiPortChamberElement"), (23, "ConstantKTeeElement")]
)
def test_a_port_against_its_declared_direction_is_named(seed: int, kind: str) -> None:
    s, r = _solve(seed)
    x = np.array(r["__x_solution__"], dtype=float)
    jct = s.network.elements["jct"]
    # Reverse the declared straight port: a merge's straight INLET, a branch's
    # straight OUTLET; arm1 (merge) / arm1 (branch) is that port's outer element.
    oid = jct._port_element_ids[jct.port_nodes.index("ps")]
    i = jct._port_outer_mdot_idx[oid]
    x[i] = -abs(x[i]) - 0.01
    hint = s._declared_direction_hints(x)
    assert f"{kind} 'jct'" in hint and "['ps']" in hint
    assert "IdealMomentumTeeElement" in hint


def test_a_converged_junction_gets_no_hint() -> None:
    s, r = _solve(0)
    assert r["__success__"], r["__message__"]
    assert s._declared_direction_hints(np.array(r["__x_solution__"])) == ""


def test_a_failed_solve_carries_the_hint_in_its_message(monkeypatch) -> None:
    monkeypatch.setattr(NetworkSolver, "_declared_direction_hints", lambda self, x: "HINT-XYZ")
    _, r = _solve(23)
    if r["__success__"]:
        pytest.skip("seed 23 converged here; the wiring applies to failed solves only")
    assert "HINT-XYZ" in r["__message__"]
