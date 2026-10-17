"""A failed solve says when a declared-direction junction is why (stress tee3 triage).

MultiPortChamberElement and ConstantKTeeElement take a flow_direction and do
not reverse a port. The 8 stress tee3 draws that stalled all solve once the
junction is an IdealMomentumTeeElement, with one port reversed: they have no
root for the declared junction. The message said only "not making good
progress"; it now names the junction and the port.
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
def test_a_stall_on_a_declared_direction_names_the_junction_and_port(seed: int, kind: str) -> None:
    _, r = _solve(seed)
    assert not r["__success__"]
    msg = r["__message__"]
    assert f"{kind} 'jct'" in msg and "against its declared direction" in msg
    assert "IdealMomentumTeeElement" in msg


def test_a_converged_junction_gets_no_hint() -> None:
    s, r = _solve(0)
    assert r["__success__"], r["__message__"]
    assert s._declared_direction_hints(np.array(r["__x_solution__"])) == ""
