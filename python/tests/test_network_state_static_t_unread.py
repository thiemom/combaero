"""No residual reads ``NetworkMixtureState.T`` (#782).

The solver carries a node's stagnation T in ``Tt``; ``T`` is filled with
the same value and was read in 79 places that meant ``Tt`` (or meant a
static T and got the stagnation one -- the physics follow-up). Every reader
now names ``Tt``, so ``T`` can become a computed static property without
moving a result. This solves one network of every stress family with the
attribute instrumented and fails on any read.

A guard, not a proof: on main's code these networks reach 45 of the 79
read sites (three seeds a family reach 51 at ten times the cost). The
conversion was checked by a full-suite run with every read logged.
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

import pytest

from combaero.network import NetworkSolver
from combaero.network.components import NetworkMixtureState

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))

from networks import FAMILIES  # noqa: E402


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_no_residual_reads_state_T(family: str, monkeypatch: pytest.MonkeyPatch) -> None:
    reads: list[str] = []
    base = object.__getattribute__

    def logged(self: NetworkMixtureState, name: str) -> object:
        if name == "T":
            frame = traceback.extract_stack(limit=2)[0]
            reads.append(f"{Path(frame.filename).name}:{frame.lineno}")
        return base(self, name)

    monkeypatch.setattr(NetworkMixtureState, "__getattribute__", logged)
    net, _ = FAMILIES[family](0)
    # Success is not the point: a failed solve has evaluated residuals too.
    NetworkSolver(net).solve(timeout=10.0)
    assert not reads, sorted(set(reads))
