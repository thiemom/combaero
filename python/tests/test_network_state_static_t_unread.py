"""No residual reads ``NetworkMixtureState.T`` (#782).

The solver carries a node's stagnation T in ``Tt``; ``T`` is the computed
STATIC temperature, which results and diagnostics report. Until 0.8 ``T``
held the stagnation value and was read in 79 places that meant ``Tt`` (or
meant a static T and got the stagnation one, #824). Every residual reader
names ``Tt`` now, so the residuals do not depend on how ``T`` is computed.
This solves one network of every stress family with the attribute
instrumented and fails on any read during a residual evaluation.

A guard, not a proof: on 0.7's code these networks reach 45 of the 79
read sites (three seeds a family reach 51 at ten times the cost). The
conversion was checked by a full-suite run with every read logged.
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path
from typing import Any

import pytest

from combaero.network import NetworkSolver
from combaero.network.components import NetworkMixtureState

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))

from networks import FAMILIES  # noqa: E402


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_no_residual_reads_state_T(family: str, monkeypatch: pytest.MonkeyPatch) -> None:
    reads: list[str] = []
    in_residuals = [0]
    base = object.__getattribute__
    evaluate = NetworkSolver._residuals_and_jacobian

    def logged(self: NetworkMixtureState, name: str) -> object:
        if name == "T" and in_residuals[0]:
            frame = traceback.extract_stack(limit=2)[0]
            reads.append(f"{Path(frame.filename).name}:{frame.lineno}")
        return base(self, name)

    def counted(self: NetworkSolver, *args: Any, **kwargs: Any) -> Any:
        in_residuals[0] += 1
        try:
            return evaluate(self, *args, **kwargs)
        finally:
            in_residuals[0] -= 1

    monkeypatch.setattr(NetworkMixtureState, "__getattribute__", logged)
    monkeypatch.setattr(NetworkSolver, "_residuals_and_jacobian", counted)
    net, _ = FAMILIES[family](0)
    # Success is not the point: a failed solve has evaluated residuals too.
    NetworkSolver(net).solve(timeout=10.0)
    assert not reads, sorted(set(reads))
