"""Every backend path configures the solver the same way.

The /solve route turned on the Newton-homotopy last retry (#642) with no cap
of its own, while NetworkRunner solved the same saved drawing without it: a
script could fail where the GUI converged. Both now build the solver through
``gui.backend.solver_setup.make_solver``.
"""

from __future__ import annotations

import asyncio
import copy
from typing import Any

import pandas as pd
import pytest
from test_network_runner import _SIMPLE

from combaero.network import NetworkSolver
from gui.backend.main import export_results, solve
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema

_MACH_LIMIT = 0.17


@pytest.fixture
def seen(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    """The solver settings in force at every top-level solve() call."""
    calls: list[dict[str, Any]] = []
    original = NetworkSolver.solve

    def spy(self: NetworkSolver, *args: Any, **kwargs: Any) -> Any:
        calls.append(
            {
                "homotopy_retry": self.homotopy_retry,
                "homotopy_budget": self.homotopy_budget,
                "bc_ramp_retry": self.bc_ramp_retry,
                "incompressible_mach_limit": self.incompressible_mach_limit,
            }
        )
        return original(self, *args, **kwargs)

    monkeypatch.setattr(NetworkSolver, "solve", spy)
    return calls


def _schema() -> dict[str, Any]:
    schema = copy.deepcopy(_SIMPLE)
    schema["solver_settings"] = {"incompressible_mach_limit": _MACH_LIMIT}
    return schema


_EXPECTED = {
    "homotopy_retry": True,
    "homotopy_budget": None,
    "bc_ramp_retry": True,
    "incompressible_mach_limit": _MACH_LIMIT,
}


def test_the_gui_route_takes_the_homotopy_retry(seen: list[dict[str, Any]]) -> None:
    asyncio.run(solve(NetworkGraphSchema(**_schema())))
    assert seen and all(call == _EXPECTED for call in seen)


def test_the_export_route_takes_the_homotopy_retry(seen: list[dict[str, Any]]) -> None:
    asyncio.run(export_results(NetworkGraphSchema(**_schema())))
    assert seen and all(call == _EXPECTED for call in seen)


def test_network_runner_solve_matches_the_gui(seen: list[dict[str, Any]]) -> None:
    NetworkRunner.from_dict(_schema()).solve()
    assert seen and all(call == _EXPECTED for call in seen)


def test_network_runner_sweep_matches_the_gui(seen: list[dict[str, Any]]) -> None:
    NetworkRunner.from_dict(_schema()).sweep(pd.DataFrame({"inlet.m_dot": [0.5, 1.0]}))
    assert len(seen) >= 2 and all(call == _EXPECTED for call in seen)


@pytest.mark.parametrize("flag", [True, False])
def test_every_backend_path_passes_the_stall_handover(
    monkeypatch: pytest.MonkeyPatch, flag: bool
) -> None:
    """The drawing's stall_handover (#864) reaches NetworkSolver.solve from
    the /solve route, the /export route and NetworkRunner alike."""
    passed: list[Any] = []
    depth = [0]
    original = NetworkSolver.solve

    def spy(self: NetworkSolver, *args: Any, **kwargs: Any) -> Any:
        if depth[0] == 0:
            passed.append(kwargs.get("stall_handover", "missing"))
        depth[0] += 1
        try:
            return original(self, *args, **kwargs)
        finally:
            depth[0] -= 1

    monkeypatch.setattr(NetworkSolver, "solve", spy)
    schema = copy.deepcopy(_SIMPLE)
    schema["solver_settings"] = {"stall_handover": flag}
    asyncio.run(solve(NetworkGraphSchema(**schema)))
    asyncio.run(export_results(NetworkGraphSchema(**schema)))
    NetworkRunner.from_dict(schema).solve()
    assert passed == [flag, flag, flag]
