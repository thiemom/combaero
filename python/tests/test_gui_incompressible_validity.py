"""The incompressible validity report reaches the GUI in canvas terms (#588)."""

from __future__ import annotations

import asyncio

import pytest
from test_gui_tee_port_step import _schema

from combaero.network import compressibility
from gui.backend.main import solve
from gui.backend.schemas import NetworkGraphSchema


def _solve():
    schema = NetworkGraphSchema.model_validate(_schema("mynard", 0.05))
    schema.solver_settings.global_regime = "incompressible"
    return asyncio.run(solve(schema))


def test_each_incompressible_element_carries_its_verdict() -> None:
    res = _solve()
    assert res.success
    row = res.element_results["ch_main"].model_dump()
    assert row["incompressible_validity"] in ("ok", "warn")
    assert row["incompressible_mach_limit"] == pytest.approx(0.2)
    assert row["incompressible_mach"] >= 0.0


def test_a_flag_is_attributed_to_its_canvas_node(monkeypatch) -> None:
    """Force a warning on a drawn channel: it is listed against that node,
    and the solve still succeeds -- a report, not a rejection."""
    original = compressibility.validity_verdicts

    def with_a_warning(network, diag, limit, node_mach=None):
        out = original(network, diag, limit, node_mach)
        out["ch_main"] = {"verdict": "warn", "mach": 0.31, "limit": limit}
        return out

    monkeypatch.setattr(compressibility, "validity_verdicts", with_a_warning)
    res = _solve()
    assert res.success
    flags = [(f.element_id, f.owner_id, f.mach, f.limit) for f in res.incompressible_warnings]
    assert flags == [("ch_main", "ch_main", 0.31, 0.2)]
