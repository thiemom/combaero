"""The second-law verdicts and the ledger reach the GUI in canvas terms (#545).

A verdict on an element the user never drew -- a tee port's sharp step, a
liner's segment -- is attributed to the node on the canvas that owns it.
"""

from __future__ import annotations

import asyncio

import pytest
from test_gui_tee_port_step import TEE, _schema

from combaero.network import second_law
from gui.backend.main import solve
from gui.backend.schemas import NetworkGraphSchema
from gui.backend.second_law_view import owner_of


@pytest.mark.parametrize(
    "element_id,owner",
    [
        ("ch_main", "ch_main"),  # drawn
        (f"__tee_step__{TEE}_straight", TEE),  # a tee's port step
        ("ln__s0", "ln"),  # a composite's sub-element
        ("orphan__x", "orphan__x"),  # no canvas owner: itself
        ("__auto_link__ch_main__ln", "ln"),  # a drawn node-to-node edge
    ],
)
def test_each_verdict_belongs_to_a_node_on_the_canvas(element_id: str, owner: str) -> None:
    gui_ids = {"ch_main", TEE, "ln"}
    steps = {f"__tee_step__{TEE}_straight": (TEE, "straight", 0.25)}
    assert owner_of(element_id, gui_ids, steps) == owner


def _solve(model: str = "mynard"):
    schema = NetworkGraphSchema.model_validate(_schema(model, 0.05))
    return asyncio.run(solve(schema))


def test_a_working_network_reports_its_balance_and_no_failures() -> None:
    res = _solve()
    assert res.success
    assert res.second_law_failures == []
    assert abs(res.ledger["energy_closure_rel"]) < 1e-6
    assert abs(res.ledger["mass_closure_rel"]) < 1e-6
    assert res.element_results["ch_main"].model_dump()["second_law"] == "pass"
    # the inlet is wired straight to a plenum: the stream names that plenum,
    # not the builder's link between them
    inlet = [row for row in res.ledger["streams"] if row["node"] == "mfb"]
    assert inlet and inlet[0]["m_dot_in"] == pytest.approx(0.1)
    assert inlet[0]["via_owner"] == "plenum" and inlet[0]["direct"] is True


def test_a_failing_hidden_step_is_flagged_on_its_tee(monkeypatch) -> None:
    step = f"__tee_step__{TEE}_straight"
    original = second_law.SecondLawReport.element_verdicts

    def with_a_failure(self):
        out = original(self)
        out[step] = {**out[step], "verdict": "fail", "sigma_norm": -1e-3, "sigma_norm_high": -1e-3}
        return out

    monkeypatch.setattr(second_law.SecondLawReport, "element_verdicts", with_a_failure)
    res = _solve()
    assert [(f.element_id, f.owner_id) for f in res.second_law_failures] == [(step, TEE)]
    # report only
    assert res.success
