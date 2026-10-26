"""A tee port meets a channel of a different bore through a sharp step (#536).

The tee's ports keep the closures' geometry (straight port = F_C), so a
straight channel narrower than the main pipe used to join a port of the main
area with nothing in between. The ideal tee, which couples on static pressure,
then raised the channel's Pt by its dynamic head and the solve was rejected;
Mynard, which couples on Pt, contracted losslessly. The GUI now inserts what
two bolted bores of different size are: a sharp contraction or expansion.

The network is the reported one, rebuilt here (the user's export is not
committed): mass-flow inlet 0.1 kg/s -> plenum -> channel D 0.1 -> dividing
tee; straight: channel D_s -> outlet; branch: channel D 0.05 -> orifice ->
channel D 0.1 -> the same outlet.
"""

from __future__ import annotations

import math
import warnings

import pytest

from combaero.network.components import AreaChangeElement
from gui.backend.graph_builder import TEE_PORT_STEP_TOL, build_network_from_schema
from gui.backend.runner import NetworkRunner
from gui.backend.schemas import NetworkGraphSchema

TEE = "tee"


def _node(i: str, t: str, **data) -> dict:
    return {"id": i, "type": t, "position": {"x": 0.0, "y": 0.0}, "data": data}


def _edge(i: str, s: str, t: str, sh: str | None = None, th: str | None = None) -> dict:
    return {"id": i, "source": s, "target": t, "sourceHandle": sh, "targetHandle": th, "data": {}}


def _schema(model: str, D_straight: float) -> dict:
    nodes = [
        _node("mfb", "mass_boundary", m_dot=0.1, Tt=300.0),
        _node("plenum", "plenum"),
        _node("ch_main", "channel", L=1.0, D=0.1),
        _node(TEE, "mpce_tee", junction_model=model, flow_direction="branch", theta_deg=90.0),
        _node("ch_str", "channel", L=1.0, D=D_straight),
        _node("ch_bra", "channel", L=1.0, D=0.05),
        _node("orf", "orifice", Cd=0.6, diameter=0.08),
        _node("ch_ret", "channel", L=1.0, D=0.1),
        _node("out", "pressure_boundary", Pt=101325.0, Tt=300.0),
    ]
    edges = [
        _edge("e1", "mfb", "plenum", "flow-source", "flow-left"),
        _edge("e2", "plenum", "ch_main", "flow-right", "flow-target"),
        _edge("e3", "ch_main", TEE, "flow-source", "port-common-target"),
        _edge("e4", TEE, "ch_str", "port-straight-source", "flow-target"),
        _edge("e5", TEE, "ch_bra", "port-branch-source", "flow-target"),
        _edge("e6", "ch_str", "out", "flow-source", "flow-target"),
        _edge("e7", "ch_bra", "orf", "flow-source", "flow-target"),
        _edge("e8", "orf", "ch_ret", "flow-source", "flow-target"),
        _edge("e9", "ch_ret", "out", "flow-source", "flow-target"),
    ]
    return {"nodes": nodes, "edges": edges}


def _steps(net) -> dict[str, AreaChangeElement]:
    return {
        eid: e
        for eid, e in net.elements.items()
        if isinstance(e, AreaChangeElement) and eid.startswith("__tee_step__")
    }


def _solve(model: str, D_straight: float):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkRunner.from_dict(_schema(model, D_straight)).solve(
            method="hybr", init_strategy="analytical_pt_prop"
        )


def test_a_narrower_straight_channel_gets_a_sharp_contraction() -> None:
    net = build_network_from_schema(NetworkGraphSchema.model_validate(_schema("mynard", 0.05)))
    steps = _steps(net)
    assert list(steps) == [f"__tee_step__{TEE}_straight"]
    step = steps[f"__tee_step__{TEE}_straight"]
    assert step.model_type == "sharp"
    # Drawn tee -> channel: from the tee's port (main area) to the channel.
    assert step.from_node == f"__tee_jct__{TEE}_straight"
    assert pytest.approx(math.pi / 4 * 0.1**2) == step.F0
    assert pytest.approx(math.pi / 4 * 0.05**2) == step.F1
    # The channel now starts on the step's far side, not on the tee's port.
    assert net.elements["ch_str"].from_node == step.to_node


@pytest.mark.parametrize("D_straight", [0.1, 0.1 * math.sqrt(1.0 - 0.9 * TEE_PORT_STEP_TOL)])
def test_no_step_where_the_bores_match(D_straight: float) -> None:
    net = build_network_from_schema(
        NetworkGraphSchema.model_validate(_schema("mynard", D_straight))
    )
    assert _steps(net) == {}
    assert net.elements["ch_str"].from_node == f"__tee_jct__{TEE}_straight"


def test_every_step_the_diagnostic_would_flag_is_inserted() -> None:
    """The builder and the junction's unmodelled-step diagnostic (#539) share
    one bound on |A_channel / A_port - 1|. An expansion just past it got no
    step while the diagnostic still flagged it, when the builder measured the
    difference against the larger area instead."""
    ratio = 1.0 + 1.02 * TEE_PORT_STEP_TOL
    net = build_network_from_schema(
        NetworkGraphSchema.model_validate(_schema("mynard", 0.1 * math.sqrt(ratio)))
    )
    assert list(_steps(net)) == [f"__tee_step__{TEE}_straight"]


def test_a_feeding_channel_steps_into_the_port() -> None:
    """Merge orientation: the channel feeds the port, so the step runs from
    the channel side (its own area) into the port."""
    raw = _schema("mynard", 0.05)
    tee = next(n for n in raw["nodes"] if n["id"] == TEE)
    tee["data"]["flow_direction"] = "merge"
    for e in raw["edges"]:
        if e["id"] == "e4":
            e.update(source="ch_str", target=TEE, sourceHandle="flow-source")
            e["targetHandle"] = "port-straight-target"
        if e["id"] == "e6":
            e.update(source="out", target="ch_str")
    net = build_network_from_schema(NetworkGraphSchema.model_validate(raw))
    step = _steps(net)[f"__tee_step__{TEE}_straight"]
    assert step.to_node == f"__tee_jct__{TEE}_straight"
    assert pytest.approx(math.pi / 4 * 0.05**2) == step.F0
    assert pytest.approx(math.pi / 4 * 0.1**2) == step.F1


@pytest.mark.parametrize("model", ["ideal_momentum", "mynard"])
def test_the_reported_network_solves_and_the_tee_reports_its_step(model: str) -> None:
    """The ideal tee was rejected here (net energy outflow, K_straight -7.4)."""
    res = _solve(model, 0.05)
    assert res.success, res.message
    assert res.consistent is True
    tee = res._element_results[TEE].model_dump()
    assert tee["step_straight_connected_area_ratio"] == pytest.approx(0.25)
    # Weisbach sharp contraction, 0.5 (1 - 0.25), referenced to the channel,
    # at M -> 0; the compressible vena contracta (#731) reads 1.2% above it
    # at this step's Mach (measured 2026-10-25: 0.3797, 0.3792).
    assert 0.375 < tee["step_straight_zeta"] < 0.375 * 1.03
    assert tee["step_straight_dPt"] > 0.0
    assert "step_common_zeta" not in tee and "step_branch_zeta" not in tee
    # The step is modelled, so the unmodelled-step warning (#539) stays off.
    assert tee["port_area_mismatch"] == 0


def test_the_closures_now_agree_on_the_split() -> None:
    """Measured, not tuned: with the step in place the parameter-free ideal tee
    and Mynard put 0.0518 and 0.0495 kg/s down the straight leg (before: the
    ideal tee rejected, Mynard 0.0526 through a lossless contraction)."""
    m = {
        model: _solve(model, 0.05)._element_results[TEE].model_dump()["m_dot_straight"]
        for model in ("ideal_momentum", "mynard")
    }
    assert m["ideal_momentum"] == pytest.approx(m["mynard"], rel=0.06)
