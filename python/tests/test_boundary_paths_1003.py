"""A path of flow-independent elements between two pressure boundaries is
refused (#1003).

Every PressureBoundary fixes its Pt, so a path made only of elements whose Pt
relation ignores their own flow (``_flow_independent_dpt``: lossless links,
fraction losses without m_ref, free vortices, fixed-ratio turbomachines) has
its pressure ratio fixed at both ends. Two outlets fed from one plenum by
plain GUI edges are such a path: with equal outlet Pt the split between them
is free (perturbed starts converged with success at four different splits,
two with an outlet flowing backwards), and with outlet Pt 0.1 % apart there
is no root (the solve fails at +-5.7e11 kg/s through the pair).

The edge-set elements couple the boundary's Pt, not its static pressure:
their rows read the boundary's Pt directly, and only channels, crossflow
segments and port discharges charge an exit head (``exit_head_lost``). A
single flow-independent link into one boundary is fine: it has one fixed end.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from test_gui_http_routes import _assert_readable_error, _edge, _node, _post
from test_lossless_loop_988 import _orifice, _pb

from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    NetworkSolver,
    PlenumNode,
    PressureLossElement,
    VortexElement,
)
from combaero.network.pressure_loss import ConstantFractionLoss


def _two_outlets(path1, path2, p2: float = 1e5) -> FlowNetwork:
    """in (2 bar) -> orifice -> p; p -> path1 -> out1; p -> path2 -> out2.
    Each path is a list of (kind, id) run from p towards its outlet."""
    g = FlowNetwork()
    g.add_node(_pb("in", 2e5))
    g.add_node(_pb("out1", 1e5))
    g.add_node(_pb("out2", p2))
    g.add_node(PlenumNode("p"))
    g.add_element(_orifice("o", "in", "p"))
    for out, path in (("out1", path1), ("out2", path2)):
        a = "p"
        for i, (kind, eid) in enumerate(path):
            b = out if i == len(path) - 1 else f"{eid}_n"
            if b != out:
                g.add_node(PlenumNode(b))
            if kind == "lossless":
                g.add_element(LosslessConnectionElement(eid, a, b))
            elif kind == "fraction":
                g.add_element(
                    PressureLossElement(eid, a, b, correlation=ConstantFractionLoss(xi=0.02))
                )
            elif kind == "vortex":
                g.add_element(VortexElement(eid, a, b, r_in=0.05, r_out=0.1))
            else:
                g.add_element(_orifice(eid, a, b))
            a = b
    return g


def _refused(g: FlowNetwork, members: list[str], ends: tuple[str, str]) -> str:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="join the pressure boundaries") as exc:
            g.validate()
    msg = str(exc.value)
    named = sorted(eid for eid in g.elements if f"{eid!r}" in msg)
    assert named == sorted(members), msg
    assert f"{ends[0]!r} and {ends[1]!r}" in msg, msg
    assert "flow-dependent loss" in msg and "merge the boundaries" in msg, msg
    return msg


def _solves(g: FlowNetwork) -> dict:
    g.validate()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__converged__"] and r["__success__"], r.get("__message__")
    return r


@pytest.mark.parametrize("p2", [1e5, 1.001e5])  # equal: split free; unequal: no root
def test_two_outlets_fed_by_lossless_links_are_refused(p2: float) -> None:
    g = _two_outlets([("lossless", "l1")], [("lossless", "l2")], p2=p2)
    msg = _refused(g, ["l1", "l2"], ("out1", "out2"))
    assert "(p -> out1, lossless)" in msg and "(p -> out2, lossless)" in msg


def test_fraction_losses_without_m_ref_between_boundaries_are_refused() -> None:
    g = _two_outlets([("fraction", "f1")], [("fraction", "f2")])
    msg = _refused(g, ["f1", "f2"], ("out1", "out2"))
    assert "ConstantFractionLoss without m_ref" in msg


def test_a_mixed_path_is_refused_in_order_and_off_path_elements_unnamed() -> None:
    """out1 <- lossless <- p -> fraction -> n -> vortex -> out2; the inlet
    orifice touches p but is not on the path."""
    g = _two_outlets([("lossless", "l1")], [("fraction", "f2"), ("vortex", "v2")])
    msg = _refused(g, ["l1", "f2", "v2"], ("out1", "out2"))
    assert msg.index("'l1'") < msg.index("'f2'") < msg.index("'v2'")


def test_a_fraction_loss_straight_between_two_boundaries_is_refused() -> None:
    g = FlowNetwork()
    g.add_node(_pb("in", 2e5))
    g.add_node(_pb("out", 1e5))
    g.add_element(PressureLossElement("f", "in", "out", correlation=ConstantFractionLoss(xi=0.1)))
    _refused(g, ["f"], ("in", "out"))


def test_two_inlets_joined_by_lossless_links_are_refused() -> None:
    """The same at the supply side: two reservoirs into one plenum."""
    g = FlowNetwork()
    g.add_node(_pb("in1", 2e5))
    g.add_node(_pb("in2", 2e5))
    g.add_node(_pb("out", 1e5))
    g.add_node(PlenumNode("p"))
    g.add_element(LosslessConnectionElement("a", "in1", "p"))
    g.add_element(LosslessConnectionElement("b", "in2", "p"))
    g.add_element(_orifice("o", "p", "out"))
    _refused(g, ["a", "b"], ("in1", "in2"))


def test_an_orifice_on_one_outlet_path_solves() -> None:
    """The flow-dependent loss on one path sets the split: p sits at out1's
    Pt, and the orifice passes what 1.0 -> 0.95 bar drives through it. A
    perturbed start lands on the same split."""

    def net() -> FlowNetwork:
        return _two_outlets([("lossless", "l1")], [("orifice", "o2")], p2=0.95e5)

    r = _solves(net())
    assert r["l1.m_dot"] > 0.0 and r["o2.m_dot"] > 0.0
    assert r["o.m_dot"] == pytest.approx(r["l1.m_dot"] + r["o2.m_dot"], rel=1e-6)
    x = np.asarray(r["__x_solution__"], dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = NetworkSolver(net()).solve(
            x0=x * (1 + 0.3 * np.random.default_rng(0).standard_normal(x.size))
        )
    assert r2["__converged__"]
    assert r2["o2.m_dot"] == pytest.approx(r["o2.m_dot"], rel=1e-5)


def test_a_single_lossless_link_into_one_boundary_solves() -> None:
    g = FlowNetwork()
    g.add_node(_pb("in", 2e5))
    g.add_node(_pb("out", 1e5))
    g.add_node(PlenumNode("p"))
    g.add_element(_orifice("o", "in", "p"))
    g.add_element(LosslessConnectionElement("l", "p", "out"))
    r = _solves(g)
    assert r["l.m_dot"] == pytest.approx(r["o.m_dot"], rel=1e-6)
    assert r["p.Pt"] == pytest.approx(1e5, rel=1e-9)  # the link carries out's Pt


def test_a_lossless_loop_through_one_boundary_keeps_the_988_message() -> None:
    """Out of a boundary and back into it is a loop of its own."""
    g = FlowNetwork()
    g.add_node(_pb("in", 2e5))
    g.add_node(_pb("out", 1e5))
    g.add_node(PlenumNode("p"))
    g.add_element(_orifice("o", "in", "p"))
    g.add_element(LosslessConnectionElement("la", "p", "out"))
    g.add_element(LosslessConnectionElement("lb", "p", "out"))
    with pytest.raises(ValueError, match="the lossless connections") as exc:
        g.validate()
    assert "'la'" in str(exc.value) and "'lb'" in str(exc.value)


def test_two_plain_gui_edges_to_two_outlets_are_refused_readably() -> None:
    """The GUI turns an untyped edge between two nodes into a lossless auto
    link: a plenum with plain edges to two outlets builds the path."""
    graph = {
        "nodes": [
            _node("inlet", "pressure_boundary", {"Pt": 2.0e5, "Tt": 300.0}),
            _node("ch", "channel", {"L": 0.5, "D": 0.05}),
            _node("p", "plenum"),
            _node("out1", "pressure_boundary", {"Pt": 1.0e5, "Tt": 300.0}),
            _node("out2", "pressure_boundary", {"Pt": 1.0e5, "Tt": 300.0}),
        ],
        "edges": [
            _edge("e1", "inlet", "ch"),
            _edge("e2", "ch", "p"),
            _edge("e3", "p", "out1"),
            _edge("e4", "p", "out2"),
        ],
        "solver_settings": {},
    }
    _assert_readable_error(
        _post("/solve", graph), "join the pressure boundaries", "p -> out1", "p -> out2"
    )
