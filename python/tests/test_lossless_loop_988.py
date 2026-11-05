"""A closed loop of lossless connections is refused (#988).

Lossless connections hold Pt equal at both ends, so around a loop made only
of them any circulating flow meets every row: the steady problem has a
one-parameter family of roots (an exact null vector of the Jacobian) and the
solver returned whichever its start lay nearest -- on stress random 5, 159
kg/s circulating while the network carried 8.7 kg/s, then 1.3/1.5 kg/s after
a seed change. In inviscid flow the circulation keeps whatever it started
with; a steady network has no history to set it, so ``validate()`` refuses the
loop and names it. A loop with any lossy element is determined and passes.
"""

from __future__ import annotations

import pytest
from test_gui_http_routes import _assert_readable_error, _edge, _node, _post

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _pb(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, 300.0, Y_AIR
    return b


def _orifice(eid: str, a: str, b: str) -> OrificeElement:
    return OrificeElement(eid, a, b, Cd=0.7, diameter=0.02, correlation="fixed")


def _between(*links: tuple[str, str, str]) -> FlowNetwork:
    """in -> orifice -> p1 ... links among p1, p2, p3 ... p2 -> orifice -> out."""
    g = FlowNetwork()
    g.add_node(_pb("in", 2e5))
    g.add_node(_pb("out", 1e5))
    used = {"p1", "p2"} | {n for _, a, b in links for n in (a, b)}
    for n in sorted(used):
        g.add_node(PlenumNode(n))
    g.add_element(_orifice("o_in", "in", "p1"))
    g.add_element(_orifice("o_out", "p2", "out"))
    for eid, a, b in links:
        if eid.startswith("o"):
            g.add_element(_orifice(eid, a, b))
        else:
            g.add_element(LosslessConnectionElement(eid, a, b))
    return g


@pytest.mark.parametrize(
    "links, loop",
    [
        # two in parallel, either direction
        ([("la", "p1", "p2"), ("lb", "p1", "p2"), ("l3", "p2", "p3"), ("l4", "p3", "p2")], None),
        (
            [("la", "p1", "p2"), ("lb", "p2", "p1"), ("l3", "p1", "p3"), ("o3", "p3", "p2")],
            ["la", "lb"],
        ),
        # a ring of three through plenums
        ([("la", "p1", "p2"), ("lb", "p2", "p3"), ("lc", "p3", "p1")], ["la", "lb", "lc"]),
    ],
)
def test_a_lossless_only_loop_is_refused_and_named(links: list, loop: list[str] | None) -> None:
    g = _between(*links)
    with pytest.raises(ValueError, match="closed loop") as exc:
        g.validate()
    msg = str(exc.value)
    named = [eid for eid, _, _ in links if f"{eid!r}" in msg]
    if loop is not None:
        assert named == loop, msg
    else:
        assert len(named) == 2, msg  # one loop, its two members
    assert "Remove one of the connections" in msg


def test_a_loop_with_a_loss_element_is_determined_and_solves() -> None:
    """Lossless in parallel with an orifice: the orifice carries none, the
    split is fixed, and nothing is refused."""
    g = _between(("la", "p1", "p2"), ("o_par", "p1", "p2"))
    g.validate()
    r = NetworkSolver(g).solve()
    assert r["__success__"]
    assert abs(r["o_par.m_dot"]) < 1e-9 * abs(r["o_in.m_dot"])


def test_a_lossless_tree_is_not_a_loop() -> None:
    """A chain of lossless links (no cycle) through a third plenum passes."""
    g = _between(("la", "p1", "p3"), ("lb", "p3", "p2"))
    g.validate()
    assert NetworkSolver(g).solve()["__success__"]


def test_two_plain_gui_edges_back_and_forth_are_refused_readably() -> None:
    """The GUI turns an untyped edge into a lossless auto link, so drawing
    p1 -> p2 and p2 -> p1 builds the loop: the solve says so."""
    graph = {
        "nodes": [
            _node("inlet", "pressure_boundary", {"Pt": 2.0e5, "Tt": 300.0}),
            _node("ch", "channel", {"L": 0.5, "D": 0.05}),
            _node("p1", "plenum"),
            _node("p2", "plenum"),
            _node("ch2", "channel", {"L": 0.5, "D": 0.05}),
            _node("outlet", "pressure_boundary", {"Pt": 1.9e5, "Tt": 300.0}),
        ],
        "edges": [
            _edge("e1", "inlet", "ch"),
            _edge("e2", "ch", "p1"),
            _edge("e3", "p1", "p2"),
            _edge("e4", "p2", "p1"),
            _edge("e5", "p2", "ch2"),
            _edge("e6", "ch2", "outlet"),
        ],
        "solver_settings": {},
    }
    _assert_readable_error(_post("/solve", graph), "closed loop", "p1 -> p2", "p2 -> p1")
