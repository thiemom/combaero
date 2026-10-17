"""#619: a compressible ejector must not admit a zero-flow root with the
mixing-plane static P_py above both supplies.

The C-D nozzle closure returned zero flow (flat in P_py) for an exit static
above the supply, so R0 and R1 held at zero flow for EVERY such P_py, and R3
(the jet-pump discharge pin, kept finite there by the #481 floors) crossed
zero at some huge P_py. The stress families found 9 of 11 "converged"
compressible draws there (P_py ~ 1e13 Pa, all flows 0, success and
consistent). Above the supply the nozzle now passes reverse flow, strictly
decreasing in P_py, so zero flow pins P_py to each supply in turn and closes
only if the supplies are equal.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.ejector_element import EjectorElement

Y = list(cb.species.dry_air_mass())
D = 0.05
A = np.pi * D * D / 4
PT = {"bp": 1.5e6, "bs": 157994.5, "bo": 302273.9}


def _net() -> FlowNetwork:
    """The stress draw (gadgets.ejector, seed 3): 15 bar primary, 1.6 bar
    secondary, 3 bar outlet, critical back pressure 2.5 bar."""
    g = FlowNetwork()
    for nid, Tt in (("bp", 300.0), ("bs", 280.0), ("bo", 300.0)):
        g.add_node(PressureBoundary(nid, Pt=PT[nid], Tt=Tt, Y=Y))
    for nid, b in (("mp", "bp"), ("ms", "bs"), ("mo", "bo")):
        mc = MomentumChamberNode(nid, area=A)
        mc.initial_guess = {f"{nid}.P": PT[b], f"{nid}.Pt": PT[b]}
        g.add_node(mc)
    for eid, a, b, L in (("cp", "bp", "mp", 0.2), ("cs", "bs", "ms", 0.2), ("co", "mo", "bo", 1.0)):
        g.add_element(ChannelElement(eid, a, b, length=L, diameter=D, regime="incompressible"))
    g.add_element(
        EjectorElement(
            id="ej",
            primary_node="mp",
            secondary_node="ms",
            outlet_node="mo",
            throat_area=6e-5,
            nozzle_exit_area=1.6343e-4,
            mixing_area=1.03445e-3,
            regime="compressible",
        )
    )
    return g


def _zero_flow_state(s: NetworkSolver, p_py: float) -> np.ndarray:
    x = np.array(s._build_x0(), dtype=float)
    for i, name in enumerate(s.unknown_names):
        if name.endswith(".m_dot"):
            x[i] = 0.0
    x[s.unknown_names.index("ej.P_jct")] = p_py
    return x


@pytest.mark.parametrize("p_py", [1.6e6, 1e8, 1e10, 5.06e13])
def test_zero_flow_above_both_supplies_is_not_a_root(p_py: float) -> None:
    s = NetworkSolver(_net())
    s.network.resolve_all_topology()
    x = _zero_flow_state(s, p_py)
    rows = [
        i for i, (owner, _k) in enumerate(s._row_kinds(x)) if owner.startswith("EjectorElement")
    ]
    res = s._residuals(x)[rows]
    # R0 (primary) and R1 (secondary) are the reverse nozzle flows at that
    # P_py: at least the throat-choked reverse flow from 1.6 MPa, ~0.5 kg/s.
    assert abs(res[0]) > 0.1
    assert abs(res[1]) > 0.1
