"""The cold start from a quadratic-resistance network solve.

Every element offering a ``_seed_resistance`` K, the start is the network of
dPt = K m|m| with mass conserved: pressure boundaries fixed, mass-flow
boundaries injecting (a sink drawing) their m_dot. A network holding any
other element keeps the propagated start. On 160 stress random draws the
cold path went from 67 to 110 converged, and its total time roughly halved.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    IdealMomentumTeeElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _two_orifices() -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(OrificeElement("o1", "a", "p", Cd=0.7, diameter=0.01, correlation="fixed"))
    g.add_element(OrificeElement("o2", "p", "b", Cd=0.6, diameter=0.02, correlation="fixed"))
    return g


def test_the_seed_is_the_resistance_network_root() -> None:
    """Two orifices in series: m = sqrt(dP / (K1 + K2)), and the plenum
    sits where the first orifice takes K1 m^2 of the drop. To 1%: the
    conductance floor (SEED_MDOT_FLOOR) shifts a start of 0.03 kg/s by 0.2%."""
    g = _two_orifices()
    g.resolve_all_topology()
    s = NetworkSolver(g)
    ref = s._infer_reference_state()
    p_map, m_map = s._resistance_network_guess(ref)
    P_ref = 1.5e5
    rho = cb.density(ref["T"], P_ref, list(cb.mass_to_mole(ref["Y"])))
    K1, K2 = (
        1.0 / (2.0 * rho * (Cd * math.pi * d * d / 4.0) ** 2)
        for Cd, d in ((0.7, 0.01), (0.6, 0.02))
    )
    m = math.sqrt(1.0e5 / (K1 + K2))
    assert m_map["o1"] == pytest.approx(m, rel=1e-2)
    assert m_map["o2"] == pytest.approx(m_map["o1"], rel=1e-9)
    assert p_map["p"] == pytest.approx(2.0e5 - K1 * m * m, rel=1e-3)


def test_a_sink_draws_its_flow() -> None:
    """A mass-flow boundary with only upstream elements takes its m_dot out."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("s", m_dot=0.05, Tt=300.0, Y=Y_AIR))
    g.add_element(OrificeElement("o", "a", "s", Cd=0.7, diameter=0.01, correlation="fixed"))
    g.resolve_all_topology()
    s = NetworkSolver(g)
    p_map, m_map = s._resistance_network_guess(s._infer_reference_state())
    assert m_map["o"] == pytest.approx(0.05, rel=1e-5)
    assert p_map["s"] < 2.0e5


def test_an_element_without_a_seed_keeps_the_propagated_start() -> None:
    """An ideal tee offers no seed resistance: no seed, and x0 is the
    propagated start bit for bit."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    for n in ("pi", "po0", "po1"):
        g.add_node(MomentumChamberNode(n, area=2e-3))
    g.add_node(PressureBoundary("b0", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("b1", Pt=1.1e5, Tt=300.0, Y=Y_AIR))
    g.add_element(OrificeElement("oi", "a", "pi", Cd=0.7, diameter=0.04, correlation="fixed"))
    g.add_element(OrificeElement("oo0", "po0", "b0", Cd=0.7, diameter=0.04, correlation="fixed"))
    g.add_element(OrificeElement("oo1", "po1", "b1", Cd=0.7, diameter=0.04, correlation="fixed"))
    g.add_element(
        IdealMomentumTeeElement(
            id="t",
            inlet_nodes=["pi"],
            outlet_nodes=["po0", "po1"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[2e-3, 2e-3, 2e-3],
        )
    )
    g.resolve_all_topology()
    s = NetworkSolver(g)
    assert s._resistance_network_guess(s._infer_reference_state()) is None
    x_on = s._build_x0()
    s.resistance_seed = False
    assert np.array_equal(x_on, s._build_x0())


@pytest.mark.slow
def test_a_stress_draw_the_propagated_start_could_not_reach() -> None:
    """Stress random 0 (24 elements): the propagated start fails with every
    retry; from the resistance seed it converges."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    out = {}
    for seeded in (False, True):
        net, _ = FAMILIES["random"](0)
        s = NetworkSolver(net)
        s.resistance_seed = seeded
        s.homotopy_retry = True
        s.bc_ramp_retry = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[seeded] = s.solve(timeout=20.0)["__converged__"]
    assert out == {False: False, True: True}
