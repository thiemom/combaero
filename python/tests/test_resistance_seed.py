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
    ChannelElement,
    FlowNetwork,
    IdealMomentumTeeElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network import solver as solver_module
from combaero.network.components import SEED_FRICTION_FACTOR

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _two_orifices(regime: str = "incompressible") -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    for eid, a, b, Cd, d in (("o1", "a", "p", 0.7, 0.01), ("o2", "p", "b", 0.6, 0.02)):
        g.add_element(
            OrificeElement(eid, a, b, Cd=Cd, diameter=d, correlation="fixed", regime=regime)
        )
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


def test_a_compressible_element_takes_its_kernels_flow(monkeypatch: pytest.MonkeyPatch) -> None:
    """The same orifices compressible: o1 chokes on the 2:1 drop. The seed
    re-fits each element past the flux gate to the flow its own kernel passes
    (#971) and lands within 1% of the solved root; the incompressible K alone
    starts 21% high. Falsified by switching the re-fit off."""
    root = NetworkSolver(_two_orifices("compressible")).solve()
    assert root["__success__"]

    def seed() -> float:
        g = _two_orifices("compressible")
        g.resolve_all_topology()
        s = NetworkSolver(g)
        return s._resistance_network_guess(s._infer_reference_state())[1]["o1"]

    assert seed() == pytest.approx(root["o1.m_dot"], rel=1e-2)
    monkeypatch.setattr(solver_module, "SEED_COMPRESSIBLE_PASSES", 0)
    assert seed() > 1.1 * root["o1.m_dot"]


def _duct_chain(n: int) -> FlowNetwork:
    """The D-7156 hot side (#971), adiabatic and in air: 1.62 kg/s at 589 K
    through a 6 cm lead duct, n 1.22 cm segments between plenums and a 1 cm
    tail into 1 atm."""
    area, dh = 0.037, 0.119
    D = math.sqrt(4.0 * area / math.pi)
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("in", m_dot=1.62, Tt=589.0, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=101325.0, Tt=589.0, Y=Y_AIR))
    for i in range(n + 1):
        g.add_node(PlenumNode(f"g{i}"))

    def duct(eid: str, a: str, b: str, length: float) -> ChannelElement:
        return ChannelElement(eid, a, b, length=length, diameter=D, Dh=dh, roughness=0.0)

    g.add_element(duct("lead", "in", "g0", 0.06))
    for i in range(1, n + 1):
        g.add_element(duct(f"h{i}", f"g{i - 1}", f"g{i}", 0.0122))
    g.add_element(duct("tail", f"g{n}", "out", 0.01))
    return g


def test_a_channel_charges_its_exit_head_only_where_its_residual_does() -> None:
    """The exit head is lost only into a static-coupled pressure boundary
    (#360): a segment between plenums is friction alone, the tail into the
    outlet friction plus one head (#971)."""
    g = _duct_chain(3)
    g.resolve_all_topology()
    rho = 0.6
    for eid in ("h1", "h2", "h3", "tail"):
        el = g.elements[eid]
        head = 1.0 if eid == "tail" else 0.0
        K = (SEED_FRICTION_FACTOR * el.length / el.Dh + head) / (2.0 * rho * el.area**2)
        assert el._seed_resistance(rho, 101325.0, 1.0) == pytest.approx(K, rel=1e-12)


def test_an_unresolved_channel_charges_its_exit_head() -> None:
    """Before topology is resolved the exit is unknown and the head is
    charged, as before #971: an estimate made while a network is built (the
    stress splitter sizes its supply from it) keeps its value."""
    g = _duct_chain(3)
    el = g.elements["h1"]
    rho = 0.6
    K = (SEED_FRICTION_FACTOR * el.length / el.Dh + 1.0) / (2.0 * rho * el.area**2)
    assert el._seed_resistance(rho, 101325.0, 1.0) == pytest.approx(K, rel=1e-12)


def test_a_chain_of_short_ducts_converges_cold() -> None:
    """Ten short segments between plenums (#971): with an exit head charged
    on every segment the seed started the chain 14 kPa high against a 13 Pa
    root, the cold solve stalled and the warm-start retry converged it. Now
    the start sits within 5% of the root's drop to the outlet (1.2% measured)
    and the cold attempt converges, without the retry."""
    s = NetworkSolver(_duct_chain(10))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = s.solve(timeout=20.0)
    assert r["__success__"]
    assert not [w for w in caught if "auto-retry" in str(w.message)]
    x0 = s._build_x0()
    i0 = s._name_to_index["g0.P"]
    assert abs(x0[i0] - r["g0.P"]) < 0.05 * (r["g0.P"] - 101325.0)


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
