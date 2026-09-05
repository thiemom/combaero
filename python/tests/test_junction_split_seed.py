"""The initial guess at a junction must conserve mass and follow the pressures.

`analytical_pt_prop`'s Bernoulli mass-flow seed only reached `ChannelElement`.
Junctions in both the validation harness and `gui/backend/graph_builder.py` are
wired with `LosslessConnectionElement`, which got nothing, so every port fell
back to the flat reference flow: on a three-port junction that is 0.1 kg/s in
against 0.2 kg/s out, with a 1:1 split whatever the boundary pressures say
(issue #271, defect 11).

Measured effect on the junction scorecard, seed off then on:

    converged                    1625 -> 1732
    rejected as inadmissible      227 ->  163
    on-point evidence            1243 -> 1254
    RMSE on the common subset   unchanged to four decimals

The seed is opt-in per element class, because a junction subclass whose port
flows follow its own internal physics must not have them overwritten -- see
`EjectorElement`.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.components import MultiPortChamberElement
from combaero.network.ejector_element import EjectorElement
from combaero.network.mpce_v2_element import MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_A = 0.01


def _net(pt_str: float, pt_bra: float, area_bra: float = _A) -> FlowNetwork:
    """Common inlet at 1 bar, two pressure-driven outlets."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=1.0e5, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_str", Pt=pt_str, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_bra", Pt=pt_bra, Tt=300.0, Y=_Y))
    net.add_node(MomentumChamberNode("port_com", area=_A))
    net.add_node(MomentumChamberNode("port_str", area=_A))
    net.add_node(MomentumChamberNode("port_bra", area=area_bra))
    net.add_element(LosslessConnectionElement("lc_com", "pb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "pb_bra"))
    net.add_element(
        MPCEv2Element(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 45.0],
            port_areas=[_A, _A, area_bra],
            flow_direction="branch",
            strict=False,
        )
    )
    return net


def _x0(net: FlowNetwork) -> dict[str, float]:
    solver = NetworkSolver(net)
    solver.network.resolve_all_topology()
    solver.network.validate()
    ref = solver._infer_reference_state()
    solver._init_overrides = solver._propagate_analytical_pt_prop(ref)
    x0 = solver._build_x0()
    return dict(zip(solver.unknown_names, x0, strict=True))


def test_the_start_state_conserves_mass_at_the_junction():
    """The defect itself: 0.1 kg/s in against 0.2 kg/s out."""
    x0 = _x0(_net(pt_str=99_980.0, pt_bra=99_950.0))

    m_in = x0["lc_com.m_dot"]
    m_out = x0["lc_str.m_dot"] + x0["lc_bra.m_dot"]
    assert m_out == pytest.approx(m_in, rel=1e-9), (
        f"x0 leaks {m_out - m_in:.4g} kg/s at the junction"
    )


def test_the_start_split_follows_the_propagated_pressures():
    """The leg with the larger propagated pressure drop gets the larger seed,
    in the Bernoulli ratio sqrt(dPt_bra / dPt_str).

    Against the PROPAGATED differences, not the boundary ones: the seed reads
    `_propagate_pressure_guess`, which walks outward from the boundaries with a
    uniform step, so the port-face estimates are not the boundary values. That
    is the contract being pinned -- the seed uses the best pressure estimate
    available at x0, whatever the propagator's own accuracy.
    """
    net = _net(pt_str=99_990.0, pt_bra=99_960.0)
    solver = NetworkSolver(net)
    solver.network.resolve_all_topology()
    solver.network.validate()
    p_guess = solver._propagate_pressure_guess(solver._infer_reference_state())
    d_str = abs(p_guess["port_com"] - p_guess["port_str"])
    d_bra = abs(p_guess["port_com"] - p_guess["port_bra"])
    assert d_str > 0.0 and d_bra > d_str, "fixture must give the legs different drops"

    x0 = _x0(net)

    ratio = x0["lc_bra.m_dot"] / x0["lc_str.m_dot"]
    assert ratio == pytest.approx(math.sqrt(d_bra / d_str), rel=1e-6)


def test_equal_pressures_fall_back_to_an_area_weighted_split():
    """With nothing to distinguish the legs, area is the incompressible
    answer -- not a 1:1 split between unequal areas."""
    x0 = _x0(_net(pt_str=99_970.0, pt_bra=99_970.0, area_bra=_A / 3.0))

    ratio = x0["lc_bra.m_dot"] / x0["lc_str.m_dot"]
    assert ratio == pytest.approx(1.0 / 3.0, rel=1e-6)


def test_every_port_is_seeded_in_its_own_canonical_direction():
    """`port_mdots[i] = port_signs[i] * outer_mdot[i]`, so a positive outer
    flow is the correct direction at that port. A negative seed would start
    Newton inside the soft-barrier region."""
    x0 = _x0(_net(pt_str=99_980.0, pt_bra=99_950.0))

    for key in ("lc_com.m_dot", "lc_str.m_dot", "lc_bra.m_dot"):
        assert x0[key] > 0.0, f"{key} seeded against its declared direction"


def test_the_total_comes_from_the_existing_propagator():
    """The seed fixes the SPLIT and the CONTINUITY; it introduces no new
    flow-scale heuristic, so the level still comes from whatever the
    topological propagator put on the common port."""
    net = _net(pt_str=99_980.0, pt_bra=99_950.0)
    solver = NetworkSolver(net)
    solver.network.resolve_all_topology()
    solver.network.validate()
    ref = solver._infer_reference_state()
    # `_propagate_mdot_guess` is keyed by ELEMENT ID, not by unknown name.
    expected = abs(solver._propagate_mdot_guess(ref)["lc_com"])

    assert _x0(net)["lc_com.m_dot"] == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# Opt-in
# ---------------------------------------------------------------------------


def test_the_chamber_junction_opts_in():
    assert MultiPortChamberElement.seeds_ports_by_pressure_split is True
    assert MPCEv2Element.seeds_ports_by_pressure_split is True


def test_the_ejector_opts_out():
    """Its port flows follow choking and entrainment, not a Bernoulli share of
    the imposed pressure differences, and it carries its own warm start. With
    the seed applied to it the cold reference network in test_gui_ejector.py
    stops converging.
    """
    assert EjectorElement.seeds_ports_by_pressure_split is False


def test_an_opted_out_junction_gets_no_port_seeds(monkeypatch):
    net = _net(pt_str=99_980.0, pt_bra=99_950.0)
    monkeypatch.setattr(
        type(net.elements["jct"]), "seeds_ports_by_pressure_split", False, raising=False
    )
    x0 = _x0(net)

    # Falls back to the flat reference flow on every port, which is the
    # pre-fix behaviour and is what the ejector needs left alone.
    assert x0["lc_str.m_dot"] == pytest.approx(x0["lc_bra.m_dot"])


def test_a_channel_wired_junction_still_gets_a_forward_seed():
    """The pre-existing ChannelElement seed must be untouched: the junction
    seed adds a case, it does not replace one. Here the straight leg is a
    channel, so it carries the older Bernoulli override, and the junction seed
    must not fight it into a reversed start."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=1.0e5, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_str", Pt=99_980.0, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_bra", Pt=99_950.0, Tt=300.0, Y=_Y))
    net.add_node(MomentumChamberNode("port_com", area=_A))
    net.add_node(MomentumChamberNode("port_str", area=_A))
    net.add_node(MomentumChamberNode("port_bra", area=_A))
    net.add_element(LosslessConnectionElement("lc_com", "pb_in", "port_com"))
    net.add_element(
        ChannelElement(
            "lc_str", "port_str", "pb_str", length=0.3, diameter=0.1, regime="incompressible"
        )
    )
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "pb_bra"))
    net.add_element(
        MPCEv2Element(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 45.0],
            port_areas=[_A, _A, _A],
            flow_direction="branch",
            strict=False,
        )
    )
    x0 = _x0(net)

    assert x0["lc_str.m_dot"] > 0.0
    assert x0["lc_bra.m_dot"] > 0.0
