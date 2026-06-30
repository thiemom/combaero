import pytest

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _get_air_Y():
    """Helper to create standard air composition."""
    Y = [0.0] * cb.num_species()
    Y[cb.species_index_from_name("N2")] = 0.79
    Y[cb.species_index_from_name("O2")] = 0.21
    return Y


def build_pressure_bounds():
    """Create standard pressure boundaries."""
    inlet = PressureBoundary("inlet", Pt=150000.0, Tt=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", Pt=100000.0, Tt=300.0, Y=_get_air_Y())
    return inlet, outlet


def test_step_1_simple_series():
    """Test: PB -> Channel -> Plenum -> Channel -> PB"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    p1 = ChannelElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    p2 = ChannelElement("p2", "n1", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1]:
        graph.add_node(n)
    for e in [p1, p2]:
        graph.add_element(e)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify solution exists and mass is conserved
    assert "p1.m_dot" in sol
    assert "p2.m_dot" in sol
    assert abs(sol["p1.m_dot"] - sol["p2.m_dot"]) < 1e-6
    assert sol["p1.m_dot"] > 0  # Positive flow from high to low pressure


def test_step_2_add_orifices():
    """Test: PB -> Channel -> n1 -> Orifice -> mc1 -> Channel -> PB"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    mc1 = MomentumChamberNode("mc1")

    p1 = ChannelElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1", Cd=0.6, diameter=0.079788, correlation="fixed")
    p2 = ChannelElement("p2", "mc1", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1, mc1]:
        graph.add_node(n)
    for e in [p1, o1, p2]:
        graph.add_element(e)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify mass conservation across all elements
    assert abs(sol["p1.m_dot"] - sol["o1.m_dot"]) < 1e-6
    assert abs(sol["o1.m_dot"] - sol["p2.m_dot"]) < 1e-6
    assert sol["p1.m_dot"] > 0


def test_step_3_full_series_no_bypass():
    """Test: Full user series without bypass"""
    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    mc1 = MomentumChamberNode("mc1")
    n2 = PlenumNode("n2")
    plen = PlenumNode("plenum")
    mc2 = MomentumChamberNode("mc2")

    p1 = ChannelElement("p1", "inlet", "n1", length=2.0, diameter=0.1, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1", Cd=0.6, diameter=0.100925, correlation="fixed")
    p2 = ChannelElement("p2", "mc1", "n2", length=2.0, diameter=0.1, roughness=1e-4)
    o2 = OrificeElement("o2", "n2", "plenum", Cd=0.6, diameter=0.100925, correlation="fixed")
    p3 = ChannelElement("p3", "plenum", "mc2", length=2.0, diameter=0.1, roughness=1e-4)
    p4 = ChannelElement("p4", "mc2", "outlet", length=2.0, diameter=0.1, roughness=1e-4)

    for n in [inlet, outlet, n1, mc1, n2, plen, mc2]:
        graph.add_node(n)
    for e in [p1, o1, p2, o2, p3, p4]:
        graph.add_element(e)

    # Set initial guesses for better convergence
    n1.initial_guess = {"n1.P": 140000.0, "n1.Pt": 140000.0}
    mc1.initial_guess = {"mc1.P": 135000.0, "mc1.Pt": 135000.0}
    n2.initial_guess = {"n2.P": 130000.0, "n2.Pt": 130000.0}
    plen.initial_guess = {"plenum.P": 120000.0, "plenum.Pt": 120000.0}
    mc2.initial_guess = {"mc2.P": 110000.0, "mc2.Pt": 110000.0}

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Verify mass conservation across entire series
    mass_flows = [sol[f"{e}.m_dot"] for e in ["p1", "o1", "p2", "o2", "p3", "p4"]]
    for i in range(len(mass_flows) - 1):
        assert abs(mass_flows[i] - mass_flows[i + 1]) < 1e-6, (
            f"Mass not conserved between elements {i} and {i + 1}"
        )

    assert all(mf > 0 for mf in mass_flows), "All mass flows should be positive"


def test_step_4_adding_bypass():
    """Bypass branch: splits at mc1, rejoins at mc2.

    mc1 (splitting) and mc2 (merging) are momentum-CV junctions built from
    MultiPortChamberElement + MomentumChamberNode port faces. The earlier
    incarnation of this test used bare MomentumChamberNode at both splits
    and was xfailed against #174 because MCN's scalar Pt = P + 0.5 rho v^2
    cannot represent merging or splitting streams in the chamber itself.
    Migration onto the momentum-CV junction (closes #174) restores the
    expected diameter-driven split: smaller-bore bypass < main branch.
    """
    import math

    D_main = 0.1
    D_bypass = 0.08
    A_main = math.pi * (D_main / 2.0) ** 2
    A_bypass = math.pi * (D_bypass / 2.0) ** 2

    graph = FlowNetwork()
    inlet, outlet = build_pressure_bounds()

    n1 = PlenumNode("n1")
    n2 = PlenumNode("n2")
    plen = PlenumNode("plenum")
    n_bypass = PlenumNode("n_bypass")

    # mc1: splitting junction (1 inflow from o1, 2 outflows to p2 and p_bypass).
    # Port-face MCNs: com (inlet) + str (straight outlet) + bra (branch outlet,
    # 90 deg). BorderCarnotLossElement carries the angled-arm loss on bra.
    mc1_com = MomentumChamberNode("mc1_com", area=A_main)
    mc1_str = MomentumChamberNode("mc1_str", area=A_main)
    mc1_bra = MomentumChamberNode("mc1_bra", area=A_bypass)
    mc1_bra_post = MomentumChamberNode("mc1_bra_post", area=A_bypass)
    loss_mc1_bra = BorderCarnotLossElement(
        "loss_mc1_bra",
        from_node="mc1_bra",
        to_node="mc1_bra_post",
        delta_geom_deg=90.0,
        area=A_bypass,
    )
    mpce_mc1 = MultiPortChamberElement(
        id="mpce_mc1",
        inlet_nodes=["mc1_com"],
        outlet_nodes=["mc1_str", "mc1_bra"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[A_main, A_main, A_bypass],
    )

    # mc2: merging junction (2 inflows from p3 and o_bypass, 1 outflow to p4).
    # Mirror layout: pre-loss node on the angled inlet feeds the bra port face.
    mc2_str = MomentumChamberNode("mc2_str", area=A_main)
    mc2_bra_pre = MomentumChamberNode("mc2_bra_pre", area=A_bypass)
    mc2_bra = MomentumChamberNode("mc2_bra", area=A_bypass)
    mc2_com = MomentumChamberNode("mc2_com", area=A_main)
    loss_mc2_bra = BorderCarnotLossElement(
        "loss_mc2_bra",
        from_node="mc2_bra_pre",
        to_node="mc2_bra",
        delta_geom_deg=90.0,
        area=A_bypass,
    )
    mpce_mc2 = MultiPortChamberElement(
        id="mpce_mc2",
        inlet_nodes=["mc2_str", "mc2_bra"],
        outlet_nodes=["mc2_com"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[A_main, A_bypass, A_main],
    )

    # Main series elements (now attach to port-face MCNs at mc1, mc2).
    p1 = ChannelElement("p1", "inlet", "n1", length=2.0, diameter=D_main, roughness=1e-4)
    o1 = OrificeElement("o1", "n1", "mc1_com", Cd=0.6, diameter=0.100925, correlation="fixed")
    p2 = ChannelElement("p2", "mc1_str", "n2", length=2.0, diameter=D_main, roughness=1e-4)
    o2 = OrificeElement("o2", "n2", "plenum", Cd=0.6, diameter=0.100925, correlation="fixed")
    p3 = ChannelElement("p3", "plenum", "mc2_str", length=2.0, diameter=D_main, roughness=1e-4)
    p4 = ChannelElement("p4", "mc2_com", "outlet", length=2.0, diameter=D_main, roughness=1e-4)

    # Bypass branch (same sizing as the original incarnation of the test).
    p_bypass = ChannelElement(
        "p_bypass", "mc1_bra_post", "n_bypass", length=5.0, diameter=D_bypass, roughness=1e-4
    )
    o_bypass = OrificeElement(
        "o_bypass", "n_bypass", "mc2_bra_pre", Cd=0.6, diameter=0.071365, correlation="fixed"
    )

    for n in [
        inlet,
        outlet,
        n1,
        n2,
        plen,
        n_bypass,
        mc1_com,
        mc1_str,
        mc1_bra,
        mc1_bra_post,
        mc2_str,
        mc2_bra_pre,
        mc2_bra,
        mc2_com,
    ]:
        graph.add_node(n)
    for e in [
        p1,
        o1,
        p2,
        o2,
        p3,
        p4,
        p_bypass,
        o_bypass,
        loss_mc1_bra,
        mpce_mc1,
        loss_mc2_bra,
        mpce_mc2,
    ]:
        graph.add_element(e)

    # Initial guesses on the named PlenumNodes (port-face MCNs default-init fine).
    n1.initial_guess = {"n1.P": 145000.0, "n1.Pt": 145000.0}
    n_bypass.initial_guess = {"n_bypass.P": 130000.0, "n_bypass.Pt": 130000.0}
    n2.initial_guess = {"n2.P": 135000.0, "n2.Pt": 135000.0}
    plen.initial_guess = {"plenum.P": 125000.0, "plenum.Pt": 125000.0}

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Mass conservation: split at mc1, rejoin at mc2.
    assert abs(sol["p1.m_dot"] - (sol["p2.m_dot"] + sol["p_bypass.m_dot"])) < 1e-6
    assert abs((sol["p3.m_dot"] + sol["p_bypass.m_dot"]) - sol["p4.m_dot"]) < 1e-6

    assert all(sol[f"{e}.m_dot"] > 0 for e in ["p1", "p2", "p3", "p4", "p_bypass"])

    # With the momentum-CV closure the diameter argument holds: smaller-bore
    # bypass (D=0.08) carries less flow than the main branch (D=0.1).
    assert sol["p_bypass.m_dot"] < sol["p2.m_dot"]


def test_network_validation_edge_cases():
    """Test additional edge cases for network validation"""

    # Test self-loop rejection
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=150000.0, Tt=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", Pt=100000.0, Tt=300.0, Y=_get_air_Y())

    graph.add_node(inlet)
    graph.add_node(outlet)

    # Should raise ValueError for self-loop
    with pytest.raises(ValueError, match="duplicate"):
        graph.add_element(
            OrificeElement(
                "self_loop", "inlet", "inlet", Cd=0.6, diameter=0.079788, correlation="fixed"
            )
        )

    # Test isolated node detection
    graph2 = FlowNetwork()
    inlet2 = PressureBoundary("inlet2", Pt=150000.0, Tt=300.0, Y=_get_air_Y())
    outlet2 = PressureBoundary("outlet2", Pt=100000.0, Tt=300.0, Y=_get_air_Y())
    isolated = PlenumNode("isolated")

    graph2.add_node(inlet2)
    graph2.add_node(outlet2)
    graph2.add_node(isolated)
    graph2.add_element(
        OrificeElement("conn", "inlet2", "outlet2", Cd=0.6, diameter=0.079788, correlation="fixed")
    )

    # Should raise ValueError for isolated node
    with pytest.raises(ValueError, match="isolated node"):
        graph2.validate()


def test_parallel_elements():
    """Test parallel elements between same nodes"""
    graph = FlowNetwork()
    inlet = PressureBoundary("inlet", Pt=150000.0, Tt=300.0, Y=_get_air_Y())
    outlet = PressureBoundary("outlet", Pt=100000.0, Tt=300.0, Y=_get_air_Y())

    graph.add_node(inlet)
    graph.add_node(outlet)

    # Two parallel orifices with different areas
    orf1 = OrificeElement(
        "orf_1", "inlet", "outlet", Cd=0.6, diameter=0.079788, correlation="fixed"
    )
    orf2 = OrificeElement(
        "orf_2", "inlet", "outlet", Cd=0.6, diameter=0.050463, correlation="fixed"
    )

    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve(method="lm")

    # Both should have positive flow
    assert sol["orf_1.m_dot"] > 0
    assert sol["orf_2.m_dot"] > 0

    # Larger orifice should carry more flow
    assert sol["orf_1.m_dot"] > sol["orf_2.m_dot"]

    # Total flow should be reasonable
    total_flow = sol["orf_1.m_dot"] + sol["orf_2.m_dot"]
    assert total_flow > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
