"""Tests for MomentumChamberNode proper implementation."""

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)


def test_momentum_chamber_dynamic_pressure():
    """Verify MomentumChamberNode enforces P_total = P + 0.5*rho*v^2."""
    graph = FlowNetwork()

    Y_air = cb.species.dry_air_mass()
    inlet = MassFlowBoundary("inlet", m_dot=0.2, T_total=600.0, Y=Y_air)
    chamber = MomentumChamberNode(
        "chamber", area=0.02
    )  # 0.02 m^2, larger area for easier convergence
    outlet = PressureBoundary("outlet", P_total=101325.0)

    orf1 = OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682)
    orf2 = OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682)

    graph.add_node(inlet)
    graph.add_node(chamber)
    graph.add_node(outlet)
    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # Get chamber state
    P = sol["chamber.P"]
    P_total = sol["chamber.P_total"]
    m_dot = 0.2  # From inlet
    T = solver._derived_states["chamber"][0]
    Y = solver._derived_states["chamber"][1]

    # Verify dynamic pressure relation: P_total = P + 0.5*rho*v^2
    X = cb.mass_to_mole(cb.normalize_fractions(Y))
    rho = cb.density(T, P, X)
    v = m_dot / (rho * chamber.area)
    q_dynamic = 0.5 * rho * v * v

    # Check that P_total = P + q_dynamic
    assert P_total == pytest.approx(P + q_dynamic, rel=1e-4)

    # Verify dynamic pressure is non-zero (unlike plenum)
    assert q_dynamic > 1.0  # Should be several Pa

    print("Momentum chamber test passed:")
    print(f"  P = {P:.1f} Pa")
    print(f"  P_total = {P_total:.1f} Pa")
    print(f"  q_dynamic = {q_dynamic:.1f} Pa")
    print(f"  v = {v:.2f} m/s")


def test_momentum_chamber_vs_plenum():
    """Compare MomentumChamberNode to PlenumNode - should have different P_total."""
    from combaero.network import PlenumNode

    Y_air = cb.species.dry_air_mass()
    m_dot = 1.0  # Higher flow for larger dynamic pressure

    # Test with momentum chamber
    graph_mom = FlowNetwork()
    inlet_mom = MassFlowBoundary("inlet", m_dot=m_dot, T_total=600.0, Y=Y_air)
    chamber_mom = MomentumChamberNode("chamber", area=0.01)
    outlet_mom = PressureBoundary("outlet", P_total=101325.0)

    graph_mom.add_node(inlet_mom)
    graph_mom.add_node(chamber_mom)
    graph_mom.add_node(outlet_mom)
    graph_mom.add_element(OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682))
    graph_mom.add_element(OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682))

    solver_mom = NetworkSolver(graph_mom)
    sol_mom = solver_mom.solve()

    # Test with plenum (P_total = P)
    graph_ple = FlowNetwork()
    inlet_ple = MassFlowBoundary("inlet", m_dot=m_dot, T_total=600.0, Y=Y_air)
    chamber_ple = PlenumNode("chamber")
    outlet_ple = PressureBoundary("outlet", P_total=101325.0)

    graph_ple.add_node(inlet_ple)
    graph_ple.add_node(chamber_ple)
    graph_ple.add_node(outlet_ple)
    graph_ple.add_element(OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682))
    graph_ple.add_element(OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682))

    solver_ple = NetworkSolver(graph_ple)
    sol_ple = solver_ple.solve()

    assert sol_mom["__success__"]
    assert sol_ple["__success__"]

    # Plenum: P_total = P
    P_ple = sol_ple["chamber.P"]
    P_total_ple = sol_ple["chamber.P_total"]
    assert P_total_ple == pytest.approx(P_ple, rel=1e-6)

    # Momentum chamber: P_total > P (due to dynamic pressure)
    P_mom = sol_mom["chamber.P"]
    P_total_mom = sol_mom["chamber.P_total"]
    assert P_total_mom > P_mom

    # Dynamic pressure should be significant
    q_dynamic = P_total_mom - P_mom
    assert q_dynamic > 10.0  # At least 10 Pa

    print("Plenum vs Momentum Chamber:")
    print(
        f"  Plenum: P={P_ple:.1f} Pa, P_total={P_total_ple:.1f} Pa, ΔP={P_total_ple - P_ple:.1f} Pa"
    )
    print(f"  Momentum: P={P_mom:.1f} Pa, P_total={P_total_mom:.1f} Pa, ΔP={q_dynamic:.1f} Pa")


def test_momentum_chamber_with_heat_exchange():
    """Verify MomentumChamberNode supports energy boundaries (heat exchange)."""
    from combaero.network import EnergyBoundary

    graph = FlowNetwork()

    inlet = MassFlowBoundary("inlet", m_dot=0.5, T_total=600.0)
    inlet.Y = cb.species.dry_air_mass()
    chamber = MomentumChamberNode("chamber", area=0.01)
    outlet = PressureBoundary("outlet")
    outlet.P_total = 101325.0
    outlet.T_total = 300.0
    outlet.Y = cb.species.dry_air_mass()

    # Add heat exchange to chamber
    chamber.add_energy_boundary(EnergyBoundary("heater", Q=5000.0))

    orf1 = OrificeElement("orf1", "inlet", "chamber", Cd=0.6, diameter=0.035682)
    orf2 = OrificeElement("orf2", "chamber", "outlet", Cd=0.6, diameter=0.035682)

    graph.add_node(inlet)
    graph.add_node(chamber)
    graph.add_node(outlet)
    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # Verify temperature increased due to heat addition
    T_chamber = solver._derived_states["chamber"][0]
    T_inlet = 600.0
    assert T_chamber > T_inlet

    # Verify dynamic pressure relation still holds
    P = sol["chamber.P"]
    P_total = sol["chamber.P_total"]
    m_dot = 0.5
    Y = solver._derived_states["chamber"][1]

    X = cb.mass_to_mole(cb.normalize_fractions(Y))
    rho = cb.density(T_chamber, P, X)
    v = m_dot / (rho * chamber.area)
    q_dynamic = 0.5 * rho * v * v

    assert P_total == pytest.approx(P + q_dynamic, rel=1e-4)

    print("Momentum chamber with heat exchange:")
    print(
        f"  T_in = {T_inlet:.1f} K → T_chamber = {T_chamber:.1f} K (ΔT = {T_chamber - T_inlet:.1f} K)"
    )
    print(f"  P_total = P + q_dynamic: {P_total:.1f} = {P:.1f} + {q_dynamic:.1f} Pa")
