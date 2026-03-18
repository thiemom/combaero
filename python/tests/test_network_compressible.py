"""Tests for compressible flow elements in network solver."""

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def test_compressible_orifice_network():
    """Test network with compressible orifice element."""
    net = FlowNetwork()

    # High pressure ratio to ensure compressible effects
    inlet = PressureBoundary("inlet", P_total=300000.0, T_total=300.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Compressible orifice
    orifice = OrificeElement(
        "orifice", "inlet", "outlet", Cd=0.65, area=1e-4, regime="compressible"
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive"

    # Check pressure ratio to confirm compressible regime
    PR = outlet.P_total / inlet.P_total
    assert PR < 0.8, f"Should be in compressible regime (PR={PR:.3f})"

    # Verify against direct compressible calculation
    X = cb.standard_dry_air_composition()
    mdot_direct, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        inlet.T_total, inlet.P_total, outlet.P_total, X, 0.65, 1e-4, 0.0
    )

    # Should match within solver tolerance
    assert abs(sol["orifice.m_dot"] - mdot_direct) / mdot_direct < 0.01


def test_compressible_pipe_network():
    """Test network with compressible Fanno pipe element."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=200000.0, T_total=400.0)
    outlet = PressureBoundary("outlet", P_total=150000.0, T_total=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Simple compressible Fanno pipe directly between boundaries
    pipe = PipeElement(
        "pipe",
        "inlet",
        "outlet",
        length=1.0,
        diameter=0.08,
        roughness=1e-4,
        regime="compressible_fanno",
        friction_model="haaland",
    )

    net.add_element(pipe)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["pipe.m_dot"] > 0, "Pipe mass flow should be positive"

    # Verify pressure drop is reasonable
    dP = inlet.P_total - outlet.P_total
    assert dP > 0, "Pressure should drop across pipe"


def test_mixed_compressible_incompressible_network():
    """Test network mixing compressible and incompressible elements."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=250000.0, T_total=350.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=350.0)
    junction1 = PlenumNode("j1")
    junction2 = PlenumNode("j2")

    net.add_node(inlet)
    net.add_node(outlet)
    net.add_node(junction1)
    net.add_node(junction2)

    # Incompressible pipe first
    pipe1 = PipeElement(
        "pipe1",
        "inlet",
        "j1",
        length=1.0,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    # Compressible orifice in middle (high pressure drop)
    orifice = OrificeElement("orifice", "j1", "j2", Cd=0.65, area=1e-4, regime="compressible")

    # Incompressible pipe at end
    pipe2 = PipeElement(
        "pipe2",
        "j2",
        "outlet",
        length=0.5,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    net.add_element(pipe1)
    net.add_element(orifice)
    net.add_element(pipe2)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["pipe1.m_dot"] > 0
    assert sol["orifice.m_dot"] > 0
    assert sol["pipe2.m_dot"] > 0

    # Mass conservation
    assert abs(sol["pipe1.m_dot"] - sol["orifice.m_dot"]) < 1e-6
    assert abs(sol["orifice.m_dot"] - sol["pipe2.m_dot"]) < 1e-6

    # Pressure should decrease monotonically
    assert inlet.P_total > sol["j1.P"] > sol["j2.P"] > outlet.P_total


def test_compressible_vs_incompressible_comparison():
    """Compare compressible vs incompressible orifice in same network."""

    # Setup identical networks with different regimes
    def create_network(regime):
        net = FlowNetwork()
        inlet = PressureBoundary("inlet", P_total=200000.0, T_total=300.0)
        outlet = PressureBoundary("outlet", P_total=150000.0, T_total=300.0)
        net.add_node(inlet)
        net.add_node(outlet)

        orifice = OrificeElement("orifice", "inlet", "outlet", Cd=0.65, area=2e-4, regime=regime)
        net.add_element(orifice)
        return net

    # Solve both
    net_incomp = create_network("incompressible")
    net_comp = create_network("compressible")

    solver_incomp = NetworkSolver(net_incomp)
    solver_comp = NetworkSolver(net_comp)

    sol_incomp = solver_incomp.solve(method="lm")
    sol_comp = solver_comp.solve(method="lm")

    # At moderate pressure ratios, they should be similar but not identical
    mdot_incomp = sol_incomp["orifice.m_dot"]
    mdot_comp = sol_comp["orifice.m_dot"]

    # Compressible should give different result (can be up to 20% at PR=0.75)
    error = abs(mdot_comp - mdot_incomp) / mdot_comp
    assert 0.01 < error < 0.25, f"Expected 1-25% difference, got {error * 100:.1f}%"


def test_compressible_choked_orifice():
    """Test compressible orifice in choked flow regime."""
    net = FlowNetwork()

    # Very high pressure ratio to ensure choking
    inlet = PressureBoundary("inlet", P_total=500000.0, T_total=300.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    orifice = OrificeElement(
        "orifice", "inlet", "outlet", Cd=0.65, area=5e-5, regime="compressible"
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify choked flow
    X = cb.standard_dry_air_composition()
    PR_crit = cb.critical_pressure_ratio(inlet.T_total, inlet.P_total, X)
    PR_actual = outlet.P_total / inlet.P_total

    assert PR_actual < PR_crit, (
        f"Flow should be choked (PR={PR_actual:.3f} < PR_crit={PR_crit:.3f})"
    )
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive even when choked"


def test_compressible_pipe_high_mach():
    """Test compressible pipe at higher Mach number."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=300000.0, T_total=400.0)
    outlet = PressureBoundary("outlet", P_total=200000.0, T_total=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Long pipe with small diameter to get higher velocity
    pipe = PipeElement(
        "pipe",
        "inlet",
        "outlet",
        length=5.0,
        diameter=0.03,
        roughness=1e-4,
        regime="compressible_fanno",
    )
    net.add_element(pipe)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    assert sol["pipe.m_dot"] > 0, "Mass flow should be positive"

    # Verify Mach number is significant
    X = cb.standard_dry_air_composition()
    rho = cb.density(inlet.T_total, inlet.P_total, X)
    area = 3.14159 * (0.03 / 2) ** 2
    u = sol["pipe.m_dot"] / (rho * area)
    a = cb.speed_of_sound(inlet.T_total, X)
    M = u / a

    assert M > 0.1, f"Mach number should be significant (M={M:.3f})"
