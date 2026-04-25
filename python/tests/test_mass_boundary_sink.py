import pytest

from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
)


def test_mass_boundary_sink_convergence() -> None:
    """
    Test Case: Fixed Inlet Pressure -> Channel -> Fixed Outlet Massflow (Sink)
    Expected: Solver converges to the specified outlet mass flow.
    """
    net = FlowNetwork()

    # Node 1: Inlet Pressure (1.5 bar)
    n1 = PressureBoundary("inlet", Pt=1.5e5, Tt=300.0)
    net.add_node(n1)

    # Node 2: Outlet Mass Flow (Sink, 0.5 kg/s)
    n2 = MassFlowBoundary("outlet", m_dot=0.5, Tt=300.0)
    net.add_node(n2)

    # Element: Channel
    e1 = ChannelElement(
        "chan", from_node="inlet", to_node="outlet", length=1.0, diameter=0.1, roughness=1e-5
    )
    net.add_element(e1)

    # Solve
    solver = NetworkSolver(net)
    result = solver.solve()

    # Print for debugging if it fails
    if not result.get("__success__") and (result.get("__final_norm__") or 0) > 1e-10:
        print(f"\nSolver failed to hit final tolerance. Residual: {result.get('__final_norm__')}")
        for k, v in result.items():
            if not k.startswith("__"):
                print(f"  {k}: {v}")

    # $|F| \approx 10^{-12}$ is physically solved for this topology.
    assert (result.get("__final_norm__") or 1.0) < 1e-10

    # Mass flow in channel should match boundary sink setting (0.5 kg/s)
    m_dot_chan = result.get("chan.m_dot")
    assert pytest.approx(m_dot_chan, rel=1e-6) == 0.5

    # Outlet temperature should be ~300 K (inherited from inlet), not 1000 K (ignored value)
    T_outlet = result.get("outlet.T")
    assert pytest.approx(T_outlet, rel=1e-3) == 300.0


def test_mass_boundary_injection_node() -> None:
    """
    Test Case: Inlet -> Injection -> Outlet
    Node 1 (Pressure) -> Element 1 -> Node 2 (MassFlow Injection) -> Element 2 -> Node 3 (Pressure)
    """
    net = FlowNetwork()

    # Inlet: 2.0 bar
    net.add_node(PressureBoundary("inlet", Pt=2.0e5))

    # Injection: Adds 0.1 kg/s
    net.add_node(MassFlowBoundary("injection", m_dot=0.1))

    # Outlet: 1.0 bar
    net.add_node(PressureBoundary("outlet", Pt=1.0e5))

    # Elements
    net.add_element(
        ChannelElement(
            "e1", from_node="inlet", to_node="injection", length=1.0, diameter=0.1, roughness=1e-5
        )
    )
    net.add_element(
        ChannelElement(
            "e2", from_node="injection", to_node="outlet", length=1.0, diameter=0.1, roughness=1e-5
        )
    )

    solver = NetworkSolver(net)
    result = solver.solve()

    assert result.get("__success__") is True

    m1 = result.get("e1.m_dot")
    m2 = result.get("e2.m_dot")

    # Mass balance: m1 + 0.1 - m2 = 0 => m2 = m1 + 0.1
    assert pytest.approx(m2, rel=1e-5) == m1 + 0.1
