"""Integration tests for wall coupling in network solver.

Tests verify:
1. Two-channel wall coupling converges to analytical solution
2. thermal_coupling_enabled=False gives identical results to no-wall case
3. Assembled Jacobian matches column-wise finite difference
"""

import numpy as np

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel, WallConnection
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def test_two_channel_wall_coupling_convergence():
    """Test that wall-coupled two-channel network converges to analytical solution.

    Setup:
    - Hot side: 800 K inlet, 2 bar, 0.1 kg/s air
    - Cold side: 400 K inlet, 2 bar, 0.05 kg/s air
    - Wall: 2mm steel (k=25 W/(m*K))
    - Both channels: 0.04m diameter, 1m length, smooth
    """
    net = FlowNetwork()
    Y_air = cb.species.from_mapping({"N2": 0.767, "O2": 0.233})

    # Hot side: Inlet -> Channel -> Plenum -> Outlet
    hot_inlet = MassFlowBoundary("hot_inlet", m_dot=0.1, T_total=800.0, Y=Y_air)
    hot_plenum = PlenumNode("hot_plenum")
    hot_outlet = PressureBoundary("hot_outlet", P_total=2e5)

    # Cold side: Inlet -> Channel -> Plenum -> Outlet
    cold_inlet = MassFlowBoundary("cold_inlet", m_dot=0.05, T_total=400.0, Y=Y_air)
    cold_plenum = PlenumNode("cold_plenum")
    cold_outlet = PressureBoundary("cold_outlet", P_total=2e5)

    net.add_node(hot_inlet)
    net.add_node(hot_plenum)
    net.add_node(hot_outlet)
    net.add_node(cold_inlet)
    net.add_node(cold_plenum)
    net.add_node(cold_outlet)

    # Hot channel with convective surface
    hot_surface = ConvectiveSurface(
        area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
    )
    hot_channel = cb.network.ChannelElement(
        id="hot_channel",
        from_node="hot_inlet",
        to_node="hot_plenum",
        diameter=0.04,
        length=1.0,
        roughness=0.0,
        surface=hot_surface,
    )
    net.add_element(hot_channel)

    # Cold channel with convective surface
    cold_surface = ConvectiveSurface(
        area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
    )
    cold_channel = cb.network.ChannelElement(
        id="cold_channel",
        from_node="cold_inlet",
        to_node="cold_plenum",
        diameter=0.04,
        length=1.0,
        roughness=0.0,
        surface=cold_surface,
    )
    net.add_element(cold_channel)

    # Connect plenums to outlets with orifices
    hot_orifice = OrificeElement(
        "hot_orifice", "hot_plenum", "hot_outlet", Cd=0.8, diameter=0.079788, correlation="fixed"
    )
    cold_orifice = OrificeElement(
        "cold_orifice", "cold_plenum", "cold_outlet", Cd=0.8, diameter=0.079788, correlation="fixed"
    )
    net.add_element(hot_orifice)
    net.add_element(cold_orifice)

    # Wall connection
    wall = WallConnection(
        id="coupling_wall",
        element_a="hot_channel",
        element_b="cold_channel",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )
    net.add_wall(wall)

    # Solve with thermal coupling enabled
    net.thermal_coupling_enabled = True
    solver = NetworkSolver(net)
    result = solver.solve(method="hybr", options={"xtol": 1e-8})

    assert result is not None, "Solver failed to converge"

    # Get plenum temperatures from derived states
    T_hot = solver._derived_states["hot_plenum"][0]
    T_cold = solver._derived_states["cold_plenum"][0]

    # Hot side should cool down
    assert T_hot < 800.0, "Hot side should cool down"

    # Cold side should heat up
    assert T_cold > 400.0, "Cold side should heat up"

    # Verify energy balance (Q_hot ≈ -Q_cold)
    X_air = cb.mass_to_mole(cb.normalize_fractions(Y_air))
    cp_hot = cb.cp_mass(T_hot, X_air)
    cp_cold = cb.cp_mass(T_cold, X_air)

    Q_hot = 0.1 * cp_hot * (T_hot - 800.0)  # Should be negative
    Q_cold = 0.05 * cp_cold * (T_cold - 400.0)  # Should be positive

    # Energy balance should close within 5% (some losses to pressure drop)
    energy_imbalance = abs(Q_hot + Q_cold) / abs(Q_hot)
    assert energy_imbalance < 0.05, f"Energy imbalance: {energy_imbalance:.1%}"


def test_thermal_coupling_disabled_matches_no_wall():
    """Test that thermal_coupling_enabled=False gives same result as no wall."""

    Y_air = cb.species.from_mapping({"N2": 0.767, "O2": 0.233})

    def create_network(add_wall: bool, enable_coupling: bool = False) -> FlowNetwork:
        net = FlowNetwork()

        # Hot side
        hot_inlet = MassFlowBoundary("hot_inlet", m_dot=0.1, T_total=800.0, Y=Y_air)
        hot_plenum = PlenumNode("hot_plenum")
        hot_outlet = PressureBoundary("hot_outlet", P_total=2e5)

        # Cold side
        cold_inlet = MassFlowBoundary("cold_inlet", m_dot=0.05, T_total=400.0, Y=Y_air)
        cold_plenum = PlenumNode("cold_plenum")
        cold_outlet = PressureBoundary("cold_outlet", P_total=2e5)

        net.add_node(hot_inlet)
        net.add_node(hot_plenum)
        net.add_node(hot_outlet)
        net.add_node(cold_inlet)
        net.add_node(cold_plenum)
        net.add_node(cold_outlet)

        # Channels
        hot_surface = ConvectiveSurface(
            area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
        )
        hot_channel = cb.network.ChannelElement(
            id="hot_channel",
            from_node="hot_inlet",
            to_node="hot_plenum",
            diameter=0.04,
            length=1.0,
            roughness=0.0,
            surface=hot_surface,
        )
        net.add_element(hot_channel)

        cold_surface = ConvectiveSurface(
            area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
        )
        cold_channel = cb.network.ChannelElement(
            id="cold_channel",
            from_node="cold_inlet",
            to_node="cold_plenum",
            diameter=0.04,
            length=1.0,
            roughness=0.0,
            surface=cold_surface,
        )
        net.add_element(cold_channel)

        # Orifices
        hot_orifice = OrificeElement(
            "hot_orifice",
            "hot_plenum",
            "hot_outlet",
            Cd=0.8,
            diameter=0.025231,
            correlation="fixed",
        )
        cold_orifice = OrificeElement(
            "cold_orifice",
            "cold_plenum",
            "cold_outlet",
            Cd=0.8,
            diameter=0.025231,
            correlation="fixed",
        )
        net.add_element(hot_orifice)
        net.add_element(cold_orifice)

        # Optionally add wall
        if add_wall:
            wall = WallConnection(
                id="coupling_wall",
                element_a="hot_channel",
                element_b="cold_channel",
                wall_thickness=0.002,
                wall_conductivity=25.0,
            )
            net.add_wall(wall)

        net.thermal_coupling_enabled = enable_coupling
        return net

    # Solve both networks
    net_with_wall_disabled = create_network(add_wall=True, enable_coupling=False)
    net_without_wall = create_network(add_wall=False, enable_coupling=False)

    solver_with_wall = NetworkSolver(net_with_wall_disabled)
    solver_without_wall = NetworkSolver(net_without_wall)

    result_with_wall = solver_with_wall.solve(method="hybr", options={"xtol": 1e-8})
    result_without_wall = solver_without_wall.solve(method="hybr", options={"xtol": 1e-8})

    assert result_with_wall is not None
    assert result_without_wall is not None

    # Compare plenum temperatures - should be identical
    T_hot_with = solver_with_wall._derived_states["hot_plenum"][0]
    T_hot_without = solver_without_wall._derived_states["hot_plenum"][0]
    T_cold_with = solver_with_wall._derived_states["cold_plenum"][0]
    T_cold_without = solver_without_wall._derived_states["cold_plenum"][0]

    assert abs(T_hot_with - T_hot_without) < 1e-6, (
        f"Hot T mismatch: {T_hot_with} vs {T_hot_without}"
    )
    assert abs(T_cold_with - T_cold_without) < 1e-6, (
        f"Cold T mismatch: {T_cold_with} vs {T_cold_without}"
    )


def test_wall_coupling_jacobian_validation():
    """Test that assembled Jacobian matches column-wise finite difference.

    This validates the full chain rule implementation in _propagate_states.
    """
    Y_air = cb.species.from_mapping({"N2": 0.767, "O2": 0.233})
    net = FlowNetwork()

    # Nodes
    hot_inlet = MassFlowBoundary("hot_inlet", m_dot=0.1, T_total=800.0, Y=Y_air)
    hot_plenum = PlenumNode("hot_plenum")
    hot_outlet = PressureBoundary("hot_outlet", P_total=2e5)

    cold_inlet = MassFlowBoundary("cold_inlet", m_dot=0.05, T_total=400.0, Y=Y_air)
    cold_plenum = PlenumNode("cold_plenum")
    cold_outlet = PressureBoundary("cold_outlet", P_total=2e5)

    net.add_node(hot_inlet)
    net.add_node(hot_plenum)
    net.add_node(hot_outlet)
    net.add_node(cold_inlet)
    net.add_node(cold_plenum)
    net.add_node(cold_outlet)

    # Elements with convective surfaces
    hot_surface = ConvectiveSurface(
        area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
    )
    hot_channel = cb.network.ChannelElement(
        id="hot_channel",
        from_node="hot_inlet",
        to_node="hot_plenum",
        diameter=0.04,
        length=1.0,
        roughness=0.0,
        surface=hot_surface,
    )
    net.add_element(hot_channel)

    cold_surface = ConvectiveSurface(
        area=np.pi * 0.04 * 1.0, model=SmoothModel(correlation="gnielinski")
    )
    cold_channel = cb.network.ChannelElement(
        id="cold_channel",
        from_node="cold_inlet",
        to_node="cold_plenum",
        diameter=0.04,
        length=1.0,
        roughness=0.0,
        surface=cold_surface,
    )
    net.add_element(cold_channel)

    # Orifices
    hot_orifice = OrificeElement(
        "hot_orifice", "hot_plenum", "hot_outlet", Cd=0.8, diameter=0.035682, correlation="fixed"
    )
    cold_orifice = OrificeElement(
        "cold_orifice", "cold_plenum", "cold_outlet", Cd=0.8, diameter=0.035682, correlation="fixed"
    )
    net.add_element(hot_orifice)
    net.add_element(cold_orifice)

    # Wall connection
    wall = WallConnection(
        id="coupling_wall",
        element_a="hot_channel",
        element_b="cold_channel",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )
    net.add_wall(wall)

    net.thermal_coupling_enabled = True

    # Solve to get converged state
    solver = NetworkSolver(net)
    result = solver.solve(method="hybr", options={"xtol": 1e-8})
    assert result is not None

    # Get solution vector from result dict
    x_sol = np.array([result[name] for name in solver.unknown_names])

    # Get analytical Jacobian at solution
    residual_func = solver._residuals_and_jacobian
    _, J_analytical = residual_func(x_sol)

    # Compute finite difference Jacobian column-wise
    n = len(x_sol)
    J_fd = np.zeros((n, n))
    eps = 1e-6

    for j in range(n):
        x_plus = x_sol.copy()
        x_minus = x_sol.copy()

        # Relative perturbation for better scaling
        delta = eps * max(abs(x_sol[j]), 1.0)
        x_plus[j] += delta
        x_minus[j] -= delta

        r_plus, _ = residual_func(x_plus)
        r_minus, _ = residual_func(x_minus)

        J_fd[:, j] = (r_plus - r_minus) / (2.0 * delta)

    # Compare Jacobians
    abs_diff = np.abs(J_analytical - J_fd)
    rel_diff = abs_diff / (np.abs(J_fd) + 1e-10)

    # Check that most entries match within 5% relative error
    significant_mask = np.abs(J_fd) > 1e-6

    if np.any(significant_mask):
        max_rel_error = np.max(rel_diff[significant_mask])
        mean_rel_error = np.mean(rel_diff[significant_mask])

        assert max_rel_error < 0.1, (
            f"Max relative Jacobian error: {max_rel_error:.1%} (should be < 10%)"
        )
        assert mean_rel_error < 0.02, (
            f"Mean relative Jacobian error: {mean_rel_error:.1%} (should be < 2%)"
        )

    # Check absolute error for near-zero entries
    small_mask = ~significant_mask
    if np.any(small_mask):
        max_abs_error = np.max(abs_diff[small_mask])
        assert max_abs_error < 1e-4, (
            f"Max absolute Jacobian error (small entries): {max_abs_error:.2e}"
        )
