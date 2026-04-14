from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)


def test_thermal_wall_multilayer_compute():
    """Test multi-layer thermal wall coupling calculation."""
    # Side A: 500 W/m2K, 1000 K
    # Side B: 2000 W/m2K, 300 K
    # Layer 1: t=1mm, k=20 W/mK -> r=5e-5
    # Layer 2: t=5mm, k=0.1 W/mK -> r=0.05
    # Area: 1.0 m^2

    L1 = WallLayer(thickness=0.001, conductivity=20.0)
    L2 = WallLayer(thickness=0.005, conductivity=0.1)

    wall = ThermalWall(id="wall1", element_a="a", element_b="b", layers=[L1, L2], contact_area=1.0)

    # h_a=500, T_a=1000, h_b=2000, T_b=300
    res = wall.compute_coupling(500.0, 1000.0, 1.0, 2000.0, 300.0, 1.0)

    # R_total = 1/500 + 5e-5 + 0.05 + 1/2000 = 0.002 + 0.00005 + 0.05 + 0.0005 = 0.05255
    # Q = A * (T_a - T_b) / R_total = 1.0 * 700 / 0.05255 = 13320.65
    assert abs(res.Q - 13320.647) < 0.1

    # T_surf_a = T0 - Q/(h_a*A) = 1000 - 13320.6 / 500 = 1000 - 26.64 = 973.36
    assert abs(res.T_wall - 973.358) < 0.1


def test_thermal_wall_fouling():
    """Test fouling resistance influence."""
    L1 = WallLayer(thickness=0.001, conductivity=50.0)

    # No fouling
    wall0 = ThermalWall(id="w0", element_a="a", element_b="b", layers=[L1], contact_area=1.0)
    res0 = wall0.compute_coupling(1000.0, 1000.0, 1.0, 1000.0, 300.0, 1.0)

    # With fouling
    wall1 = ThermalWall(
        id="w1", element_a="a", element_b="b", layers=[L1], contact_area=1.0, R_fouling=0.01
    )
    res1 = wall1.compute_coupling(1000.0, 1000.0, 1.0, 1000.0, 300.0, 1.0)

    assert res1.Q < res0.Q
    # Q_fouled = 700 / (1/1000 + 0.001/50 + 0.01 + 1/1000) = 700 / 0.01202 = 58236.27
    assert abs(res1.Q - 58236.27) < 10.0


def test_solver_interface_temperature():
    """Test that solver populates T_interface correctly with real heat flow."""
    net = FlowNetwork()

    # Path A: Hot (1000K -> 950K)
    in_a = PressureBoundary("in_a", P_total=200000, T_total=1000)
    out_a = PressureBoundary("out_a", P_total=150000, T_total=1000)
    net.add_node(in_a)
    net.add_node(out_a)
    p_hot = ChannelElement(
        "p_hot", from_node="in_a", to_node="out_a", length=1.0, diameter=0.05, roughness=1e-5
    )
    p_hot.surface = ConvectiveSurface(area=0.1)
    net.add_element(p_hot)

    # Path B: Cold (300K -> 350K)
    in_b = PressureBoundary("in_b", P_total=200000, T_total=300)
    out_b = PressureBoundary("out_b", P_total=150000, T_total=300)
    net.add_node(in_b)
    net.add_node(out_b)
    p_cold = ChannelElement(
        "p_cold", from_node="in_b", to_node="out_b", length=1.0, diameter=0.05, roughness=1e-5
    )
    p_cold.surface = ConvectiveSurface(area=0.1)
    net.add_element(p_cold)

    # 3-layer wall: Hot -> L1 -> L2 -> L3 -> Cold
    L1 = WallLayer(thickness=0.001, conductivity=50.0)  # r=2e-5
    L2 = WallLayer(thickness=0.010, conductivity=0.1)  # r=0.1 (insulator)
    L3 = WallLayer(thickness=0.001, conductivity=50.0)  # r=2e-5

    wall = ThermalWall(
        id="wall", element_a="p_hot", element_b="p_cold", layers=[L1, L2, L3], contact_area=0.1
    )
    net.add_wall(wall)

    solver = NetworkSolver(net)
    res = solver.solve()

    assert res["__success__"]
    assert "wall.T_interface" in res
    profile = res["wall.T_interface"]

    # 3 layers -> 4 temperatures [T_surf_hot, T12, T23, T_surf_cold]
    assert len(profile) == 4
    # Check strict monotonicity (Hot > Cold)
    assert profile[0] > profile[1] > profile[2] > profile[3]
    # Check that T12 and T23 are between surface temperatures
    assert profile[0] > profile[1]
    assert profile[2] > profile[3]
    # Insulator should take the bulk of the dT
    dT_L2 = profile[1] - profile[2]
    dT_total = profile[0] - profile[3]
    assert dT_L2 / dT_total > 0.9  # R2 is much larger than R1, R3
