import combaero as cb
from combaero.network import (
    CombustorNode,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def test_combustor_pressure_loss():
    # 1. Setup Network
    net = FlowNetwork()

    # Inlet Air: 10 kg/s, 700 K, 5 bar
    Y_air = cb.species.dry_air_mass()
    inlet = MassFlowBoundary("Inlet", m_dot=10.0, T_total=700.0, Y=list(Y_air))
    net.add_node(inlet)

    # Fuel Inlet: CH4 at 300K
    Y_fuel = cb.species.pure_species("CH4")
    fuel_inlet = MassFlowBoundary("Fuel", m_dot=0.5, T_total=300.0, Y=list(Y_fuel))
    net.add_node(fuel_inlet)

    # Combustor: zeta = 0.03 + 0.005 * theta
    def p_loss(ctx):
        # Return (zeta, d_zeta_d_theta)
        return 0.03 + 0.005 * ctx.theta, 0.005

    combustor = CombustorNode("MainBurner", pressure_loss=p_loss)
    net.add_node(combustor)

    # Intermediate node
    pre_outlet = PlenumNode("PreOutlet")
    net.add_node(pre_outlet)

    # Outlet: 4.5 bar
    outlet = PressureBoundary("Outlet", P_total=450000.0)
    net.add_node(outlet)

    # Connections
    net.add_element(LosslessConnectionElement("Stream1", from_node="Inlet", to_node="MainBurner"))
    net.add_element(LosslessConnectionElement("Stream2", from_node="Fuel", to_node="MainBurner"))
    net.add_element(
        LosslessConnectionElement("Stream3", from_node="MainBurner", to_node="PreOutlet")
    )

    # Large area orifice
    orifice = OrificeElement("Orifice", from_node="PreOutlet", to_node="Outlet", Cd=0.8, area=0.1)
    net.add_element(orifice)

    solver = NetworkSolver(net)

    # 2. Check Squareness
    # Trigger unknown initialization
    x0_auto = solver._build_x0()
    f0, _ = solver._residuals_and_jacobian(x0_auto)

    print(f"Total Unknowns: {len(solver.unknown_names)}")
    print(f"Total Residuals: {len(f0)}")

    if len(solver.unknown_names) != len(f0):
        print("WARNING: Network is not square!")
        # Print equation counts
        # We know node residuals are collected first, then elements.

    # 3. Solve
    print("\nSolving network...")
    res = solver.solve()

    # 4. Verify
    state = solver._get_node_state(combustor, solver._last_x)
    print(f"Convergence successful? {res.get('__success__')}")
    print(f"Combustor Outlet T: {state.T:.2f} K")
    print(f"Combustor Outlet P: {state.P:.2f} Pa")

    # Residuals
    f_res, _ = solver._residuals_and_jacobian(solver._last_x)
    print("\nFinal residuals:")
    for i, val in enumerate(f_res):
        name = solver.unknown_names[i] if i < len(solver.unknown_names) else f"(Extra) eq {i}"
        print(f"  {name}: {val:.2e}")


if __name__ == "__main__":
    test_combustor_pressure_loss()
