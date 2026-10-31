"""Combustor pressure loss that grows with heat release.

Pressure loss is not a property of the combustor node: it sits on the edge
next to it, a PressureLossElement whose correlation reads the combustor's
temperature-rise ratio Theta = T_burned / T_unburned - 1. Here the liner loses
a total-pressure fraction dPt / Pt_in = xi0 + k * Theta (cold-flow 3%, plus
0.5% per unit Theta), and the example checks the solved loss against that
formula at the solved Theta.
"""

import combaero as cb
from combaero.network import (
    CombustorNode,
    FlowNetwork,
    LinearThetaFractionLoss,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
)

XI0 = 0.03  # cold-flow loss fraction (Theta = 0)
K_THETA = 0.005  # added loss fraction per unit Theta


def main() -> None:
    net = FlowNetwork()

    # Air: 10 kg/s at 700 K; fuel: CH4 at 300 K. The outlet is at 4.5 bar.
    net.add_node(MassFlowBoundary("Inlet", m_dot=10.0, Tt=700.0, Y=list(cb.species.dry_air_mass())))
    net.add_node(
        MassFlowBoundary("Fuel", m_dot=0.5, Tt=300.0, Y=list(cb.species.pure_species("CH4")))
    )
    net.add_node(CombustorNode("MainBurner"))
    net.add_node(PlenumNode("PreOutlet"))
    net.add_node(PressureBoundary("Outlet", Pt=450000.0))

    net.add_element(LosslessConnectionElement("Air", from_node="Inlet", to_node="MainBurner"))
    net.add_element(LosslessConnectionElement("FuelLine", from_node="Fuel", to_node="MainBurner"))
    # The liner loss, on the combustor's outlet edge.
    net.add_element(
        PressureLossElement(
            "Liner",
            from_node="MainBurner",
            to_node="PreOutlet",
            correlation=LinearThetaFractionLoss(k=K_THETA, xi0=XI0),
        )
    )
    net.add_element(
        OrificeElement(
            "Exit",
            from_node="PreOutlet",
            to_node="Outlet",
            Cd=0.8,
            correlation="fixed",
            diameter=0.356825,
        )
    )

    solver = NetworkSolver(net)
    res = solver.solve()
    if not res["__success__"]:
        raise SystemExit(f"solve failed: {res['__message__']}")

    burner = net.nodes["MainBurner"]
    theta = burner._T_burned / burner._T_unburned - 1.0
    pt_burner = res["MainBurner.Pt"]
    pt_after = res["PreOutlet.Pt"]
    xi_solved = 1.0 - pt_after / pt_burner
    xi_formula = XI0 + K_THETA * theta

    print("Combustor pressure loss linear in Theta")
    print(f"  T_unburned {burner._T_unburned:7.1f} K, T_burned {burner._T_burned:7.1f} K")
    print(f"  Theta = {theta:.3f}")
    print(f"  Pt: burner {pt_burner / 1e5:.4f} bar -> after liner {pt_after / 1e5:.4f} bar")
    print(f"  loss fraction solved {xi_solved:.5f}, formula xi0 + k Theta = {xi_formula:.5f}")
    if abs(xi_solved - xi_formula) > 1e-6:
        raise SystemExit("the solved loss does not match the correlation")


if __name__ == "__main__":
    main()
