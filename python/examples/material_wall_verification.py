from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)


def verify_material_scaling() -> None:
    """Verifies that ThermalWall conductivity scales with temperature for database materials."""
    print("--- Verifying Material Conductivity Scaling ---")

    def solve_at_temp(T_hot: float) -> float:
        net = FlowNetwork()

        # Nodes
        net.add_node(PressureBoundary("hot_in", Pt=2e6, Tt=T_hot))
        net.add_node(PlenumNode("hot_p"))
        net.add_node(PressureBoundary("hot_out", Pt=1e6, Tt=T_hot))

        net.add_node(PressureBoundary("cold_in", Pt=2e6, Tt=300.0))
        net.add_node(PlenumNode("cold_p"))
        net.add_node(PressureBoundary("cold_out", Pt=1e6, Tt=300.0))

        # Elements (Channels with convective surfaces)
        # We use large area to ensure coupling is sensitive
        surf_area = 0.01
        net.add_element(
            ChannelElement(
                "hot_link",
                "hot_in",
                "hot_p",
                length=0.1,
                diameter=0.01,
                roughness=1e-5,
                surface=ConvectiveSurface(area=surf_area),
            )
        )
        net.add_element(
            ChannelElement(
                "hot_exit", "hot_p", "hot_out", length=0.1, diameter=0.01, roughness=1e-5
            )
        )

        net.add_element(
            ChannelElement(
                "cold_link",
                "cold_in",
                "cold_p",
                length=0.1,
                diameter=0.01,
                roughness=1e-5,
                surface=ConvectiveSurface(area=surf_area),
            )
        )
        net.add_element(
            ChannelElement(
                "cold_exit", "cold_p", "cold_out", length=0.1, diameter=0.01, roughness=1e-5
            )
        )

        # Thermal Wall with Inconel 718
        # Inconel 718 k increases with T (approx 11 W/mK at 300K, 25 W/mK at 1200K)
        wall = ThermalWall(
            id="test_wall",
            element_a="hot_link",
            element_b="cold_link",
            layers=[WallLayer(thickness=0.005, conductivity=11.0, material="inconel718")],
        )
        net.add_wall(wall)

        solver = NetworkSolver(net)
        # Solve with many steps to ensure lagged properties converge
        solver.solve(lambda_steps=30)

        # Get solve results
        k_final = wall.layers[0].conductivity
        print(f"Solved at T_hot={T_hot} K: Conductivity = {k_final:.3f} W/mK")
        return k_final

    # Solve at two points
    k_300 = solve_at_temp(300.0)
    k_1000 = solve_at_temp(1000.0)

    # Basic assertions - Inconel 718 k increases with T
    assert k_1000 > k_300, (
        f"Conductivity should increase with temperature for Inconel 718: k(300)={k_300:.3f}, k(1000)={k_1000:.3f}"
    )
    assert 10.0 <= k_300 <= 15.0, f"Expected k ~ 11.5 at 300K, got {k_300:.3f}"
    assert 18.0 <= k_1000 <= 30.0, f"Expected k ~ 20-25 with 1000K boundary, got {k_1000:.3f}"

    print("Success: Material conductivity scaled correctly!")


if __name__ == "__main__":
    verify_material_scaling()
