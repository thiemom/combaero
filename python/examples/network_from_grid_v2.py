import math
import sys
import time

import numpy as np

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)


def main():
    m_dot_air_input = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    print(f"--- Running Combustor Grid v2 with m_dot_air = {m_dot_air_input} kg/s ---")

    # Inlet composition and states
    T_air = 560.0
    P_inest = 500000.0  # Estimated inlet pressure for property calculations
    X_air = cb.humid_air_composition(288.16, 101325, 0.6)
    Y_air = cb.mole_to_mass(X_air)

    Y_fuel = [0.0] * cb.num_species()
    Y_fuel[cb.species_index_from_name("H2")] = 1.0
    X_fuel = cb.mass_to_mole(Y_fuel)

    # ---------------------------------------------------------
    # 1. Determine Required Fuel Flow for Phi = 1.0
    # ---------------------------------------------------------
    air_stream = cb.Stream()
    air_stream.set_T(T_air).set_P(P_inest).set_X(X_air).set_mdot(m_dot_air_input)

    fuel_stream = cb.Stream()
    fuel_stream.set_T(300.0).set_P(P_inest).set_X(X_fuel)

    fuel_stream = cb.set_fuel_stream_for_phi(1.0, fuel_stream, air_stream)
    m_dot_fuel = fuel_stream.mdot
    print(f"Targeting Phi=1.0. Calculated required fuel mass flow: {m_dot_fuel:.6f} kg/s")

    # ---------------------------------------------------------
    # 2. Build the Network
    # ---------------------------------------------------------
    net = FlowNetwork()

    # Boundaries
    inlet_air = MassFlowBoundary("inlet_air", m_dot=m_dot_air_input, T_total=T_air, Y=Y_air)
    inlet_fuel = MassFlowBoundary("inlet_fuel", m_dot=m_dot_fuel, T_total=300.0, Y=Y_fuel)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=Y_air)

    net.add_node(inlet_air)
    net.add_node(inlet_fuel)
    net.add_node(outlet)

    # Mixing Plenum & Combustor
    mixing_plenum = PlenumNode("mixing_plenum")
    combustor = CombustorNode("combustor", method="complete")

    net.add_node(mixing_plenum)
    net.add_node(combustor)

    # N Orifices from air to plenum (parallel branches to reduce high mach)
    n_air_orifices = 4
    for i in range(n_air_orifices):
        cd_val = 0.8
        d_orifice = 30.0 / 1000.0  # 30 mm
        A_orifice = 0.25 * np.pi * d_orifice**2
        net.add_element(
            OrificeElement(f"ori_air_{i}", "inlet_air", "mixing_plenum", Cd=cd_val, area=A_orifice)
        )

    # 1 Orifice for fuel
    d_fuel_ori = 10.0 / 1000.0
    A_fuel_ori = 0.25 * np.pi * d_fuel_ori**2
    net.add_element(
        OrificeElement("ori_fuel", "inlet_fuel", "mixing_plenum", Cd=0.6, area=A_fuel_ori)
    )

    # Orifice from plenum to combustor
    d_comb_in = 60.0 / 1000.0
    A_comb_in = 0.25 * np.pi * d_comb_in**2
    net.add_element(
        OrificeElement("ori_comb_in", "mixing_plenum", "combustor", Cd=0.8, area=A_comb_in)
    )

    # N Series Momentum Chambers (large area to keep Mach < 0.3)
    n_series = 5
    D_pipe = 0.4  # 40 cm diameter for very low Mach < 0.3
    L_pipe = 1.0
    A_pipe = 0.25 * np.pi * D_pipe**2

    chambers = []
    for i in range(n_series):
        chamber = MomentumChamberNode(f"chamber_{i}", area=A_pipe)
        chambers.append(chamber)
        net.add_node(chamber)

        from_node = "combustor" if i == 0 else f"chamber_{i - 1}"
        net.add_element(
            PipeElement(
                f"pipe_{i}",
                from_node,
                f"chamber_{i}",
                length=L_pipe,
                diameter=D_pipe,
                roughness=2e-5,
            )
        )

    # Final orifice to outlet
    net.add_element(
        OrificeElement("ori_exit", f"chamber_{n_series - 1}", "outlet", Cd=0.8, area=A_pipe)
    )

    print(f"Network built: {len(net.nodes)} nodes, {len(net.elements)} elements")

    # ---------------------------------------------------------
    # 3. Solve the Network
    # ---------------------------------------------------------
    solver_options = {
        "maxfev": 2000,
        "xtol": 1e-4,
        "ftol": 1e-8,
    }

    solver = NetworkSolver(net)
    start_time = time.time()

    # Use the robust multipass solver mechanism
    result = solver.solve(
        method="lm",
        timeout=120,
        options=solver_options,
        robust=True,
        init_strategy="incompressible_warmstart",
    )

    wall_time = time.time() - start_time
    print(f"Wall time: {wall_time:.2f} seconds")
    print(f"Solver success: {result.get('__success__', False)}")
    print(f"Message: {result.get('__message__', '')}")

    # ---------------------------------------------------------
    # 4. Diagnostics
    # ---------------------------------------------------------
    print("\n=== Diagnostics ===")
    print(f"Air Inlet P_total: {result.get('inlet_air.P_total', float('nan')) / 1e5:.3f} bar")
    print(f"Fuel Inlet P_total: {result.get('inlet_fuel.P_total', float('nan')) / 1e5:.3f} bar")
    print(f"Plenum P_total: {result.get('mixing_plenum.P_total', float('nan')) / 1e5:.3f} bar")
    print(f"Combustor T_total (T_ad): {result.get('combustor.T_total', float('nan')):.1f} K")
    print(f"Outlet P_total: {result.get('outlet.P_total', float('nan')) / 1e5:.3f} bar")

    mach_numbers = []
    for i in range(n_series):
        # Calculate mach in pipes
        m_dot = result.get(f"pipe_{i}.m_dot", float("nan"))
        from_node = "combustor" if i == 0 else f"chamber_{i - 1}"
        from_node = "combustor" if i == 0 else f"chamber_{i - 1}"

        # We need composition for Mach. The network solver resolves this. We can use extract_complete_states
        pass

    complete_states = solver.extract_complete_states(result)
    for i in range(n_series):
        from_node = "combustor" if i == 0 else f"chamber_{i - 1}"
        state = complete_states.get(from_node)
        m_dot = result.get(f"pipe_{i}.m_dot", float("nan"))
        if state and not math.isnan(m_dot):
            v = m_dot / (state.thermo.rho * A_pipe)
            # Reconstruct X for mach calculation
            # Actually, standard mach is `state.thermo.a` in the complete state!
            a = state.thermo.a
            mach = v / a
            mach_numbers.append(mach)
            print(f"Debug: pipe_{i}.m_dot = {m_dot:.3f} kg/s, Mach = {mach:.4f}")

    if mach_numbers:
        print(
            f"\nMach in downstream pipes - Min: {min(mach_numbers):.4f}, Mean: {sum(mach_numbers) / len(mach_numbers):.4f}, Max: {max(mach_numbers):.4f}"
        )


if __name__ == "__main__":
    main()
