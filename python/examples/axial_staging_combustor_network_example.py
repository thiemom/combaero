import numpy as np

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PressureBoundary,
)


def main():
    print("==============================================")
    print("CombAero - Multi-Stream Combustor Network Demo")
    print("==============================================")

    # Axial Staging Demo:
    # 1. Pilot stage (Air + Pilot Fuel)
    # 2. Main stage (Pilot Products + Air + Main Fuel)
    # 3. Axial stage (First Stage Products + Axial Fuel)
    net = FlowNetwork()

    # --- Boundary Parameters ---
    T_combustor_target = 1700.0  # [K]
    P_turbine_inlet = 20e5  # [Pa]

    # 1. Define Air Stream
    X_air = cb.species.humid_air(288.15, 101325.0, 0.6)
    air = cb.Stream()
    air.set_T(700.0).set_X(X_air).set_mdot(2.0)

    # 2. Define Fuel Stream and calculate stoichiometry
    X_fuel = cb.species.pure_species("CH4")
    fuel = cb.Stream()
    fuel.set_T(298.15).set_X(X_fuel)
    # Determine total fuel required to hit T_target for the given air stream
    fuel.set_mdot(cb.set_fuel_stream_for_Tad(T_combustor_target, fuel, air).mdot)

    # 3. Fuel Splitting
    AFF = 0.2  # Axial Fuel Fraction (of total fuel)
    PFF = 0.33  # Pilot Fuel Fraction (of the remaining 80%)

    m_dot_pilot = PFF * (1.0 - AFF) * fuel.mdot
    m_dot_main = (1.0 - PFF) * (1.0 - AFF) * fuel.mdot
    m_dot_axial = AFF * fuel.mdot

    # 4. Air Splitting
    m_dot_air_pilot = 0.3 * air.mdot
    m_dot_air_main = 0.7 * air.mdot

    print(f"Total Air:   {air.mdot:.3f} kg/s")
    print(f"Total Fuel:  {fuel.mdot:.4f} kg/s")
    print(f"  Pilot:     {m_dot_pilot:.4f} kg/s")
    print(f"  Main:      {m_dot_main:.4f} kg/s")
    print(f"  Axial:     {m_dot_axial:.4f} kg/s")

    # --- Nodes ---
    # Boundaries
    Y_air = air.state.Y
    Y_fuel = fuel.state.Y

    air_pilot_in = MassFlowBoundary("air_pilot_in", m_dot=m_dot_air_pilot, T_total=700.0, Y=Y_air)
    air_main_in = MassFlowBoundary("air_main_in", m_dot=m_dot_air_main, T_total=700.0, Y=Y_air)

    f_pilot_in = MassFlowBoundary("f_pilot_in", m_dot=m_dot_pilot, T_total=300.0, Y=Y_fuel)
    f_main_in = MassFlowBoundary("f_main_in", m_dot=m_dot_main, T_total=300.0, Y=Y_fuel)
    f_axial_in = MassFlowBoundary("f_axial_in", m_dot=m_dot_axial, T_total=300.0, Y=Y_fuel)

    turbine_back_pressure = PressureBoundary(
        "turbine_bc", P_total=P_turbine_inlet, T_total=T_combustor_target, Y=Y_air
    )

    # Internal Mixing/Combustion Nodes
    pilot_zone = CombustorNode("pilot_zone", method="complete")
    main_zone = CombustorNode("main_zone", method="complete")
    axial_zone = CombustorNode("axial_zone", method="complete")

    # Transition node to preserve momentum/velocity between stages
    intermediate_exit = MomentumChamberNode("intermediate_exit", area=0.05)

    # --- Elements ---
    # Stage 1: Pilot
    e_air_p = LosslessConnectionElement("e_air_p", "air_pilot_in", "pilot_zone")
    e_fuel_p = LosslessConnectionElement("e_fuel_p", "f_pilot_in", "pilot_zone")

    # Stage 1b: Main
    e_pilot_to_main = LosslessConnectionElement("e_pilot_to_main", "pilot_zone", "main_zone")
    e_air_m = LosslessConnectionElement("e_air_m", "air_main_in", "main_zone")
    e_fuel_m = LosslessConnectionElement("e_fuel_m", "f_main_in", "main_zone")

    # Transition
    e_main_to_intermediate = LosslessConnectionElement(
        "e_main_to_intermediate", "main_zone", "intermediate_exit"
    )

    # Stage 2: Axial
    e_inter_to_axial = LosslessConnectionElement(
        "e_inter_to_axial", "intermediate_exit", "axial_zone"
    )
    e_fuel_a = LosslessConnectionElement("e_fuel_a", "f_axial_in", "axial_zone")

    # Final Discharge
    nozzle = OrificeElement(
        "nozzle", "axial_zone", "turbine_bc", Cd=0.85, area=0.02, regime="compressible"
    )

    # --- Assemble Network ---
    net.add_node(air_pilot_in)
    net.add_node(air_main_in)
    net.add_node(f_pilot_in)
    net.add_node(f_main_in)
    net.add_node(f_axial_in)
    net.add_node(turbine_back_pressure)
    net.add_node(pilot_zone)
    net.add_node(main_zone)
    net.add_node(intermediate_exit)
    net.add_node(axial_zone)

    net.add_element(e_air_p)
    net.add_element(e_fuel_p)
    net.add_element(e_pilot_to_main)
    net.add_element(e_air_m)
    net.add_element(e_fuel_m)
    net.add_element(e_main_to_intermediate)
    net.add_element(e_inter_to_axial)
    net.add_element(e_fuel_a)
    net.add_element(nozzle)

    # --- Solve ---
    solver = NetworkSolver(net)
    print("\nSolving network (no initial guesses)...")
    try:
        sol = solver.solve(method="lm", use_jac=True)
    except Exception as e:
        print(f"Solver failed: {e}")
        return

    # --- Results ---
    print("\n" + "=" * 50)
    print("AXIAL STAGING COMPONENT RESULTS")
    print("=" * 50)
    print(f"{'Node':<25} {'T [K]':>10} {'P_tot [kPa]':>12}")
    print("-" * 50)
    for node_id in ["pilot_zone", "main_zone", "intermediate_exit", "axial_zone"]:
        T = sol[f"{node_id}.T"]
        P = sol[f"{node_id}.P_total"] / 1e3
        print(f"{node_id:<25} {T:10.2f} {P:12.2f}")

    print("\n--- Mass Balance ---")
    m_in_total = air.mdot + fuel.mdot
    m_out = sol["nozzle.m_dot"]
    print(f"Total Mass Flow In:  {m_in_total:10.6f} kg/s")
    print(f"Total Mass Flow Out: {m_out:10.6f} kg/s")
    print(f"Error:               {(m_in_total - m_out):10.2e} kg/s")

    print("\n--- Composition at Exit ---")
    mole_indices = {name: i for i, name in enumerate(cb.species.names)}
    important_species = ["CO2", "H2O", "O2", "N2", "CO", "H2"]

    Y_out = [sol[f"axial_zone.Y[{i}]"] for i in range(cb.species.num_species)]
    X_out = cb.species.to_mole(np.array(Y_out))

    for spec in important_species:
        idx = mole_indices[spec]
        print(f"  {spec:<4}: {X_out[idx] * 100:8.4f}%")

    print("\nDone.")


if __name__ == "__main__":
    main()
