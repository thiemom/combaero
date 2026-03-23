import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    LosslessConnectionElement,
    MassFlowBoundary,
    OrificeElement,
    PressureBoundary,
)


def main():
    print("==============================================")
    print("CombAero - Multi-Stream Combustor Network Demo")
    print("==============================================")

    # Example: 1 Air stream, 2 Fuel streams (pilot and main), discharging through an orifice
    net = FlowNetwork()

    # --- Boundaries ---
    # Air stream (standard air)
    X_air = cb.species.dry_air()
    Y_air = cb.species.dry_air_mass()
    m_dot_air = 2.0
    air_in = MassFlowBoundary("air_in", m_dot=m_dot_air, T_total=600.0, Y=Y_air)

    # Fuel/air ratio
    FAR = 0.03

    # Pilot fuel fraction
    PFF = 0.33

    # Fuel stream 1 (Pilot, pure CH4)
    X_ch4 = [0.0] * cb.num_species()
    X_ch4[cb.species_index_from_name("CH4")] = 1.0
    Y_ch4 = cb.mole_to_mass(X_ch4)
    m_dot_pilot = PFF * m_dot_air * FAR
    fuel_pilot = MassFlowBoundary("fuel_pilot", m_dot=m_dot_pilot, T_total=300.0, Y=Y_ch4)

    # Fuel stream 2 (Main, pure CH4)
    m_dot_main = (1 - PFF) * m_dot_air * FAR
    fuel_main = MassFlowBoundary("fuel_main", m_dot=m_dot_main, T_total=300.0, Y=Y_ch4)

    # Back pressure (discharge atmosphere)
    p_out = PressureBoundary("exhaust", P_total=101325.0, T_total=300.0, Y=Y_air)

    # --- Inner Nodes ---
    combustor = CombustorNode("combustor", method="complete")

    # Provide an initial guess to avoid dP=0 singularity across the nozzle at x0
    combustor.initial_guess = {
        "combustor.P_total": 150000.0,
        "combustor.P": 140000.0,
        "combustor.T": 1500.0,
    }

    # --- Elements ---
    # Connect everything directly to the combustor
    e_air = LosslessConnectionElement("e_air", "air_in", "combustor")
    e_air.initial_guess = {"e_air.m_dot": 2.0}

    e_pilot = LosslessConnectionElement("e_pilot", "fuel_pilot", "combustor")
    e_pilot.initial_guess = {"e_pilot.m_dot": 0.015}

    e_main = LosslessConnectionElement("e_main", "fuel_main", "combustor")
    e_main.initial_guess = {"e_main.m_dot": 0.045}

    # Connect to exhaust via nozzle
    nozzle = OrificeElement("nozzle", "combustor", "exhaust", Cd=0.8, area=0.03)
    nozzle.initial_guess = {"nozzle.m_dot": 2.06}

    net.add_node(air_in)
    net.add_node(fuel_pilot)
    net.add_node(fuel_main)
    net.add_node(combustor)
    net.add_node(p_out)

    net.add_element(e_air)
    net.add_element(e_pilot)
    net.add_element(e_main)
    net.add_element(nozzle)

    # Solve
    solver = NetworkSolver(net)
    print("Solving network...")
    sol = solver.solve(method="lm", use_jac=True)

    # Print results
    print("\n--- Combustor Results ---")
    print(f"P_total: {sol['combustor.P_total'] / 1000.0:.1f} kPa")
    print(f"T_out:   {sol['combustor.T']:.2f} K")

    print("\n--- Network Mass Balance ---")
    print(f"Mass Flow In:  {m_dot_air + m_dot_pilot + m_dot_main:.4f} kg/s")
    print(f"Mass Flow Out: {sol['nozzle.m_dot']:.4f} kg/s")

    print("\n--- Network Energy Balance ---")
    h_air, _ = cb._core.enthalpy_and_jacobian(600.0, X_air)
    h_ch4, _ = cb._core.enthalpy_and_jacobian(300.0, X_ch4)
    H_in = m_dot_air * h_air + (m_dot_pilot + m_dot_main) * h_ch4

    Y_out = [sol[f"combustor.Y[{i}]"] for i in range(cb.num_species())]
    X_out = cb.mass_to_mole(Y_out)
    h_out, _ = cb._core.enthalpy_and_jacobian(sol["combustor.T"], X_out)
    H_out = sol["nozzle.m_dot"] * h_out

    print(f"Energy Flux In:  {H_in / 1e6:.4f} MW")
    print(f"Energy Flux Out: {H_out / 1e6:.4f} MW")

    print(f"\nDischarge CO2 Mol Frac: {X_out[cb.species_index_from_name('CO2')]:.4f}")
    print(f"Discharge H2O Mol Frac: {X_out[cb.species_index_from_name('H2O')]:.4f}")


if __name__ == "__main__":
    main()
