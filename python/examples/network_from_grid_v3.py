import math
import time

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    ConvectiveSurface,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
    SmoothModel,
    WallConnection,
)


def build_and_solve_network(
    m_dot_air, T_ad, n_serial, n_parallel, total_area, total_ht_area, timeout=120, maxfev=10000
):
    T_air = 560.0
    P_ref = 101325.0  # Standard reference pressure for initial stream properties
    X_air = cb.humid_air_composition(288.16, P_ref, 0.6)
    Y_air = cb.mole_to_mass(X_air)

    Y_fuel = [0.0] * cb.num_species()
    Y_fuel[cb.species_index_from_name("H2")] = 1.0
    X_fuel = cb.mass_to_mole(Y_fuel)

    air_stream = cb.Stream()
    air_stream.set_T(T_air).set_P(P_ref).set_X(X_air).set_mdot(m_dot_air)

    fuel_stream = cb.Stream()
    fuel_stream.set_T(300.0).set_P(P_ref).set_X(X_fuel)

    fuel_stream = cb.set_fuel_stream_for_Tad(T_ad, fuel_stream, air_stream)
    m_dot_fuel = fuel_stream.mdot

    net = FlowNetwork()

    inlet_air = MassFlowBoundary("inlet_air", m_dot=m_dot_air, T_total=T_air, Y=Y_air)
    inlet_fuel = MassFlowBoundary("inlet_fuel", m_dot=m_dot_fuel, T_total=300.0, Y=Y_fuel)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=Y_air)

    net.add_node(inlet_air)
    net.add_node(inlet_fuel)
    net.add_node(outlet)

    distributor = PlenumNode("distributor")
    collector = PlenumNode("collector")
    net.add_node(distributor)
    net.add_node(collector)

    net.add_element(
        OrificeElement(
            "ori_air_in", "inlet_air", "distributor", Cd=0.8, area=total_area, regime="compressible"
        )
    )

    A_branch = total_area / n_parallel
    D_branch = math.sqrt(4 * A_branch / math.pi)

    A_ht_pipe = (total_ht_area / n_parallel) / max(1, n_serial)
    L_pipe = 1.0

    distributor_fuel = PlenumNode("distributor_fuel")
    net.add_node(distributor_fuel)
    net.add_element(
        OrificeElement(
            "ori_fuel_in",
            "inlet_fuel",
            "distributor_fuel",
            Cd=0.6,
            area=A_branch * 0.1 * n_parallel,
            regime="compressible",
        )
    )

    for p in range(n_parallel):
        mixing = PlenumNode(f"mixing_plenum_{p}")
        combustor = CombustorNode(f"combustor_{p}", method="complete")
        net.add_node(mixing)
        net.add_node(combustor)

        net.add_element(
            OrificeElement(
                f"ori_fuel_{p}",
                "distributor_fuel",
                f"mixing_plenum_{p}",
                Cd=0.6,
                area=A_branch * 0.1,
                regime="compressible",
            )
        )

        for i in range(n_serial + 1):
            net.add_node(MomentumChamberNode(f"cold_chamber_{p}_{i}", area=A_branch))
            net.add_node(MomentumChamberNode(f"hot_chamber_{p}_{i}", area=A_branch))

        net.add_element(
            OrificeElement(
                f"ori_cold_in_{p}",
                "distributor",
                f"cold_chamber_{p}_0",
                Cd=0.8,
                area=A_branch,
                regime="compressible",
            )
        )
        for i in range(n_serial):
            net.add_element(
                PipeElement(
                    f"cold_pipe_{p}_{i}",
                    f"cold_chamber_{p}_{i}",
                    f"cold_chamber_{p}_{i + 1}",
                    length=L_pipe,
                    diameter=D_branch,
                    roughness=2e-5,
                    regime="compressible_fanno",
                    surface=ConvectiveSurface(area=A_ht_pipe, model=SmoothModel())
                    if A_ht_pipe > 0
                    else None,
                )
            )

        net.add_element(
            OrificeElement(
                f"ori_comb_in_{p}",
                f"cold_chamber_{p}_{n_serial}",
                f"mixing_plenum_{p}",
                Cd=0.8,
                area=A_branch,
                regime="compressible",
            )
        )
        net.add_element(
            OrificeElement(
                f"ori_mix_to_comb_{p}",
                f"mixing_plenum_{p}",
                f"combustor_{p}",
                Cd=0.8,
                area=A_branch,
                regime="compressible",
            )
        )
        net.add_element(
            OrificeElement(
                f"ori_comb_out_{p}",
                f"combustor_{p}",
                f"hot_chamber_{p}_0",
                Cd=0.8,
                area=A_branch,
                regime="compressible",
            )
        )

        for i in range(n_serial):
            net.add_element(
                PipeElement(
                    f"hot_pipe_{p}_{i}",
                    f"hot_chamber_{p}_{i}",
                    f"hot_chamber_{p}_{i + 1}",
                    length=L_pipe,
                    diameter=D_branch,
                    roughness=2e-5,
                    regime="compressible_fanno",
                    surface=ConvectiveSurface(area=A_ht_pipe, model=SmoothModel())
                    if A_ht_pipe > 0
                    else None,
                )
            )

            if A_ht_pipe > 0:
                wall_cold_pipe_idx = n_serial - 1 - i
                net.add_wall(
                    WallConnection(
                        id=f"wall_{p}_{i}",
                        element_a=f"hot_pipe_{p}_{i}",
                        element_b=f"cold_pipe_{p}_{wall_cold_pipe_idx}",
                        wall_thickness=0.005,
                        wall_conductivity=20.0,
                        contact_area=A_ht_pipe,
                    )
                )

        net.add_element(
            OrificeElement(
                f"ori_hot_out_{p}",
                f"hot_chamber_{p}_{n_serial}",
                "collector",
                Cd=0.8,
                area=A_branch,
                regime="compressible",
            )
        )

    net.add_element(
        OrificeElement(
            "ori_exhaust", "collector", "outlet", Cd=0.8, area=total_area, regime="compressible"
        )
    )

    solver_options = {
        "maxfev": maxfev,
        "xtol": 1e-4,
        "ftol": 1e-6,
    }

    solver = NetworkSolver(net)
    start_time = time.time()
    result = solver.solve(
        method="lm",
        timeout=timeout,
        options=solver_options,
        robust=True,
        init_strategy="incompressible_warmstart",
    )
    wall_time = time.time() - start_time

    complete_states = solver.extract_complete_states(result)
    max_hot_mach = 0.0

    for p in range(n_parallel):
        for i in range(n_serial + 1):
            m_dot = result.get(f"ori_hot_out_{p}.m_dot", float("nan"))
            state = complete_states.get(f"hot_chamber_{p}_{i}")
            if state and not math.isnan(m_dot):
                v = m_dot / (state.thermo.rho * A_branch)
                mach = v / state.thermo.a
                max_hot_mach = max(max_hot_mach, mach)

    return {
        "success": result.get("__success__", False),
        "message": result.get("__message__", ""),
        "time": wall_time,
        "n_nodes": len(net.nodes),
        "n_elems": len(net.elements),
        "P_inlet_air": result.get("inlet_air.P_total", float("nan")) / 1e5,
        "T_comb": result.get("combustor_0.T_total", float("nan")),
        "T_hot_out": result.get("collector.T_total", float("nan")),
        "max_hot_mach": max_hot_mach,
        "residual": result.get("__f_norm__", float("nan")),
    }


def print_result(params, res):
    s = "PASS" if res["success"] else "FAIL"
    print(
        f"{s} | m={params['m_dot_air']:4.1f} nP={params['n_parallel']} nS={params['n_serial']} HT={params['total_ht_area']:4.0f} Tad={params['T_ad']:4.0f} | "
        + f"t={res['time']:5.1f}s P_in={res['P_inlet_air']:5.1f}b Tcmb={res['T_comb']:6.1f}K Thout={res['T_hot_out']:6.1f}K Mach={res['max_hot_mach']:4.2f}"
    )


def main():
    total_area = 0.25 * math.pi * (0.8) ** 2  # Even safer Mach numbers, max total diameter 0.8m
    tests = [
        {"m_dot_air": 10.0, "T_ad": 1700.0, "n_parallel": 1, "n_serial": 1, "total_ht_area": 50.0},
        {"m_dot_air": 10.0, "T_ad": 1700.0, "n_parallel": 2, "n_serial": 1, "total_ht_area": 50.0},
        {"m_dot_air": 10.0, "T_ad": 1700.0, "n_parallel": 1, "n_serial": 5, "total_ht_area": 50.0},
        {"m_dot_air": 10.0, "T_ad": 1700.0, "n_parallel": 2, "n_serial": 5, "total_ht_area": 50.0},
        {"m_dot_air": 20.0, "T_ad": 1700.0, "n_parallel": 2, "n_serial": 2, "total_ht_area": 50.0},
        {"m_dot_air": 20.0, "T_ad": 2200.0, "n_parallel": 2, "n_serial": 2, "total_ht_area": 50.0},
        {
            "m_dot_air": 10.0,
            "T_ad": 1700.0,
            "n_parallel": 1,
            "n_serial": 10,
            "total_ht_area": 100.0,
        },
        {"m_dot_air": 10.0, "T_ad": 1700.0, "n_parallel": 4, "n_serial": 4, "total_ht_area": 100.0},
        {
            "m_dot_air": 10.0,
            "T_ad": 1700.0,
            "n_parallel": 10,
            "n_serial": 10,
            "total_ht_area": 200.0,
        },
        {
            "m_dot_air": 10.0,
            "T_ad": 1700.0,
            "n_parallel": 20,
            "n_serial": 5,
            "total_ht_area": 200.0,
        },
        {
            "m_dot_air": 10.0,
            "T_ad": 1700.0,
            "n_parallel": 5,
            "n_serial": 20,
            "total_ht_area": 200.0,
        },
    ]

    print(f"--- Running Combustor Full Coupling Limits ({len(tests)} cases) ---")
    for params in tests:
        try:
            res = build_and_solve_network(**params, total_area=total_area, timeout=300)
            print_result(params, res)
            if not res["success"]:
                print(f"    FAIL REASON: {res['message']}")
        except Exception as e:
            print(
                f"ERROR | m={params['m_dot_air']:4.1f} nP={params['n_parallel']} nS={params['n_serial']} HT={params['total_ht_area']:4.0f} Tad={params['T_ad']:4.0f} | Exception: {str(e)}"
            )


if __name__ == "__main__":
    main()
