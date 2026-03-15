import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    EffectiveAreaConnectionElement,
    MassFlowBoundary,
    OrificeElement,
    PressureBoundary,
)


def debug():
    net = FlowNetwork()
    Y_air = cb.standard_dry_air_composition()
    Y_ch4 = [0.0] * 14
    Y_ch4[5] = 1.0  # CH4 index

    m_dot_air = 2.0
    m_dot_fuel1 = 0.02
    m_dot_fuel2 = 0.04

    air_in = MassFlowBoundary("air", m_dot=m_dot_air, T_total=600.0, Y=Y_air)
    fuel1 = MassFlowBoundary("f1", m_dot=m_dot_fuel1, T_total=300.0, Y=Y_ch4)
    fuel2 = MassFlowBoundary("f2", m_dot=m_dot_fuel2, T_total=300.0, Y=Y_ch4)
    p_out = PressureBoundary("out", P_total=101325.0, T_total=300.0)
    combustor = CombustorNode("combustor", method="complete", pressure_loss_frac=0.05)

    net.add_node(air_in)
    net.add_node(fuel1)
    net.add_node(fuel2)
    net.add_node(combustor)
    net.add_node(p_out)

    e_air = EffectiveAreaConnectionElement("e_air", "air", "combustor", effective_area=0.05)
    e_f1 = EffectiveAreaConnectionElement("e_f1", "f1", "combustor", effective_area=0.005)
    e_f2 = EffectiveAreaConnectionElement("e_f2", "f2", "combustor", effective_area=0.005)
    nozzle = OrificeElement("nozzle", "combustor", "out", Cd=0.8, area=0.02)

    net.add_element(e_air)
    net.add_element(e_f1)
    net.add_element(e_f2)
    net.add_element(nozzle)

    solver = NetworkSolver(net)
    sol = solver.solve(method="hybr", use_jac=True)

    print(f"Success: {sol['__success__']}")
    print(f"Message: {sol['__message__']}")
    print(f"m_dot_air (inlet): {m_dot_air}")
    print(f"e_air.m_dot (sol): {sol['e_air.m_dot']}")
    print(f"nozzle.m_dot (sol): {sol['nozzle.m_dot']}")
    print(f"Inlet sum: {m_dot_air + m_dot_fuel1 + m_dot_fuel2}")

    # Check individual residuals
    x = solver._build_x0()
    # update x with sol
    for i, name in enumerate(solver.unknown_names):
        x[i] = sol[name]

    res = solver._residuals(x)
    print("\nResiduals:")
    for i, r in enumerate(res):
        if abs(r) > 1e-6:
            print(f"  [{i}] residual too high: {r}")


if __name__ == "__main__":
    debug()
