"""Integration tests to prove correct energy and mass conservation across network boundaries."""

import numpy as np
import pytest

import combaero as cb
from combaero.network import ConstantFractionLoss, FlowNetwork, NetworkSolver
from combaero.network.components import (
    CombustorNode,
    EffectiveAreaConnectionElement,
    MassFlowBoundary,
    OrificeElement,
    PressureBoundary,
)


def test_combustor_network_conservation():
    """Validates that a CombustorNode embedded in a network conserves mass and energy precisely."""
    net = FlowNetwork()

    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    m_dot_air = 2.0
    m_dot_fuel1 = 0.02
    m_dot_fuel2 = 0.04

    air_in = MassFlowBoundary("air", m_dot=m_dot_air, T_total=600.0, Y=Y_air)
    air_in.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel1 = MassFlowBoundary("f1", m_dot=m_dot_fuel1, T_total=300.0, Y=Y_ch4)
    fuel1.initial_guess = {"f1.P_total": 160000.0, "f1.P": 160000.0}
    fuel2 = MassFlowBoundary("f2", m_dot=m_dot_fuel2, T_total=300.0, Y=Y_ch4)
    fuel2.initial_guess = {"f2.P_total": 160000.0, "f2.P": 160000.0}

    p_out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    # complete combustion
    combustor = CombustorNode("combustor", method="complete")

    # Provide an initial guess to avoid dP=0 singularity across the nozzle at x0
    combustor.initial_guess = {
        "combustor.P_total": 150000.0,
        "combustor.P": 140000.0,
        "combustor.T": 1500.0,
    }

    net.add_node(air_in)
    net.add_node(fuel1)
    net.add_node(fuel2)
    net.add_node(combustor)
    net.add_node(p_out)

    # Direct connections to combustor (Add some area to drop pressure)
    e_air = EffectiveAreaConnectionElement("e_air", "air", "combustor", effective_area=0.05)
    e_air.initial_guess = {"e_air.m_dot": 2.0}

    e_f1 = EffectiveAreaConnectionElement("e_f1", "f1", "combustor", effective_area=0.005)
    e_f1.initial_guess = {"e_f1.m_dot": 0.02}

    e_f2 = EffectiveAreaConnectionElement("e_f2", "f2", "combustor", effective_area=0.005)
    e_f2.initial_guess = {"e_f2.m_dot": 0.04}

    net.add_element(e_air)
    net.add_element(e_f1)
    net.add_element(e_f2)

    # Discharge orifice
    nozzle = OrificeElement(
        "nozzle", "combustor", "out", Cd=0.8, diameter=0.159577, correlation="fixed"
    )
    nozzle.initial_guess = {"nozzle.m_dot": 2.06}
    net.add_element(nozzle)

    solver = NetworkSolver(net)
    sol = solver.solve(method="hybr", use_jac=True)

    # 1. Mass Conservation
    m_dot_out = sol["nozzle.m_dot"]
    m_dot_in_total = m_dot_air + m_dot_fuel1 + m_dot_fuel2

    assert np.isclose(m_dot_out, m_dot_in_total, rtol=1e-3), (
        f"Mass not conserved! In: {m_dot_in_total}, Out: {m_dot_out}"
    )

    # 2. Energy Conservation
    h_air, _ = cb._core.enthalpy_and_jacobian(600.0, cb.mass_to_mole(Y_air))
    h_ch4, _ = cb._core.enthalpy_and_jacobian(300.0, cb.mass_to_mole(Y_ch4))

    H_in = m_dot_air * h_air + (m_dot_fuel1 + m_dot_fuel2) * h_ch4

    T_out = sol["combustor.T"]
    X_out = [sol[f"combustor.Y[{i}]"] for i in range(14)]

    h_out, _ = cb._core.enthalpy_and_jacobian(T_out, cb.mass_to_mole(X_out))
    H_out = m_dot_out * h_out

    # Combustor calculates energy balance via combustion_residuals -> h_in_total = h_out_chamber
    # This proves the adiabatic solver effectively converged on energy
    assert np.isclose(H_in, H_out, rtol=1e-4), f"Energy not conserved! In: {H_in}, Out: {H_out}"

    # 3. Species Validity
    assert X_out[3] > 0.01, f"Expected CO2 in exhaust, got {X_out[3]}"
    assert X_out[4] > 0.01, f"Expected H2O in exhaust, got {X_out[4]}"
    assert np.isclose(sum(X_out), 1.0, rtol=1e-3), "Exhaust molar fractions don't sum to 1"


def _make_simple_combustor_net(xi: float | None):
    """Helper: simple network with mass-flow in -> combustor -> pressure-out."""
    net = FlowNetwork()
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    air_in = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air_in.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel_in = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel_in.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}
    p_out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    loss = ConstantFractionLoss(xi=xi) if xi is not None else None
    combustor = CombustorNode("comb", method="complete", pressure_loss=loss)
    combustor.initial_guess = {"comb.P_total": 150000.0, "comb.P": 140000.0}

    net.add_node(air_in)
    net.add_node(fuel_in)
    net.add_node(combustor)
    net.add_node(p_out)

    e_air = EffectiveAreaConnectionElement("e_air", "air", "comb", effective_area=0.05)
    e_air.initial_guess = {"e_air.m_dot": 1.0}
    e_fuel = EffectiveAreaConnectionElement("e_fuel", "fuel", "comb", effective_area=0.005)
    e_fuel.initial_guess = {"e_fuel.m_dot": 0.03}
    nozzle = OrificeElement("nozzle", "comb", "out", Cd=0.8, diameter=0.12, correlation="fixed")
    nozzle.initial_guess = {"nozzle.m_dot": 1.03}
    net.add_element(e_air)
    net.add_element(e_fuel)
    net.add_element(nozzle)

    return net, combustor


def test_pressure_loss_applied_in_network():
    """Regression: ConstantFractionLoss must apply exactly P_out = P_in * (1 - xi).

    The key invariant: combustor.P_total / upstream_P_total == (1 - xi).
    The upstream pressure floats (MassFlowBoundary); verify the ratio is exact.
    """
    xi = 0.05
    net_loss, combustor = _make_simple_combustor_net(xi=xi)
    sol = NetworkSolver(net_loss).solve(method="hybr", use_jac=True)

    mr = combustor._last_mix_res
    P_total_out = mr.P_total_mix

    # upstream P_total (MassFlowBoundary satisfies P_total = P, both floats together)
    P_in_air = sol["air.P"]  # air boundary P_total == P
    P_in_fuel = sol["fuel.P"]
    mdot_air, mdot_fuel = 1.0, 0.03
    P_in_combustor = (mdot_air * P_in_air + mdot_fuel * P_in_fuel) / (mdot_air + mdot_fuel)

    ratio = P_total_out / P_in_combustor
    assert ratio == pytest.approx(1 - xi, rel=1e-4), (
        f"Loss ratio wrong: P_out/P_in={ratio:.6f}, expected {1 - xi:.6f}"
    )
    assert P_total_out < P_in_combustor, "P_total must drop across combustor"
