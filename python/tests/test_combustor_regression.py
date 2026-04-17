"""Regression tests for CombustorNode bugs discovered during solver debugging.

These tests guard against:
1. Mass balance / solver divergence when CombustorNode is fed by MassFlowBoundary
   inlets via LosslessConnectionElement (the original symptom was comb.P being
   unconstrained, leading to divergence and a ~2% mass imbalance).
2. Pressure loss being silently ignored when using LosslessConnectionElement
   inlets (marked xfail until the structural fix for the P / P_total decoupling
   is implemented).
"""

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConstantFractionLoss,
    ConstantHeadLoss,
    FlowNetwork,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    LosslessConnectionElement,
    NetworkSolver,
    OrificeElement,
)
from combaero.network.components import (
    CombustorNode,
    EffectiveAreaConnectionElement,
    MassFlowBoundary,
    PressureBoundary,
)


def _make_lossless_inlet_network(xi=None):
    """Build network: MassFlowBoundary x2 -> Lossless -> Combustor -> Channel -> PressureBoundary."""
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    net = FlowNetwork()
    air = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}

    loss = ConstantFractionLoss(xi=xi) if xi is not None else None
    comb = CombustorNode("comb", method="complete", pressure_loss=loss)
    comb.initial_guess = {"comb.P_total": 150000.0, "comb.P": 150000.0}

    out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    for n in (air, fuel, comb, out):
        net.add_node(n)

    link_air = LosslessConnectionElement("l_air", "air", "comb")
    link_air.initial_guess = {"l_air.m_dot": 1.0}
    link_fuel = LosslessConnectionElement("l_fuel", "fuel", "comb")
    link_fuel.initial_guess = {"l_fuel.m_dot": 0.03}
    nozzle = OrificeElement("nozzle", "comb", "out", Cd=0.8, diameter=0.12, correlation="fixed")
    nozzle.initial_guess = {"nozzle.m_dot": 1.03}
    for e in (link_air, link_fuel, nozzle):
        net.add_element(e)

    return net, comb


def test_combustor_lossless_inlets_mass_balance():
    """Regression: CombustorNode fed by MassFlowBoundary via LosslessConnectionElement.

    Previously the comb.P unknown was unconstrained when the node residual was
    P_total - P_total_mix = 0, causing solver divergence and a ~2% mass imbalance.
    The fix (P_total - P = 0) ties both pressure unknowns together.
    """
    net, _comb = _make_lossless_inlet_network(xi=None)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)

    assert sol["__success__"], f"Solver failed to converge: |F|={sol['__final_norm__']:.2e}"
    assert sol["__final_norm__"] < 1e-6

    # Mass balance: l_air.m_dot + l_fuel.m_dot == nozzle.m_dot == 1.03
    m_in = sol["l_air.m_dot"] + sol["l_fuel.m_dot"]
    m_out = sol["nozzle.m_dot"]
    assert abs(m_in - 1.03) < 1e-6, f"Inlet mass wrong: {m_in:.6f} != 1.03"
    assert abs(m_out - 1.03) < 1e-6, f"Outlet mass wrong: {m_out:.6f} != 1.03"
    assert abs(m_in - m_out) < 1e-8, f"Mass imbalance: {m_in - m_out:.2e}"

    # P and P_total must be consistent (combustor is low-velocity volume).
    assert abs(sol["comb.P"] - sol["comb.P_total"]) < 1.0


def _make_effective_area_inlet_network(pressure_loss=None):
    """Build network: MassFlowBoundary x2 -> EffectiveArea -> Combustor -> Orifice -> PressureBoundary."""
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    net = FlowNetwork()
    air = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}
    comb = CombustorNode("comb", method="complete", pressure_loss=pressure_loss)
    comb.initial_guess = {"comb.P_total": 150000.0, "comb.P": 140000.0}
    out = PressureBoundary("out", P_total=101325.0, T_total=300.0)
    for n in (air, fuel, comb, out):
        net.add_node(n)

    e_air = EffectiveAreaConnectionElement("e_air", "air", "comb", effective_area=0.05)
    e_air.initial_guess = {"e_air.m_dot": 1.0}
    e_fuel = EffectiveAreaConnectionElement("e_fuel", "fuel", "comb", effective_area=0.005)
    e_fuel.initial_guess = {"e_fuel.m_dot": 0.03}
    nozzle = OrificeElement("nozzle", "comb", "out", Cd=0.8, diameter=0.12, correlation="fixed")
    nozzle.initial_guess = {"nozzle.m_dot": 1.03}
    for e in (e_air, e_fuel, nozzle):
        net.add_element(e)

    return net, comb


def test_combustor_pressure_loss_computes_mix_result():
    """Sanity check: each loss model produces a non-trivial P_total_mix in mix_res.

    This does NOT prove the loss affects the solution — just that the mixer is
    wired and produces the expected P_total_mix value. See the parametric
    ``test_combustor_pressure_loss_affects_solution`` for the enforcement check.
    """
    xi = 0.05
    net, comb = _make_effective_area_inlet_network(pressure_loss=ConstantFractionLoss(xi=xi))
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)
    assert sol["__success__"]
    mix_res = comb._last_mix_res
    P_in_weighted = (1.0 * sol["air.P"] + 0.03 * sol["fuel.P"]) / 1.03
    assert mix_res.P_total_mix == pytest.approx((1.0 - xi) * P_in_weighted, rel=1e-4)


LOSS_MODELS = [
    pytest.param(
        ConstantFractionLoss(xi=0.05), "ConstantFractionLoss(xi=0.05)", id="constant_fraction"
    ),
    pytest.param(
        LinearThetaFractionLoss(k=0.05, xi0=0.02),
        "LinearThetaFractionLoss(k=0.05, xi0=0.02)",
        id="linear_theta_fraction",
    ),
    pytest.param(
        ConstantHeadLoss(zeta=5.0, area=0.05), "ConstantHeadLoss(zeta=5.0)", id="constant_head"
    ),
    pytest.param(
        LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05),
        "LinearThetaHeadLoss(k=1.0, zeta0=3.0)",
        id="linear_theta_head",
    ),
]


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Known bug: pressure loss models compute P_total_mix correctly in the "
        "mixer but the CombustorNode residual (P_total - P = 0) never references "
        "P_total_mix — so the loss has zero effect on the Newton solution. "
        "Upstream boundary pressures are identical with or without loss. "
        "Affects all loss variants: ConstantFractionLoss, LinearThetaFractionLoss, "
        "ConstantHeadLoss, LinearThetaHeadLoss."
    ),
)
@pytest.mark.parametrize("loss_model,label", LOSS_MODELS)
def test_combustor_pressure_loss_affects_solution(loss_model, label):
    """The solution with loss MUST differ from the lossless baseline.

    Physically: with fixed outlet P_total and fixed inlet mass flows, introducing
    a combustor pressure loss must raise the upstream MassFlowBoundary P_total so
    the nozzle still passes the required flow against the same back pressure.
    """
    # Lossless baseline
    net0, _ = _make_effective_area_inlet_network(pressure_loss=None)
    sol0 = NetworkSolver(net0).solve(method="hybr", use_jac=True)
    assert sol0["__success__"], f"Lossless baseline failed: |F|={sol0['__final_norm__']:.2e}"

    # Lossy case
    net, _ = _make_effective_area_inlet_network(pressure_loss=loss_model)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)
    assert sol["__success__"], f"{label} solver failed: |F|={sol['__final_norm__']:.2e}"

    # Upstream pressure MUST rise to overcome the loss. Require at least 0.1% difference.
    rel_diff = abs(sol["air.P_total"] - sol0["air.P_total"]) / sol0["air.P_total"]
    assert rel_diff > 1e-3, (
        f"{label}: loss not enforced — air.P_total={sol['air.P_total']:.2f} vs "
        f"lossless={sol0['air.P_total']:.2f} (rel diff {rel_diff:.2e})"
    )


@pytest.mark.xfail(
    reason=(
        "Known bug: with LosslessConnectionElement inlets, the P_total - P = 0 residual "
        "leaves comb.P_total constrained only by upstream.P_total = comb.P_total (from "
        "the link), which has no dependency on P_total_mix. The pressure loss is computed "
        "but never enforced. Fix requires a second residual on P_total tied to P_total_mix, "
        "which currently over-constrains the system."
    ),
    strict=True,
)
def test_combustor_lossless_inlets_pressure_loss_enforced():
    """Regression target: pressure loss must be enforced even with lossless inlets.

    Expected physical behavior: comb.P_total = (1 - xi) * upstream.P_total, forcing
    Newton to raise the MassFlowBoundary pressures so the nozzle still passes 1.03 kg/s.
    """
    xi = 0.05
    net, comb = _make_lossless_inlet_network(xi=xi)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)
    assert sol["__success__"]

    # Upstream P_total must be higher than comb.P_total by the loss fraction.
    P_in_weighted = (1.0 * sol["air.P_total"] + 0.03 * sol["fuel.P_total"]) / 1.03
    ratio = sol["comb.P_total"] / P_in_weighted
    assert ratio == pytest.approx(1.0 - xi, rel=1e-4), (
        f"Loss not enforced: comb.P_total/P_in = {ratio:.5f}, expected {1.0 - xi:.5f}"
    )


def test_combustor_channel_outlet_network_converges():
    """Regression: the original failing network — MassFlowBoundary inlets, CombustorNode,
    ChannelElement, PressureBoundary outlet — must converge with correct mass balance."""
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    net = FlowNetwork()
    air = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}
    comb = CombustorNode("comb", method="complete")
    comb.initial_guess = {"comb.P_total": 150000.0, "comb.P": 150000.0}
    out = PressureBoundary("out", P_total=101325.0, T_total=300.0)
    for n in (air, fuel, comb, out):
        net.add_node(n)

    link_air = LosslessConnectionElement("l_air", "air", "comb")
    link_air.initial_guess = {"l_air.m_dot": 1.0}
    link_fuel = LosslessConnectionElement("l_fuel", "fuel", "comb")
    link_fuel.initial_guess = {"l_fuel.m_dot": 0.03}
    chan = ChannelElement("chan", "comb", "out", length=0.5, diameter=0.15, roughness=1e-5)
    chan.initial_guess = {"chan.m_dot": 1.03}
    for e in (link_air, link_fuel, chan):
        net.add_element(e)

    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)
    assert sol["__success__"], f"|F|={sol['__final_norm__']:.2e}"
    m_total = sol["l_air.m_dot"] + sol["l_fuel.m_dot"]
    assert abs(m_total - sol["chan.m_dot"]) < 1e-6
    assert abs(m_total - 1.03) < 1e-6
