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
    PressureLossElement,
)
from combaero.network.components import (
    CombustorNode,
    EffectiveAreaConnectionElement,
    MassFlowBoundary,
    PlenumNode,
    PressureBoundary,
)


def _make_lossless_inlet_network(xi=None):
    """Build network: MassFlowBoundary x2 -> Lossless -> Combustor -> [loss?] -> Orifice -> PressureBoundary.

    When ``xi`` is provided, a :class:`PressureLossElement` with a
    :class:`ConstantFractionLoss` correlation sits between the combustor and
    a ``post`` plenum; the orifice is downstream of ``post``.  When ``xi`` is
    None the combustor connects to the orifice via a lossless link.
    """
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    net = FlowNetwork()
    air = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}

    comb = CombustorNode("comb", method="complete")
    comb.initial_guess = {"comb.P_total": 150000.0, "comb.P": 150000.0}
    post = PlenumNode("post")
    post.initial_guess = {"post.P_total": 140000.0, "post.P": 140000.0}

    out = PressureBoundary("out", P_total=101325.0, T_total=300.0)

    for n in (air, fuel, comb, post, out):
        net.add_node(n)

    link_air = LosslessConnectionElement("l_air", "air", "comb")
    link_air.initial_guess = {"l_air.m_dot": 1.0}
    link_fuel = LosslessConnectionElement("l_fuel", "fuel", "comb")
    link_fuel.initial_guess = {"l_fuel.m_dot": 0.03}

    if xi is not None:
        loss = PressureLossElement("loss", "comb", "post", correlation=ConstantFractionLoss(xi=xi))
        loss.initial_guess = {"loss.m_dot": 1.03}
        net.add_element(loss)
    else:
        mid = LosslessConnectionElement("loss", "comb", "post")
        mid.initial_guess = {"loss.m_dot": 1.03}
        net.add_element(mid)

    nozzle = OrificeElement("nozzle", "post", "out", Cd=0.8, diameter=0.12, correlation="fixed")
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


def _make_effective_area_inlet_network(loss_correlation=None):
    """Build network: MassFlowBoundary x2 -> EffectiveArea -> Combustor -> [loss?] -> Orifice -> PressureBoundary.

    When ``loss_correlation`` is provided, a :class:`PressureLossElement` with
    that correlation sits between the combustor and a new ``post`` plenum; the
    orifice is downstream of ``post``.  Otherwise the combustor connects to the
    orifice via a lossless link.
    """
    Y_air = cb.species.dry_air_mass()
    Y_ch4 = cb.species.pure_species("CH4")

    net = FlowNetwork()
    air = MassFlowBoundary("air", m_dot=1.0, T_total=600.0, Y=Y_air)
    air.initial_guess = {"air.P_total": 160000.0, "air.P": 160000.0}
    fuel = MassFlowBoundary("fuel", m_dot=0.03, T_total=300.0, Y=Y_ch4)
    fuel.initial_guess = {"fuel.P_total": 160000.0, "fuel.P": 160000.0}
    comb = CombustorNode("comb", method="complete")
    comb.initial_guess = {"comb.P_total": 150000.0, "comb.P": 140000.0}
    post = PlenumNode("post")
    post.initial_guess = {"post.P_total": 135000.0, "post.P": 135000.0}
    out = PressureBoundary("out", P_total=101325.0, T_total=300.0)
    for n in (air, fuel, comb, post, out):
        net.add_node(n)

    e_air = EffectiveAreaConnectionElement("e_air", "air", "comb", effective_area=0.05)
    e_air.initial_guess = {"e_air.m_dot": 1.0}
    e_fuel = EffectiveAreaConnectionElement("e_fuel", "fuel", "comb", effective_area=0.005)
    e_fuel.initial_guess = {"e_fuel.m_dot": 0.03}

    if loss_correlation is not None:
        loss = PressureLossElement("loss", "comb", "post", correlation=loss_correlation)
        loss.initial_guess = {"loss.m_dot": 1.03}
        net.add_element(loss)
    else:
        mid = LosslessConnectionElement("loss", "comb", "post")
        mid.initial_guess = {"loss.m_dot": 1.03}
        net.add_element(mid)

    nozzle = OrificeElement("nozzle", "post", "out", Cd=0.8, diameter=0.12, correlation="fixed")
    nozzle.initial_guess = {"nozzle.m_dot": 1.03}
    for e in (e_air, e_fuel, nozzle):
        net.add_element(e)

    return net, comb


def test_combustor_pressure_loss_applied_at_post_plenum():
    """Sanity check: PressureLossElement enforces ``post/comb = (1 - xi)`` exactly."""
    xi = 0.05
    net, _ = _make_effective_area_inlet_network(loss_correlation=ConstantFractionLoss(xi=xi))
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"]
    ratio = sol["post.P_total"] / sol["comb.P_total"]
    assert ratio == pytest.approx(1.0 - xi, rel=1e-4)


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


@pytest.mark.parametrize("loss_model,label", LOSS_MODELS)
def test_combustor_pressure_loss_affects_solution(loss_model, label):
    """The solution with loss MUST differ from the lossless baseline.

    Physically: with fixed outlet P_total and fixed inlet mass flows, introducing
    a combustor pressure loss must raise the upstream MassFlowBoundary P_total so
    the nozzle still passes the required flow against the same back pressure.
    """
    # Lossless baseline
    net0, _ = _make_effective_area_inlet_network(loss_correlation=None)
    sol0 = NetworkSolver(net0).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol0["__success__"], f"Lossless baseline failed: |F|={sol0['__final_norm__']:.2e}"

    # Lossy case
    net, _ = _make_effective_area_inlet_network(loss_correlation=loss_model)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"], f"{label} solver failed: |F|={sol['__final_norm__']:.2e}"

    # Upstream pressure MUST rise to overcome the loss. Require at least 0.1% difference.
    rel_diff = abs(sol["air.P_total"] - sol0["air.P_total"]) / sol0["air.P_total"]
    assert rel_diff > 1e-3, (
        f"{label}: loss not enforced - air.P_total={sol['air.P_total']:.2f} vs "
        f"lossless={sol0['air.P_total']:.2f} (rel diff {rel_diff:.2e})"
    )


# ---------------------------------------------------------------------------
# Exact-physics assertions: once the plan lands, these must pass with tight
# tolerances, not just "loss has some effect". They guard against a sloppy
# implementation that only partially enforces the loss.
# ---------------------------------------------------------------------------


def test_combustor_constant_fraction_loss_exact_pressure_rise():
    """ConstantFractionLoss(xi=0.05): upstream P_total must rise by 1/(1-xi) - 1 = 5.263%.

    With fixed outlet P_total and fixed inlet m_dot, the nozzle (orifice) dictates
    a required upstream-of-orifice stagnation pressure. Introducing a 5% combustor
    total-pressure loss forces the MassFlowBoundary pressures up by exactly the
    inverse factor 1/(1-xi) to deliver the same flow against the same back pressure.
    """
    xi = 0.05

    # Lossless baseline
    net0, _ = _make_effective_area_inlet_network(loss_correlation=None)
    sol0 = NetworkSolver(net0).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol0["__success__"]

    # Lossy case
    net, _ = _make_effective_area_inlet_network(loss_correlation=ConstantFractionLoss(xi=xi))
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"], f"|F|={sol['__final_norm__']:.2e}"

    # Upstream pressure rise: approximately 1/(1-xi) - 1 = 5.26% for xi=0.05.
    # The exact factor depends on the pressure-flow law of the upstream element;
    # allow a ~1% tolerance on the rise magnitude.
    expected_rise = 1.0 / (1.0 - xi) - 1.0
    actual_rise = (sol["air.P_total"] - sol0["air.P_total"]) / sol0["air.P_total"]
    assert actual_rise == pytest.approx(expected_rise, rel=1e-2), (
        f"Upstream P_total rise {actual_rise * 100:.4f}% != "
        f"expected {expected_rise * 100:.4f}% (1/(1-xi)-1)"
    )

    # Post-loss plenum P_total must be (1-xi) * combustor P_total exactly.
    ratio = sol["post.P_total"] / sol["comb.P_total"]
    assert ratio == pytest.approx(1.0 - xi, rel=1e-4), (
        f"post.P_total / comb.P_total = {ratio:.5f} != {1.0 - xi:.5f}"
    )


def test_combustor_constant_head_loss_exact_pressure_rise():
    """ConstantHeadLoss(zeta): the loss fraction xi = zeta * 0.5*rho*v^2 / P_in.

    For the same network with fixed inlet m_dot and outlet P, the solver must
    satisfy post.P_total = (1 - xi_converged) * comb.P_total where
    xi_converged is the dynamic-head fraction at the converged flow state.
    """
    zeta = 5.0
    area = 0.05

    # Lossy case with zeta-based head loss
    net, _ = _make_effective_area_inlet_network(
        loss_correlation=ConstantHeadLoss(zeta=zeta, area=area)
    )
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"], f"|F|={sol['__final_norm__']:.2e}"

    # xi at converged state = 1 - post/comb (total-pressure ratio)
    xi_eff = 1.0 - sol["post.P_total"] / sol["comb.P_total"]
    assert xi_eff > 0.001, f"Converged xi too small: {xi_eff:.2e}"

    # Upstream P_total rise relative to lossless baseline matches 1/(1-xi_eff) - 1
    net0, _ = _make_effective_area_inlet_network(loss_correlation=None)
    sol0 = NetworkSolver(net0).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol0["__success__"]
    expected_rise = 1.0 / (1.0 - xi_eff) - 1.0
    actual_rise = (sol["air.P_total"] - sol0["air.P_total"]) / sol0["air.P_total"]
    assert actual_rise == pytest.approx(expected_rise, rel=5e-3), (
        f"Zeta loss: upstream P_total rise {actual_rise * 100:.4f}% != "
        f"expected {expected_rise * 100:.4f}% (based on converged xi={xi_eff:.5f})"
    )


def test_combustor_lossless_inlets_pressure_loss_enforced():
    """Regression: pressure loss downstream of a lossless-inlet combustor is enforced.

    Topology: MassFlowBoundary x2 -> Lossless links -> Combustor -> PressureLossElement
    -> post plenum -> Orifice -> PressureBoundary. The loss ``post/comb = 1 - xi``
    must hold exactly after convergence.
    """
    xi = 0.05
    net, _ = _make_lossless_inlet_network(xi=xi)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"]

    ratio = sol["post.P_total"] / sol["comb.P_total"]
    assert ratio == pytest.approx(1.0 - xi, rel=1e-4), (
        f"Loss not enforced: post/comb = {ratio:.5f}, expected {1.0 - xi:.5f}"
    )


# ---------------------------------------------------------------------------
# PressureLossElement-specific validation and override behaviours.
# ---------------------------------------------------------------------------


def test_pressure_loss_cold_flow_warning_without_theta_source():
    """Linear-theta correlation on a non-combustor edge emits UserWarning and
    falls back to the cold-flow coefficient (theta = 0 -> xi = xi0)."""
    import warnings

    net = FlowNetwork()
    # Two generic plenums, no combustor anywhere.
    p_in = PressureBoundary("p_in", P_total=200000.0, T_total=500.0)
    p_out = PressureBoundary("p_out", P_total=101325.0, T_total=500.0)
    p_mid = PlenumNode("mid")
    p_mid.initial_guess = {"mid.P_total": 150000.0, "mid.P": 150000.0}
    for n in (p_in, p_mid, p_out):
        net.add_node(n)

    # Placement: loss element between mid (plenum) and p_out.
    # Neither endpoint has_theta -> cold-flow fallback + warning.
    xi0 = 0.10
    loss = PressureLossElement(
        "loss", "mid", "p_out", correlation=LinearThetaFractionLoss(k=0.8, xi0=xi0)
    )
    loss.initial_guess = {"loss.m_dot": 0.5}
    # Provide an inflow element to reach the plenum.
    feed = OrificeElement("feed", "p_in", "mid", Cd=0.8, diameter=0.1, correlation="fixed")
    feed.initial_guess = {"feed.m_dot": 0.5}
    net.add_element(feed)
    net.add_element(loss)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})

    msgs = [str(wi.message) for wi in w if issubclass(wi.category, UserWarning)]
    assert any("linear-theta" in m and "loss" in m for m in msgs), (
        f"Expected cold-flow warning; got: {msgs!r}"
    )
    assert sol["__success__"]
    # Cold-flow: theta = 0 -> xi = xi0. Ratio = P_out_boundary / mid.P_total.
    P_out_fixed = 101325.0
    ratio = P_out_fixed / sol["mid.P_total"]
    assert ratio == pytest.approx(1.0 - xi0, rel=1e-4)


def test_pressure_loss_ambiguous_both_endpoints_theta_raises():
    """PressureLossElement between two has_theta=True nodes raises at validate."""
    net = FlowNetwork()
    p_in = PressureBoundary("p_in", P_total=200000.0, T_total=500.0)
    p_out = PressureBoundary("p_out", P_total=101325.0, T_total=500.0)
    c1 = CombustorNode("c1", method="complete")
    c1.initial_guess = {"c1.P_total": 180000.0, "c1.P": 180000.0}
    c2 = CombustorNode("c2", method="complete")
    c2.initial_guess = {"c2.P_total": 160000.0, "c2.P": 160000.0}
    for n in (p_in, c1, c2, p_out):
        net.add_node(n)

    feed = LosslessConnectionElement("feed", "p_in", "c1")
    loss = PressureLossElement("loss", "c1", "c2", correlation=ConstantFractionLoss(xi=0.05))
    tail = OrificeElement("tail", "c2", "p_out", Cd=0.8, diameter=0.1, correlation="fixed")
    for e in (feed, loss, tail):
        net.add_element(e)

    with pytest.raises(ValueError, match="both endpoints"):
        NetworkSolver(net).solve(method="hybr", use_jac=True)


def test_pressure_loss_explicit_theta_source_override():
    """theta_source='<id>' disambiguates when both endpoints are combustors."""
    net = FlowNetwork()
    p_in = PressureBoundary("p_in", P_total=400000.0, T_total=500.0)
    p_out = PressureBoundary("p_out", P_total=101325.0, T_total=500.0)
    c1 = CombustorNode("c1", method="complete")
    c1.initial_guess = {"c1.P_total": 380000.0, "c1.P": 380000.0}
    c2 = CombustorNode("c2", method="complete")
    c2.initial_guess = {"c2.P_total": 340000.0, "c2.P": 340000.0}
    for n in (p_in, c1, c2, p_out):
        net.add_node(n)

    feed = LosslessConnectionElement("feed", "p_in", "c1")
    # Explicit: use c1's theta (not c2's).  Correlation is constant-fraction so
    # theta value is irrelevant for xi but the resolution must not raise.
    loss = PressureLossElement(
        "loss",
        "c1",
        "c2",
        correlation=ConstantFractionLoss(xi=0.05),
        theta_source="c1",
    )
    tail = OrificeElement("tail", "c2", "p_out", Cd=0.8, diameter=0.1, correlation="fixed")
    for e in (feed, loss, tail):
        net.add_element(e)

    sol = NetworkSolver(net).solve(method="hybr", use_jac=True, options={"xtol": 1e-12})
    assert sol["__success__"]
    # Loss must still be applied exactly.
    ratio = sol["c2.P_total"] / sol["c1.P_total"]
    assert ratio == pytest.approx(0.95, rel=1e-4)


def test_pressure_loss_invalid_theta_source_raises():
    """theta_source pointing to a has_theta=False node raises at resolve."""
    net = FlowNetwork()
    p_in = PressureBoundary("p_in", P_total=200000.0, T_total=500.0)
    p_out = PressureBoundary("p_out", P_total=101325.0, T_total=500.0)
    mid = PlenumNode("mid")
    for n in (p_in, mid, p_out):
        net.add_node(n)

    feed = OrificeElement("feed", "p_in", "mid", Cd=0.8, diameter=0.1, correlation="fixed")
    # 'mid' is a plenum, has_theta=False -> must raise.
    loss = PressureLossElement(
        "loss",
        "mid",
        "p_out",
        correlation=ConstantFractionLoss(xi=0.05),
        theta_source="mid",
    )
    for e in (feed, loss):
        net.add_element(e)

    with pytest.raises(ValueError, match="does not provide theta"):
        NetworkSolver(net).solve(method="hybr", use_jac=True)


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


def test_combustor_diagnostics_exposes_m_dot():
    """CombustorNode.diagnostics must expose the total mass flow as ``m_dot``.

    The GUI's Live Telemetry panel renders `state.m_dot` in a prominent slot;
    without this field combustor nodes would silently display ṁ = 0.
    """
    net, _ = _make_lossless_inlet_network(xi=None)
    sol = NetworkSolver(net).solve(method="hybr", use_jac=True)
    assert sol["__success__"], f"|F|={sol['__final_norm__']:.2e}"

    # Flat diagnostic key exists and equals the expected mass flow.
    assert "comb.m_dot" in sol, "CombustorNode.diagnostics must expose 'm_dot'"
    m_comb = sol["comb.m_dot"]
    m_in = sol["l_air.m_dot"] + sol["l_fuel.m_dot"]
    assert abs(m_comb - m_in) < 1e-6, f"Combustor m_dot {m_comb:.6f} != inlet sum {m_in:.6f}"
    assert abs(m_comb - 1.03) < 1e-6
