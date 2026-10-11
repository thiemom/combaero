"""IdealMomentumTeeElement: the ideal momentum junction kernel in a network (#517).

The kernel's closed forms are held in test_ideal_junction.py/.cpp. What is new
here is the wiring: the solver's unknowns reach the kernel in the right slots,
the angle convention maps the tee's port angles onto outward axes, ``P_jct`` is
a static pressure, and the solve converges -- including through a reversed
port, which the closure models and the Mynard tee refuses.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.optimize._numdiff import approx_derivative

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    IdealMomentumTeeElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from validation.junction import joining_variants as jv

_D = 0.05
_A = math.pi * (_D / 2.0) ** 2
_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _branch_net(
    m_in: float = 0.05,
    theta: float = 90.0,
    regime: str = "incompressible",
    pt_bra: float = 1.9e5,
) -> FlowNetwork:
    """MFB-fed dividing tee into two PB sinks through long arm channels."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=m_in, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_str", Pt=1.9e5, Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_bra", Pt=pt_bra, Tt=300.0, Y=_Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(
        ChannelElement("ch_in", "src", "mc_com", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_str", "mc_str", "pb_str", length=5.0, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_bra", "mc_bra", "pb_bra", length=5.0, diameter=_D, regime=regime)
    )
    net.add_element(
        IdealMomentumTeeElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, theta],
        )
    )
    return net


def _merge_net(
    m: tuple[float, float] = (0.03, 0.02), theta: float = 90.0, regime: str = "incompressible"
) -> FlowNetwork:
    """Two MFB suppliers joining into one PB sink."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_str", m_dot=m[0], Tt=300.0, Y=_Y))
    net.add_node(MassFlowBoundary("mfb_bra", m_dot=m[1], Tt=300.0, Y=_Y))
    net.add_node(PressureBoundary("pb_out", Pt=1.9e5, Tt=300.0, Y=_Y))
    for nid in ("mc_com", "mc_str", "mc_bra"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(
        ChannelElement("ch_str", "mfb_str", "mc_str", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_bra", "mfb_bra", "mc_bra", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_out", "mc_com", "pb_out", length=0.5, diameter=_D, regime=regime)
    )
    net.add_element(
        IdealMomentumTeeElement(
            id="jct",
            inlet_nodes=["mc_str", "mc_bra"],
            outlet_nodes=["mc_com"],
            inlet_angles_deg=[0.0, theta],
            outlet_angles_deg=[0.0],
        )
    )
    return net


def _solve(net: FlowNetwork) -> tuple[dict, NetworkSolver]:
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    return sol, solver


def _rho(sol: dict, node: str) -> float:
    return float(sol[f"{node}.rho"])


def _receiver_identity(sol: dict, node: str, m_out: float, V_n: float) -> float:
    """P_j - p_jn + rho_j u_j (u_j - V.n_j), in Pa: the receiver row."""
    u = m_out / (_rho(sol, node) * _A)
    return sol[f"{node}.P"] - sol["jct.P_jct"] + _rho(sol, node) * u * (u - V_n)


# --- the rows hold on the converged network ---------------------------------


@pytest.mark.parametrize("theta", [45.0, 90.0, 135.0])
def test_dividing_rows_hold_at_the_converged_state(theta: float) -> None:
    """Supplier at p_jn; each receiver's momentum balance with V = v_supplier.

    Exact up to the role-smoothing bias, which is tiny next to the head.
    """
    sol, _ = _solve(_branch_net(theta=theta))
    m_c, m_s, m_b = sol["ch_in.m_dot"], sol["ch_str.m_dot"], sol["ch_bra.m_dot"]
    assert m_c > 0.0 and m_s > 0.0 and m_b > 0.0
    head = sol["mc_com.Pt"] - sol["mc_com.P"]
    assert sol["mc_com.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    # Supplier velocity along the main axis: the straight receiver sees all of
    # it, the branch receiver its projection on the branch axis.
    u_c = m_c / (_rho(sol, "mc_com") * _A)
    assert _receiver_identity(sol, "mc_str", m_s, u_c) == pytest.approx(0.0, abs=1e-4 * head)
    V_n = u_c * math.cos(math.radians(theta))
    assert _receiver_identity(sol, "mc_bra", m_b, V_n) == pytest.approx(0.0, abs=1e-4 * head)


def test_dividing_K_is_vector_borda_at_low_mach() -> None:
    """K_straight = (1 - q_s)^2, K_branch = 1 + q_b^2 - 2 q_b cos(theta)
    (equal areas), against the port nodes' own Pt. M ~ 0.03, so the
    isentropic port heads differ from rho u^2/2 by ~2e-4: the tolerance."""
    theta = 60.0
    sol, _ = _solve(_branch_net(theta=theta))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    q_s = sol["ch_str.m_dot"] / sol["ch_in.m_dot"]
    q_b = sol["ch_bra.m_dot"] / sol["ch_in.m_dot"]
    K_s = (sol["mc_com.Pt"] - sol["mc_str.Pt"]) / q_c
    K_b = (sol["mc_com.Pt"] - sol["mc_bra.Pt"]) / q_c
    assert K_s == pytest.approx((1.0 - q_s) ** 2, abs=2e-3)
    assert K_b == pytest.approx(1.0 + q_b**2 - 2.0 * q_b * math.cos(math.radians(theta)), abs=2e-3)


@pytest.mark.parametrize("theta", [45.0, 90.0, 120.0])
def test_merge_K_is_blaisdell_manson_theory_at_low_mach(theta: float) -> None:
    """The converged network's K11/K12 are the scored ``BM theory`` variant."""
    sol, _ = _solve(_merge_net(theta=theta))
    q_c = sol["mc_com.Pt"] - sol["mc_com.P"]
    q_lat = sol["ch_bra.m_dot"] / sol["ch_out.m_dot"]
    U, A, th = jv.tee(q_lat, 1.0, theta)
    ref = jv.MODELS["BM theory"](U, A, th)
    assert (sol["mc_str.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[1], abs=2e-3)
    assert (sol["mc_bra.Pt"] - sol["mc_com.Pt"]) / q_c == pytest.approx(ref[2], abs=2e-3)


# --- solve behaviour ----------------------------------------------------------


def _reversed_branch_solution() -> dict:
    """Declared dividing, but the branch sink sits ABOVE the junction, so the
    branch flows back in and the tee runs as a merge."""
    sol = NetworkSolver(_branch_net(pt_bra=1.95e5)).solve(timeout=60.0)
    assert sol["ch_bra.m_dot"] < 0.0  # back into the junction
    return sol


def test_a_reversed_converged_port_is_demoted() -> None:
    """The closure's rows hold through the reversal -- both suppliers sit at
    p_jn -- but the network cannot represent the reversed port, so the
    answer is reported as inconsistent rather than returned as converged."""
    sol = _reversed_branch_solution()
    head = sol["mc_str.Pt"] - sol["mc_str.P"]
    assert sol["mc_com.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    assert sol["mc_bra.P"] == pytest.approx(sol["jct.P_jct"], abs=1e-4 * head)
    assert sol["__consistent__"] is False
    assert not sol["__success__"]


@pytest.mark.xfail(
    strict=True,
    reason="junction-port state propagation follows the declared direction: a "
    "reversed port node sees zero flow and drops its head (Pt = P)",
)
def test_a_reversed_port_node_keeps_its_dynamic_head() -> None:
    """Flips when the solver feeds junction ports by the sign of their flow.
    Then the reversal demotion in verify_solution_consistent can go."""
    sol = _reversed_branch_solution()
    m = sol["ch_bra.m_dot"]
    rho = sol["mc_bra.rho"]
    head = 0.5 * rho * (m / (rho * _A)) ** 2
    assert sol["mc_bra.Pt"] - sol["mc_bra.P"] == pytest.approx(head, rel=1e-2)


@pytest.mark.parametrize("builder", [_branch_net, _merge_net])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_network_jacobian_matches_central_differences(builder, regime: str) -> None:
    """At a converged state perturbed by 1e-3, NOT at x0.

    x0 seeds every pressure equal, which sits a compressible channel exactly
    on its direction switch, where its own derivative is one-sided; that is
    the channel's, not this element's, and a check there measures the kink.
    Compared as J_ij |x_j| against the row's largest such term.
    """
    net = builder(regime=regime)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, use_jac=False)
    assert sol["__success__"], sol.get("__message__")
    x = np.array([sol[n] for n in solver.unknown_names], dtype=float)
    x = x * (1.0 + 1e-3 * np.random.default_rng(0).standard_normal(x.size))
    _, J = solver._residuals_and_jacobian(x)
    J = J.toarray()
    step = np.maximum(np.abs(x) * 1e-7, 1e-10)
    J_fd = approx_derivative(solver._residuals, x, method="3-point", abs_step=step)
    scaled_err = np.abs(J - J_fd) * np.abs(x)
    row_ref = np.max(np.abs(J_fd) * np.abs(x), axis=1) + 1e-30
    assert float(np.max(scaled_err.max(axis=1) / row_ref)) < 1e-5


def test_compressible_merge_converges_with_the_analytic_jacobian_from_a_warm_start() -> None:
    """From the FD solve's root perturbed by 1%, the analytic Jacobian walks
    back: the element's (f, J) is good away from the x0 tie."""
    net = _merge_net(m=(0.18, 0.12), theta=45.0, regime="compressible")
    solver = NetworkSolver(net)
    ref = solver.solve(timeout=60.0, use_jac=False)
    assert ref["__success__"]
    x0 = np.array([ref[n] for n in solver.unknown_names], dtype=float) * 1.01
    sol = NetworkSolver(_merge_net(m=(0.18, 0.12), theta=45.0, regime="compressible")).solve(
        timeout=60.0, x0=x0
    )
    assert sol["__success__"], sol.get("__message__")
    assert sol["jct.P_jct"] == pytest.approx(ref["jct.P_jct"], rel=1e-8)


def test_static_seed_moves_P_jct_below_the_common_Pt() -> None:
    """``analytical_pt_prop`` seeds P_jct from the common port's Pt; this
    closure's P_jct is static, so its seed sits one common head lower."""
    net = _branch_net(m_in=0.3)
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    ref = solver._infer_reference_state()
    overrides = solver._propagate_analytical_pt_prop(ref)
    pt = solver._propagate_pressure_guess(ref)["mc_com"]
    m = overrides["ch_in.m_dot"]
    rho = max(pt / (287.0 * float(ref.get("T", 300.0))), 1e-3)
    assert overrides["jct.P_jct"] == pytest.approx(pt - m * m / (2.0 * rho * _A * _A), rel=1e-12)
    assert overrides["jct.P_jct"] < pt


# --- element-level contracts --------------------------------------------------


def test_axes_follow_the_tee_angle_convention() -> None:
    """Outlet: the angle. Inlet: pi - angle (its flow arrives along the angle)."""
    el = IdealMomentumTeeElement(
        "j",
        inlet_nodes=["a", "b"],
        outlet_nodes=["c"],
        inlet_angles_deg=[0.0, 45.0],
        outlet_angles_deg=[30.0],
        port_areas=[1.0, 1.0, 1.0],
    )
    assert el.axis_rad() == pytest.approx(
        [math.pi, math.pi - math.radians(45.0), math.radians(30.0)]
    )


def test_more_ports_than_the_kernel_takes_is_refused() -> None:
    with pytest.raises(ValueError, match="at most 6"):
        IdealMomentumTeeElement("j", inlet_nodes=["a", "b", "c"], outlet_nodes=["d", "e", "f", "g"])


def test_four_port_cross_converges_with_both_side_arms_receiving() -> None:
    """N > 3 is in scope (the Mynard tee refuses it): a cross with one
    supply, a straight run and two opposed 90-deg side arms splits the side
    flow evenly by symmetry."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("src", m_dot=0.05, Tt=300.0, Y=_Y))
    for nid in ("pb_s", "pb_l", "pb_r"):
        net.add_node(PressureBoundary(nid, Pt=1.9e5, Tt=300.0, Y=_Y))
    for nid in ("mc_c", "mc_s", "mc_l", "mc_r"):
        net.add_node(MomentumChamberNode(nid, area=_A))
    net.add_element(ChannelElement("ch_in", "src", "mc_c", length=0.5, diameter=_D))
    for arm in ("s", "l", "r"):
        net.add_element(
            ChannelElement(f"ch_{arm}", f"mc_{arm}", f"pb_{arm}", length=5.0, diameter=_D)
        )
    net.add_element(
        IdealMomentumTeeElement(
            "jct",
            inlet_nodes=["mc_c"],
            outlet_nodes=["mc_s", "mc_l", "mc_r"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0, -90.0],
        )
    )
    sol, _ = _solve(net)
    assert sol["ch_l.m_dot"] == pytest.approx(sol["ch_r.m_dot"], rel=1e-8)
    assert sol["ch_s.m_dot"] > sol["ch_l.m_dot"] > 0.0


def test_energy_verifier_rejects_a_junction_that_creates_flow_work() -> None:
    net = _merge_net()
    net.resolve_all_topology()
    jct = net.elements["jct"]
    flows = {"ch_str.m_dot": 0.03, "ch_bra.m_dot": 0.02, "ch_out.m_dot": 0.05}
    good = {**flows, "mc_str.Pt": 1.9010e5, "mc_bra.Pt": 1.9012e5, "mc_com.Pt": 1.9005e5}
    bad = {**flows, "mc_str.Pt": 1.9010e5, "mc_bra.Pt": 1.9012e5, "mc_com.Pt": 1.9020e5}
    assert jct.verify_solution_consistent(good) is True
    assert jct.verify_solution_consistent(bad) is False
    # A dead port is fine, a reversed one is not (see the xfail above).
    dead = {**good, "ch_bra.m_dot": -1e-9, "ch_str.m_dot": 0.05}
    reversed_ = {**good, "ch_bra.m_dot": -0.01, "ch_str.m_dot": 0.06}
    assert jct.verify_solution_consistent(dead) is True
    assert jct.verify_solution_consistent(reversed_) is False


def test_diagnostics_name_the_regime_from_the_converged_flows() -> None:
    sol, solver = _solve(_branch_net(theta=90.0))
    net = solver.network
    jct = net.elements["jct"]
    states = [solver._get_node_state(net.nodes[n], solver._last_x) for n in jct.port_nodes]
    mdots = [jct._port_signs[i] * sol[f"{jct._port_element_ids[i]}.m_dot"] for i in range(3)]
    diag = jct.diagnostics(states, sol["jct.P_jct"], mdots)
    assert diag["p_jn"] == sol["jct.P_jct"]
    assert diag["mass_flow_ratio"] == pytest.approx(sol["ch_bra.m_dot"] / sol["ch_in.m_dot"])
    q_b = diag["mass_flow_ratio"]
    assert diag["K_branch"] == pytest.approx(1.0 + q_b * q_b, abs=2e-3)
    assert "K11" not in diag
