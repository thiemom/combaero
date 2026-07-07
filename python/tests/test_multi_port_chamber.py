"""
Network-level tests for MultiPortChamberElement (momentum-CV junction).

The element + companion BorderCarnotLossElement land the PDF-spec successor
to the K-closure tee (docs/junction/momentum cv implementation guide.pdf).
These tests cover the assembly side: the junction wires correctly into
combaero's solver, the port-MCN mass-row skip works, sign conventions
resolve cleanly from topology, and the system converges to a mass-conserving
root on a simple low-Mach 3-port network.

Tier-4 cross-check vs TeeJunctionElement and Tier-1/2 validation against
Hager / Bassett / Perez-Garcia / Wang reference data live in separate files.
"""

import math

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _three_port_net(
    Pt_in: float = 2.1e5,
    Pt_str: float = 2.05e5,
    Pt_bra: float = 2.00e5,
    Tt: float = 300.0,
    diameter: float = 0.05,
    length: float = 0.3,
) -> FlowNetwork:
    """3-port momentum-CV junction: 1 pressure inlet, 2 pressure outlets, all
    via finite-length compressible channels with realistic friction.

    Topology:
       pb_in  --ch_in-->  mc_com \\
                                  >- jct -< mc_str --ch_str--> pb_str
                                            mc_bra --ch_bra--> pb_bra

    All three ports are MomentumChamberNodes (PDF requires per-port (P, Pt)
    state at the junction face). The junction is lossless; channel friction
    sets up the Pt gradient that drives flow.
    """
    Y = _DRY_AIR_Y
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=Pt_in, Tt=Tt, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=Tt, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_bra, Tt=Tt, Y=Y))
    net.add_node(MomentumChamberNode("mc_com", area=math.pi * (diameter / 2.0) ** 2))
    net.add_node(MomentumChamberNode("mc_str", area=math.pi * (diameter / 2.0) ** 2))
    net.add_node(MomentumChamberNode("mc_bra", area=math.pi * (diameter / 2.0) ** 2))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_in", "mc_com", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_str", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra", "mc_bra", "pb_bra", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
        )
    )
    return net


@pytest.fixture(scope="module")
def three_port_solution() -> tuple[NetworkSolver, dict]:
    """One shared solve of the default 3-port net for the test_three_port_*
    assertions below.

    The solve is seeded via init_strategy='outletref_warmstart' to pin
    the PHYSICAL root: this net sits in the ejector regime (the 2.05 bar
    straight sink is entrained by the 2.1 bar feed and everything exits
    through the 2.0 bar branch -- Bassett 2001 Section 3), and the
    sign-symmetric v1 impulse rows admit an exact MIRROR root as well
    (flow entering from the LOWEST-Pt sink and exiting to the highest --
    energetically impossible for a passive junction; the known v1
    joining-flow inconsistency). The plain cold start converges onto
    that mirror root; the v1 energy-consistency hook now demotes it and
    the auto-retry self-heals onto the physical root (see
    test_three_port_default_path_self_heals_from_mirror_root), but the
    explicit seed keeps this shared fixture on the fast direct path.
    Since the 2026-07-07 residual-scale fix all init paths converge in
    seconds (this net was previously the ~187 s "doomed primary"
    exemplar).
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=120.0, init_strategy="outletref_warmstart")
    return solver, sol


def test_three_port_network_converges(three_port_solution):
    """Solver assembles and converges on the basic 3-port momentum-CV network."""
    solver, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Residual norm at the converged state is small.
    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert max(abs(r) for r in res) < 1e-3, f"residuals not small: {res}"


def test_residual_scales_classify_mpce_impulse_rows_as_pressure():
    """Regression for the 2026-07-07 scaling fix: on the 3-port net the
    ONLY ref_mdot-scaled row is the junction sum-mass row (port-MCN mass
    rows are skipped; MCN closures, channel rows, and the N impulse rows
    are all pressure-magnitude). Before the fix all N+1 MPCE rows were
    ref_mdot, over-weighting the impulse rows by ref_p/ref_mdot (~1e5)
    in the scaled system -- the root cause of the doomed-primary stall
    class on MPCE networks.
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = solver._build_x0()
    ref_p, ref_mdot = solver._reference_scales(x0)
    scales = solver._build_residual_scales(x0)
    assert len(scales) == len(x0)
    n_mdot_rows = int(np.sum(scales == ref_mdot))
    n_p_rows = int(np.sum(scales == ref_p))
    assert n_mdot_rows == 1, f"expected 1 ref_mdot row (junction sum-mass), got {n_mdot_rows}"
    assert n_p_rows == len(scales) - 1


def _resolved_junction():
    net = _three_port_net()
    net.resolve_all_topology()
    return net.elements["jct"]


def _sol_dict(m_in: float, m_str: float, m_bra: float, pt_com: float, pt_str: float, pt_bra: float):
    return {
        "ch_in.m_dot": m_in,
        "ch_str.m_dot": m_str,
        "ch_bra.m_dot": m_bra,
        "mc_com.Pt": pt_com,
        "mc_str.Pt": pt_str,
        "mc_bra.Pt": pt_bra,
    }


def test_v1_energy_check_rejects_mirror_root():
    """Suppliers at ~2.0 bar feeding a ~2.11 bar collector: the exact
    sign-flipped image of the physical root, energetically impossible
    for a passive junction."""
    jct = _resolved_junction()
    mirror = _sol_dict(-0.43, 0.29, -0.72, 2.11e5, 2.05e5, 1.99e5)
    assert jct.verify_solution_consistent(mirror) is False


def test_v1_energy_check_accepts_physical_ejector_root():
    """The 2.1 bar primary entrains the 2.05 bar straight arm and all
    flow exits at the 2.0 bar branch: reversal at a port is physical
    (Bassett 2001 Section 3). Two suppliers => joining mode => the
    check is out of scope (and the state is energetically fine anyway)."""
    jct = _resolved_junction()
    ejector = _sol_dict(0.36, -0.21, 0.57, 2.09e5, 2.04e5, 2.02e5)
    assert jct.verify_solution_consistent(ejector) is True


def test_v1_energy_check_skips_joining_mode():
    """>= 2 suppliers (joining) is out of scope even when a collector
    exceeds every supplier Pt: the v1 impulse rows themselves can
    manufacture flow work for joining flow (the documented deficiency
    that motivated MPCEv2), so policing energy here would reject the
    v1 merge model's own legitimate solutions."""
    jct = _resolved_junction()
    # Suppliers: the common port (element-frame +0.4, junction-convention
    # -0.4) and the reversed straight arm (-0.25); collector: the branch
    # at 2.2e5 Pa, above both supplier Pts (2.05e5 / 2.0e5).
    joining = _sol_dict(0.4, -0.25, 0.65, 2.05e5, 2.0e5, 2.2e5)
    assert jct.verify_solution_consistent(joining) is True


def test_v1_energy_check_zero_flow_is_vacuous():
    jct = _resolved_junction()
    idle = _sol_dict(1e-12, -1e-12, 1e-12, 1.5e5, 1.5e5, 1.5e5)
    assert jct.verify_solution_consistent(idle) is True


def test_v1_energy_check_missing_keys_pass():
    jct = _resolved_junction()
    partial = {"ch_in.m_dot": -0.43, "mc_com.Pt": 2.11e5}
    assert jct.verify_solution_consistent(partial) is True


def test_v1_energy_check_tolerates_near_equal_pt():
    """A collector within rel_tol of the max supplier Pt is numerical
    noise, not manufactured work."""
    jct = _resolved_junction()
    near = _sol_dict(0.4, 0.2, 0.2, 2.0e5, 2.0e5 * (1.0 + 0.5e-4), 1.9e5)
    assert jct.verify_solution_consistent(near) is True


def test_three_port_default_path_self_heals_from_mirror_root():
    """End-to-end payoff of the v1 energy hook: the plain cold solve on
    this net converges onto the mirror root, the post-solve verification
    demotes it, and the outlet-ref auto-retry then reaches the physical
    ejector root -- success with forward common inflow, no manual
    strategy selection needed.
    """
    net = _three_port_net()
    solver = NetworkSolver(net)
    with pytest.warns(UserWarning):
        sol = solver.solve(timeout=120.0)
    assert sol["__success__"], sol.get("__message__")
    assert sol["ch_in.m_dot"] > 0.0, f"inlet reversed: {sol['ch_in.m_dot']}"
    ssu = solver._diagnostic_data["solver_settings_used"]
    assert ssu["init_strategy"] == "outletref_warmstart"
    assert ssu.get("auto_retry") is True


def test_stall_wiring_forced_detector(monkeypatch):
    """Deterministic wiring test for the stall machinery (its natural
    fixture, this net's doomed primary, was cured by the residual-scale
    fix): force the detector to fire in both phases and assert the
    hybr handover and the LM-phase abort route correctly.
    """
    from combaero.network import solver as solver_module

    monkeypatch.setattr(solver_module, "_stall_detected", lambda *a, **k: True)
    net = _three_port_net()
    solver = NetworkSolver(net)
    with pytest.warns(UserWarning, match="did not converge"):
        sol = solver.solve(timeout=30.0, auto_retry=False)
    assert not sol["__success__"]
    diag = solver._diagnostic_data
    assert diag["stall_handover"] is True
    assert diag["lm_stall"] is True
    assert "LM fallback also stalled" in str(sol.get("__message__"))


def test_three_port_mass_conservation(three_port_solution):
    """All three port flows balance at the junction (sum-mass residual = 0).

    Junction convention: ch_in is the declared inlet (sign -1: +ch_in is flow
    INTO junction), ch_str and ch_bra are outlets (sign +1: +ch_str / +ch_bra
    are flow OUT of junction). Conservation:
        -ch_in.m_dot + ch_str.m_dot + ch_bra.m_dot = 0
    """
    _, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    imbalance = -m_in + m_str + m_bra
    assert abs(imbalance) < 1e-6, (
        f"mass imbalance {imbalance} kg/s: m_in={m_in}, m_str={m_str}, m_bra={m_bra}"
    )


def test_three_port_converges_with_conservation(three_port_solution):
    """Network converges and conserves mass.

    Directional behaviour at asymmetric pressure BCs depends on the cross-
    coupling: at sufficiently asymmetric collector Pt the junction acts as
    an ejector and the lower-Pt branch may source flow INTO the junction
    (water-jet-pump effect). This is correct physics per Bassett 2001
    Section 3 (and the basis of the K6 < 0 region at low q). We assert
    convergence + mass conservation only; direction is BC-dependent.
    """
    _, sol = three_port_solution
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    # Common inlet always forward at these BCs (common Pt is highest).
    assert m_in > 0.0, f"inlet reversed: m_in={m_in}"
    # Mass conservation at the junction.
    imbalance = m_in - m_str - m_bra
    assert abs(imbalance) < 1e-3, (
        f"mass imbalance {imbalance}: m_in={m_in}, m_str={m_str}, m_bra={m_bra}"
    )


def test_three_port_p_jct_in_physical_range(three_port_solution):
    """P_jct should lie between the lowest collector Pt and the source Pt."""
    # Fixture defaults (see _three_port_net signature).
    Pt_in = 2.1e5
    Pt_str = 2.05e5
    Pt_bra = 2.00e5
    _, sol = three_port_solution
    assert sol["__success__"]

    P_jct = sol["jct.P_jct"]
    # Sanity: P_jct should sit between the lowest collector P and the source.
    # Use ~10% slack: friction can push P_jct outside the Pt envelope.
    lo, hi = min(Pt_str, Pt_bra) * 0.9, Pt_in * 1.1
    assert lo < P_jct < hi, f"P_jct={P_jct} outside physical range [{lo}, {hi}]"


def test_border_carnot_loss_element_in_series_with_junction():
    """Bolt a BorderCarnotLossElement onto the branch port.

    Smoke test: the network still solves with the loss element in series and
    series mass conservation holds along the lateral path (loss element and
    its downstream channel carry the same m_dot).
    """
    Y = _DRY_AIR_Y
    diameter = 0.05
    length = 0.3
    area = math.pi * (diameter / 2.0) ** 2

    # With loss element on branch:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_in", Pt=2.1e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=2.05e5, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_bra", Pt=2.00e5, Tt=300.0, Y=Y))
    net.add_node(MomentumChamberNode("mc_com", area=area))
    net.add_node(MomentumChamberNode("mc_str", area=area))
    net.add_node(MomentumChamberNode("mc_bra", area=area))
    # Intermediate node between BorderCarnot and the branch channel.
    net.add_node(MomentumChamberNode("mc_bra_post", area=area))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_in", "mc_com", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_str", length=length, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        BorderCarnotLossElement(
            "loss_bra", from_node="mc_bra", to_node="mc_bra_post", delta_geom_deg=90.0, area=area
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra",
            "mc_bra_post",
            "pb_bra",
            length=length,
            diameter=diameter,
            regime="compressible",
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
        )
    )

    sol = NetworkSolver(net).solve()
    assert sol["__success__"], f"did not converge: {sol.get('__message__')}"

    # Series mass conservation: ch_bra (downstream channel) carries the same
    # flow as the loss element (no junction between them, no mass-row skip).
    m_loss = sol["loss_bra.m_dot"]
    m_ch_bra = sol["ch_bra.m_dot"]
    assert abs(m_loss - m_ch_bra) < 1e-6, (
        f"loss/channel series mismatch: m_loss={m_loss}, m_ch_bra={m_ch_bra}"
    )
