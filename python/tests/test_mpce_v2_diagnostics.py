"""
End-to-end test that ``MPCEv2Element.diagnostics`` emits K_straight,
K_branch (separating) / K11, K12 (joining), and mass_flow_ratio so the
GUI can display them and validation cross-checks can use them.

The values must match a direct call to ``junction_loss_coefficient`` at
the converged operating point -- this is the audit-against-Bassett /
Idelchik / Hager pipeline plan step #4.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from validation.junction.models.mpce_v2_network import MPCEv2Network
from validation.junction.models.mynard2010 import junction_loss_coefficient


def _mynard_K_direct(
    q: float, psi: float, theta_deg: float, alpha: float, flow_dir: str
) -> tuple[float, float]:
    """Direct Mynard K_straight/K_branch (separating) or K11/K12 (joining)."""
    A_str = 1.0
    A_bra = 1.0 / psi
    A_com = 1.0
    m_str = (1.0 - q) * 1.0
    m_bra = q * 1.0
    m_com = 1.0
    if flow_dir == "branch":  # separating
        # com supplier (+U), str+bra collectors (-U)
        U = np.array([m_com / A_com, -m_str / A_str, -m_bra / A_bra])
        theta = np.array([math.pi, 0.0, math.radians(theta_deg)])
    else:  # joining
        # str+bra suppliers (+U), com collector (-U)
        U = np.array([m_str / A_str, m_bra / A_bra, -m_com / A_com])
        theta = np.array([0.0, math.radians(theta_deg), math.pi])
    A = np.array([A_com, A_str, A_bra]) if flow_dir == "branch" else np.array([A_str, A_bra, A_com])
    r = junction_loss_coefficient(U, A, theta, joining_etransfer_alpha=alpha)
    assert r.K is not None and len(r.K) == 2
    return float(r.K[0]), float(r.K[1])


def test_separating_emits_K_straight_K_branch_and_matches_direct_mynard():
    """The K values the element emits must match a direct Mynard call at
    the same (q, psi, theta)."""
    net = MPCEv2Network(joining_etransfer_alpha=0.0)  # alpha=0 for separating test
    # imposed_q at theta=90, psi=1, q=0.3 — a sane separating operating point
    q, psi, theta_deg = 0.3, 1.0, 90.0
    r = net.evaluate_network(
        "bassett2001", "K6", q, psi, math.radians(theta_deg), topology="imposed_q"
    )
    assert r.converged, f"expected converged solve, got: {r.message!r}"

    # Direct Mynard at the same operating point
    K_str_ref, K_bra_ref = _mynard_K_direct(q, psi, theta_deg, 0.0, "branch")

    # solve_and_extract uses K = (Pt_com - Pt_branch)/q_dyn_com convention,
    # which equals Mynard's K for the collector ports. Tolerance: ~5% to
    # account for Pt at port-MCN vs Pt at far-field after the channel.
    assert r.K_straight is not None and r.K_lateral is not None
    assert r.K_straight == pytest.approx(K_str_ref, abs=0.05)
    assert r.K_lateral == pytest.approx(K_bra_ref, abs=0.05)


class _State:
    """Minimal stand-in for NetworkMixtureState (Pt + density used)."""

    def __init__(self, Pt: float, rho: float = 1.16) -> None:
        self.Pt = Pt
        self.P = Pt - 1000.0  # arbitrary; not used by diagnostics K logic
        self.T = 300.0
        self._rho = rho

    def density(self) -> float:
        return self._rho


def test_diagnostics_element_level_separating_emits_K_named_aliases():
    """Unit test: call diagnostics() directly on a separating element and
    check the new fields are present + match direct Mynard call.

    This verifies the new emission path independently from solve_and_extract."""
    from combaero.network.mpce_v2_element import MPCEv2Element

    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_com"],
        outlet_nodes=["port_str", "port_bra"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="branch",
        joining_etransfer_alpha=0.0,  # disable correction for cleaner mapping
    )
    # Mark _port_element_ids so the parent class doesn't error; values
    # don't affect the K math, only the dict keys downstream.
    e._port_element_ids = ["lc0", "lc1", "lc2"]

    # Operating point: q=0.3, psi=1, theta=90, m_com=0.1 kg/s
    m_com = 0.1
    q = 0.3
    port_mdots = [-m_com, +(1.0 - q) * m_com, +q * m_com]  # com inlet, str+bra outlets
    states = [_State(1e5), _State(1e5), _State(1e5)]
    diag = e.diagnostics(states, Pt_jct=1e5, port_mdots=port_mdots)

    # Expected from direct Mynard
    K_str_ref, K_bra_ref = _mynard_K_direct(q, 1.0, 90.0, 0.0, "branch")

    assert "K_straight" in diag
    assert "K_branch" in diag
    assert "mass_flow_ratio" in diag
    assert "Pt_jct" in diag
    assert diag["K_straight"] == pytest.approx(K_str_ref, abs=1e-6)
    assert diag["K_branch"] == pytest.approx(K_bra_ref, abs=1e-6)
    assert diag["mass_flow_ratio"] == pytest.approx(q, abs=1e-9)
    # Topology-aware mdot aliases let runner.py expose ElementResult.m_dot
    # as the junction throughput instead of the meaningless 0.0 fallback.
    assert diag["m_dot_com"] == pytest.approx(m_com, abs=1e-12)
    assert diag["m_dot_branch"] == pytest.approx(q * m_com, abs=1e-12)
    assert diag["m_dot_straight"] == pytest.approx((1.0 - q) * m_com, abs=1e-12)


def test_diagnostics_element_level_joining_emits_K11_K12():
    """Same shape for joining: K11/K12 + mass_flow_ratio."""
    from combaero.network.mpce_v2_element import MPCEv2Element

    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_str", "port_bra"],
        outlet_nodes=["port_com"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
        joining_etransfer_alpha=0.2,
    )
    e._port_element_ids = ["lc0", "lc1", "lc2"]

    m_com = 0.1
    q = 0.4
    # str+bra inlets (negative mdot), com outlet (positive)
    port_mdots = [-(1.0 - q) * m_com, -q * m_com, +m_com]
    states = [_State(1e5), _State(1e5), _State(1e5)]
    diag = e.diagnostics(states, Pt_jct=1e5, port_mdots=port_mdots)

    K11_ref, K12_ref = _mynard_K_direct(q, 1.0, 90.0, 0.2, "merge")
    assert "K11" in diag
    assert "K12" in diag
    assert "mass_flow_ratio" in diag
    assert diag["K11"] == pytest.approx(K11_ref, abs=1e-6)
    assert diag["K12"] == pytest.approx(K12_ref, abs=1e-6)
    assert diag["mass_flow_ratio"] == pytest.approx(q, abs=1e-9)
    assert diag["m_dot_com"] == pytest.approx(m_com, abs=1e-12)
    assert diag["m_dot_branch"] == pytest.approx(q * m_com, abs=1e-12)
    assert diag["m_dot_straight"] == pytest.approx((1.0 - q) * m_com, abs=1e-12)


def test_diagnostics_without_port_mdots_returns_parent_fields_only():
    """Without port_mdots (legacy callers), diagnostics returns the parent
    class fields only — no K, no aliases. Ensures backward compatibility."""
    from combaero.network.mpce_v2_element import MPCEv2Element

    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_com"],
        outlet_nodes=["port_str", "port_bra"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="branch",
    )
    e._port_element_ids = ["lc0", "lc1", "lc2"]
    states = [_State(1e5), _State(1e5), _State(1e5)]
    diag = e.diagnostics(states, Pt_jct=1e5, port_mdots=None)
    assert "P_jct" in diag
    assert "n_ports" in diag
    assert "K_straight" not in diag
    assert "K11" not in diag


def test_joining_emits_K11_K12_and_matches_direct_mynard():
    """Joining flow: K11 (str->com), K12 (bra->com) match direct Mynard
    + the alpha=0.2 correction."""
    net = MPCEv2Network()  # default alpha=0.2
    q, psi, theta_deg = 0.5, 1.0, 90.0
    r = net.evaluate_network(
        "bassett2001", "K12", q, psi, math.radians(theta_deg), topology="imposed_q"
    )
    assert r.converged, f"expected converged solve, got: {r.message!r}"

    K11_ref, K12_ref = _mynard_K_direct(q, psi, theta_deg, 0.2, "merge")

    assert r.K_straight is not None and r.K_lateral is not None
    assert r.K_straight == pytest.approx(K11_ref, abs=0.05)
    assert r.K_lateral == pytest.approx(K12_ref, abs=0.05)
