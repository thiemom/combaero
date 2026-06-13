"""
Tests for the constrained-topology ``flow_direction`` field on MPCEv2Element.

The element accepts ``flow_direction = "merge" | "branch"`` and raises
``ValueError`` at solve time if the observed flow pattern at residual
evaluation disagrees with the declaration.
"""

from __future__ import annotations

import math

import pytest

from combaero.network.mpce_v2_element import MPCEv2Element


def _branch_element() -> MPCEv2Element:
    """Canonical 3-port separating T at theta=90deg, equal areas."""
    return MPCEv2Element(
        id="jct",
        inlet_nodes=["port_com"],
        outlet_nodes=["port_str", "port_bra"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="branch",
    )


def _merge_element() -> MPCEv2Element:
    """Canonical 3-port joining T at theta=90deg, equal areas."""
    return MPCEv2Element(
        id="jct",
        inlet_nodes=["port_str", "port_bra"],
        outlet_nodes=["port_com"],
        inlet_angles_deg=[0.0, 90.0],
        outlet_angles_deg=[0.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="merge",
    )


class _State:
    """Minimal stand-in for NetworkMixtureState (only Pt + density used)."""

    def __init__(self, Pt: float, rho: float = 1.16) -> None:
        self.Pt = Pt
        self._rho = rho

    def density(self) -> float:
        return self._rho


def test_invalid_flow_direction_rejected_at_construction() -> None:
    with pytest.raises(ValueError, match="flow_direction"):
        MPCEv2Element(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            port_areas=[0.01, 0.01, 0.01],
            flow_direction="auto",  # type: ignore[arg-type]
        )


def test_branch_with_separating_flow_succeeds() -> None:
    """branch + (1 supplier, 2 collectors) is the canonical separating case."""
    e = _branch_element()
    states = [_State(1e5), _State(1e5), _State(1e5)]
    port_mdots = [-0.1, +0.05, +0.05]  # com in, str+bra out
    residuals, jac = e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)
    assert len(residuals) == 4  # 3 port residuals + mass conservation
    assert math.isfinite(residuals[0])


def test_branch_with_joining_flow_raises() -> None:
    """branch element with joining flow at runtime must refuse."""
    e = _branch_element()
    states = [_State(1e5), _State(1e5), _State(1e5)]
    port_mdots = [+0.1, -0.05, -0.05]  # com out, str+bra in -- joining
    with pytest.raises(ValueError, match="declared flow_direction='branch'"):
        e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)


def test_merge_with_joining_flow_succeeds() -> None:
    """merge + (2 suppliers, 1 collector) is the canonical joining case."""
    e = _merge_element()
    states = [_State(1e5), _State(1e5), _State(1e5)]
    port_mdots = [-0.05, -0.05, +0.1]  # str+bra in, com out
    residuals, jac = e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)
    assert len(residuals) == 4
    assert math.isfinite(residuals[0])


def test_merge_with_separating_flow_raises() -> None:
    """merge element with separating flow at runtime must refuse."""
    e = _merge_element()
    states = [_State(1e5), _State(1e5), _State(1e5)]
    port_mdots = [+0.05, +0.05, -0.1]  # str+bra out, com in -- separating
    with pytest.raises(ValueError, match="declared flow_direction='merge'"):
        e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)


def test_degenerate_flow_does_not_trip_constraint() -> None:
    """All-zero mdots fall through to continuity residual and do not raise."""
    e = _branch_element()
    states = [_State(1e5), _State(1e5), _State(1e5)]
    port_mdots = [0.0, 0.0, 0.0]
    residuals, jac = e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)
    assert len(residuals) == 4


def test_soft_mode_does_not_raise_on_wrong_direction() -> None:
    """strict=False replaces the raise with a soft-barrier residual."""
    e = MPCEv2Element(
        id="jct",
        inlet_nodes=["port_com"],
        outlet_nodes=["port_str", "port_bra"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01, 0.01, 0.01],
        flow_direction="branch",
        strict=False,
    )
    states = [_State(1e5), _State(1e5), _State(1e5)]
    # Wrong direction for branch: 2 suppliers + 1 collector
    port_mdots = [+0.1, -0.05, -0.05]
    residuals, jac = e.residuals(states, Pt_jct=1e5, port_mdots=port_mdots)
    assert len(residuals) == 4
    assert all(math.isfinite(r) for r in residuals)


def test_verify_solution_consistent_accepts_physical() -> None:
    """outer_mdot > 0 for all connecting elements -> physical."""
    e = _branch_element()
    e._port_element_ids = ["lc_a", "lc_b", "lc_c"]
    sol = {"lc_a.m_dot": 0.1, "lc_b.m_dot": 0.05, "lc_c.m_dot": 0.05}
    assert e.verify_solution_consistent(sol) is True


def test_verify_solution_consistent_rejects_wrong_basin() -> None:
    """Any negative outer_mdot -> reject (soft mode landed off-physical)."""
    e = _branch_element()
    e._port_element_ids = ["lc_a", "lc_b", "lc_c"]
    sol = {"lc_a.m_dot": 0.1, "lc_b.m_dot": -0.05, "lc_c.m_dot": 0.05}
    assert e.verify_solution_consistent(sol) is False


def test_verify_solution_consistent_handles_missing_keys() -> None:
    """Connecting elements without m_dot unknowns are skipped, not rejected."""
    e = _branch_element()
    e._port_element_ids = ["pb_in", "lc_b", "lc_c"]
    sol = {"lc_b.m_dot": 0.05, "lc_c.m_dot": 0.05}  # pb_in has no m_dot key
    assert e.verify_solution_consistent(sol) is True
