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
