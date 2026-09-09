"""Degenerate iterates in MultiPortChamberElement: what the solver is handed.

Three faces of one question, measured on the scorecard before any of this was
written (issue #271, step 2):

* the pre-check fallback returned a lossless residual with an EMPTY Jacobian
  on 145 iterates from 26 solves -- and 0 of those 26 converged;
* a port at exactly zero flow was excluded from MultiPortChamberElement masks but classified
  as a supplier by Mynard, so K came back mis-shaped and was silently zeroed
  -- 384 times per scorecard, a lossless junction at initialization;
* a 4-port junction returned finite residuals with an all-zero loss Jacobian.

These tests pin the replacements. Ground truth is what a Newton solver needs
-- a residual with the Jacobian that belongs to it -- and, for the dead-port
case, the FlowRatio -> 0 limit proven in step 1.4.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb
from combaero.network import mpce_element as v2
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_element import ConstantKTeeElement, MultiPortChamberElement

_Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _element(flow_direction: str = "branch", strict: bool = False) -> MultiPortChamberElement:
    element = MultiPortChamberElement.__new__(MultiPortChamberElement)
    element.id = "jct"
    element.N = 3
    element.port_nodes = ["p0", "p1", "p2"]
    element.port_areas = [0.01, 0.01, 0.01]
    element.port_angles_deg = [180.0, 0.0, 90.0]
    element.flow_direction = flow_direction
    element._port_signs = [-1.0, 1.0, 1.0] if flow_direction == "branch" else [-1.0, -1.0, 1.0]
    element._port_element_ids = ["e0", "e1", "e2"]
    element.strict = strict
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    element.eta_scale = 1.0
    return element


def _states(pt=(100_300.0, 100_100.0, 100_050.0)):
    return [
        NetworkMixtureState(P=1.0e5, Pt=pt[i], T=300.0, Tt=300.5, m_dot=0.0, Y=_Y) for i in range(3)
    ]


# ---------------------------------------------------------------------------
# (a) the pre-check fallback: never an empty Jacobian
# ---------------------------------------------------------------------------


def test_cold_start_all_zero_iterate_gets_continuity_with_its_jacobian():
    """Every port ~zero is a cold start; the continuity residual is right,
    but it has a Jacobian and the solver must receive it."""
    residuals, jac = _element().residuals(_states(), 100_300.0, [0.0, 0.0, 0.0])

    assert len(residuals) == 4
    for i in range(3):
        assert jac[i][f"p{i}.Pt"] == 1.0
        assert jac[i]["jct.P_jct"] == -1.0
    assert jac[3] == {"e0.m_dot": -1.0, "e1.m_dot": 1.0, "e2.m_dot": 1.0}


@pytest.mark.parametrize(
    "label, mdots",
    [
        ("every port flowing out", [0.0021, 0.0011, 0.0048]),
        ("every port flowing in", [-0.05, -0.02, -0.03]),
    ],
)
def test_same_sign_iterate_is_a_wrong_direction_state_with_a_jacobian(label, mdots):
    """All ports flowing the same way violates mass conservation and is a
    wrong-direction state for at least one port. It must reach the soft
    barrier -- which carries the derivative that pulls the offending port back
    toward zero -- not a lossless residual with nothing for Newton to use."""
    residuals, jac = _element().residuals(_states(), 100_300.0, mdots)

    assert len(residuals) == 4
    assert all(i in jac for i in range(4)), f"{label}: rows missing from the Jacobian"
    penalised = [i for i in range(3) if any(k.endswith(".m_dot") for k in jac[i])]
    assert penalised, f"{label}: no port is being pulled back -- this is the old lossless fallback"


def test_same_sign_iterate_raises_when_strict():
    """Under strict=True a wrong-direction state raises, as it always did for
    a single wrong port; the all-same-sign case is no longer exempt."""
    with pytest.raises(ValueError, match="wrong flow direction"):
        _element(strict=True).residuals(_states(), 100_300.0, [0.0021, 0.0011, 0.0048])


# ---------------------------------------------------------------------------
# (b) a port at exactly zero: both classifications must agree
# ---------------------------------------------------------------------------


def _spy_outer_mdot(monkeypatch, element, states, pt_jct, mdots):
    """The velocities the CLOSURE ends up seeing, observed at the kernel call.

    The snapped state used to be checked by spying on the Python closure. Since
    the whole-element (f, J) moved to C++ (#271 step 4.4) the residual path
    hands the kernel mass flows and the kernel derives U itself, so the same
    quantity is read back from the flows: U_i = -port_sign_i * outer_i /
    (rho_i A_i), which is sign-equivalent to -port_sign_i * outer_i.
    """
    seen = {}
    real = v2._core.mpce_residuals_and_jacobian

    def spy(p_static, p_total, rho, drho_dp, outer_mdot, pt, geom):
        seen["U"] = np.array(
            [
                -float(geom.port_sign[i])
                * float(outer_mdot[i])
                / (float(rho[i]) * float(geom.area[i]))
                for i in range(len(outer_mdot))
            ]
        )
        return real(p_static, p_total, rho, drho_dp, outer_mdot, pt, geom)

    monkeypatch.setattr(v2._core, "mpce_residuals_and_jacobian", spy)
    element.residuals(states, pt_jct, list(mdots))
    assert "U" in seen, "the kernel was never reached"
    return seen["U"]


def test_exact_zero_outlet_is_classified_as_a_dead_collector(monkeypatch):
    """[-0.1, -0.0, 0.125]: the straight outlet at -0.0 (an initial guess).

    Mynard's `Q < 0` puts -0.0 with the suppliers; MultiPortChamberElement excluded it. The
    closure must now be handed a tiny flow in the port's DECLARED direction,
    so it sees one supplier and two collectors.
    """
    U = _spy_outer_mdot(monkeypatch, _element(), _states(), 100_300.0, [-0.1, -0.0, 0.125])

    assert U[1] < 0.0, "the dead outlet must reach the closure as a collector"
    assert (U > 0).sum() == 1 and (U < 0).sum() == 2


def test_exact_zero_outlet_carries_the_dead_branch_loss():
    """FlowRatio -> 0 gives K -> 1 for a dead outlet (step 1.4). Its residual
    row must therefore carry a loss term; it used to be pure continuity."""
    states = _states()
    residuals, _ = _element().residuals(states, 100_300.0, [-0.1, -0.0, 0.125])

    pure_continuity = states[1].Pt - 100_300.0
    assert residuals[1] != pytest.approx(pure_continuity), (
        "row 1 is bare continuity: K for the dead outlet is still being zeroed"
    )


def test_exact_zero_inlet_in_a_merge_is_a_dead_supplier(monkeypatch):
    U = _spy_outer_mdot(monkeypatch, _element("merge"), _states(), 100_300.0, [-0.1, 0.0, 0.1])

    assert U[1] > 0.0, "a dead inlet must reach the closure as a supplier"


# ---------------------------------------------------------------------------
# (c) N > 3 is rejected where it can be seen
# ---------------------------------------------------------------------------


def _four_port(cls):
    return cls(
        id="jct",
        inlet_nodes=["a"],
        outlet_nodes=["b", "c", "d"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0, 45.0],
        port_areas=[0.01] * 4,
    )


@pytest.mark.parametrize("cls", [MultiPortChamberElement, ConstantKTeeElement])
def test_four_port_junction_is_rejected_at_construction(cls):
    """Mynard's K conversion is defined for three branches only. A 4-port
    junction used to return finite residuals with an all-zero loss Jacobian,
    silently. It is now refused with a message that says why."""
    with pytest.raises(ValueError, match="3-branch"):
        _four_port(cls)


def test_three_port_junction_is_still_accepted():
    element = MultiPortChamberElement(
        id="jct",
        inlet_nodes=["a"],
        outlet_nodes=["b", "c"],
        inlet_angles_deg=[0.0],
        outlet_angles_deg=[0.0, 90.0],
        port_areas=[0.01] * 3,
    )
    assert element.N == 3
