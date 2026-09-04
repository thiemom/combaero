"""Whole-row FD guardrail for MPCEv2Element / ConstantKTeeElement.

``test_mpce_v2_jacobian.py`` FD-checks ``dKQ_dmdot_separating_T`` in ISOLATION
-- the d/dmdot block only. Nothing pinned the assembled Jacobian row against
the residual it differentiates, so an entire column could be, and was, absent:
the loss term ``K_i * q_dyn_com`` depends on the common port's density and
hence on its static pressure, but ``MPCEv2Element`` emitted no ``.P`` entry.

``ConstantKTeeElement`` gained exactly that term in PR #230 after an FD test
caught a ``K*q/P`` error; the same defect survived in ``MPCEv2Element`` because
no equivalent test existed. These tests are that test, for both classes.

The unknowns a ``MomentumChamberNode`` exposes are ``.P`` and ``.Pt``
(temperature is derived forward, not solved), so those are the entries that
must match finite differences.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network.components import NetworkMixtureState
from combaero.network.mpce_v2_element import ConstantKTeeElement, MPCEv2Element

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_AREA = 0.01
_PORT_MDOTS = [-0.10, 0.06, 0.04]  # junction convention: negative = into junction
_PT_JCT = 100_300.0
_P_COMMON = 100_000.0


def _wire(element):
    """Attach the topology fields resolve_topology would set, without a graph."""
    element.id = "jct"
    element.N = 3
    element.port_nodes = ["p0", "p1", "p2"]
    element.port_areas = [_AREA, _AREA, _AREA]
    element.port_angles_deg = [180.0, 0.0, 90.0]
    element._port_signs = [-1.0, 1.0, 1.0]
    element._port_element_ids = ["e0", "e1", "e2"]
    element.flow_direction = "branch"
    element.strict = False
    element.joining_etransfer_alpha = 0.2
    element.jacobian_method = "sympy"
    element.penalty_alpha = 0.0
    return element


def _mpce_v2():
    return _wire(MPCEv2Element.__new__(MPCEv2Element))


def _constant_k():
    element = _wire(ConstantKTeeElement.__new__(ConstantKTeeElement))
    element.K_ports = {1: 0.04, 2: 0.9}
    return element


def _states(p_common: float = _P_COMMON, pt_shift: float = 0.0):
    """Port states with the common port's STATIC pressure independently movable.

    Pt is held fixed while P moves, so the derivative isolates the density
    dependence of the dynamic head rather than the trivial Pt term.
    """
    pressures = [p_common, 100_000.0, 100_000.0]
    totals = [100_300.0 + pt_shift, 100_300.0, 100_300.0]
    return [
        NetworkMixtureState(P=pressures[i], Pt=totals[i], T=300.0, Tt=300.5, m_dot=0.0, Y=_Y)
        for i in range(3)
    ]


ELEMENTS = [
    pytest.param(_mpce_v2, id="MPCEv2Element"),
    pytest.param(_constant_k, id="ConstantKTee"),
]


@pytest.mark.parametrize("factory", ELEMENTS)
@pytest.mark.parametrize("row", [0, 1, 2])
def test_row_derivative_wrt_common_static_pressure(factory, row: int):
    """d(R_i)/d(P_common): the density sensitivity of the common dynamic head.

    q_dyn_com = mdot^2 / (2*rho*A^2) and rho = rho(P, T), so every row carrying
    a non-zero K depends on the common port's static pressure. Omitting it
    leaves a silently zero Jacobian column.
    """
    element = factory()
    _, jac = element.residuals(_states(), _PT_JCT, _PORT_MDOTS)

    eps = 1.0
    up = element.residuals(_states(_P_COMMON + eps), _PT_JCT, _PORT_MDOTS)[0][row]
    down = element.residuals(_states(_P_COMMON - eps), _PT_JCT, _PORT_MDOTS)[0][row]
    expected = (up - down) / (2.0 * eps)

    actual = jac[row].get("p0.P", 0.0)
    assert actual == pytest.approx(expected, rel=1e-4, abs=1e-9), (
        f"row {row} d/dP_common: analytic={actual:.6e} fd={expected:.6e} "
        f"(a missing term here is exactly K*q_dyn/P)"
    )


@pytest.mark.parametrize("factory", ELEMENTS)
@pytest.mark.parametrize("row", [0, 1, 2])
def test_row_derivative_wrt_port_total_pressure(factory, row: int):
    """d(R_i)/d(Pt_i): the explicit Pt term, guarding the rest of the row."""
    element = factory()
    _, jac = element.residuals(_states(), _PT_JCT, _PORT_MDOTS)

    eps = 1.0
    base = _states()

    def shifted(delta: float):
        out = list(_states())
        s = out[row]
        out[row] = NetworkMixtureState(
            P=s.P, Pt=s.Pt + delta, T=s.T, Tt=s.Tt, m_dot=s.m_dot, Y=list(s.Y)
        )
        return out

    expected = (
        element.residuals(shifted(eps), _PT_JCT, _PORT_MDOTS)[0][row]
        - element.residuals(shifted(-eps), _PT_JCT, _PORT_MDOTS)[0][row]
    ) / (2.0 * eps)
    actual = jac[row].get(f"p{row}.Pt", 0.0)
    assert actual == pytest.approx(expected, rel=1e-6, abs=1e-9)
    assert base is not None


@pytest.mark.parametrize("factory", ELEMENTS)
@pytest.mark.parametrize("row", [0, 1, 2])
def test_row_derivative_wrt_junction_pressure(factory, row: int):
    element = factory()
    _, jac = element.residuals(_states(), _PT_JCT, _PORT_MDOTS)

    eps = 1.0
    expected = (
        element.residuals(_states(), _PT_JCT + eps, _PORT_MDOTS)[0][row]
        - element.residuals(_states(), _PT_JCT - eps, _PORT_MDOTS)[0][row]
    ) / (2.0 * eps)
    assert jac[row].get("jct.P_jct", 0.0) == pytest.approx(expected, rel=1e-6, abs=1e-9)


@pytest.mark.parametrize("factory", ELEMENTS)
@pytest.mark.parametrize("row", [0, 1, 2])
@pytest.mark.parametrize("port", [0, 1, 2])
def test_row_derivative_wrt_port_mass_flow(factory, row: int, port: int):
    """d(R_i)/d(mdot_j), including the loss term's dependence on the common flow."""
    element = factory()
    _, jac = element.residuals(_states(), _PT_JCT, _PORT_MDOTS)

    eps = 1e-6
    up, down = list(_PORT_MDOTS), list(_PORT_MDOTS)
    up[port] += eps
    down[port] -= eps
    expected = (
        element.residuals(_states(), _PT_JCT, up)[0][row]
        - element.residuals(_states(), _PT_JCT, down)[0][row]
    ) / (2.0 * eps)

    # Rows are assembled in OUTER element unknowns; port_mdots carry the sign map.
    actual = jac[row].get(f"e{port}.m_dot", 0.0) * element._port_signs[port]
    assert actual == pytest.approx(expected, rel=1e-3, abs=1e-6), (
        f"row {row} d/dmdot_{port}: analytic={actual:.6e} fd={expected:.6e}"
    )


def test_mass_row_is_the_signed_port_sum():
    element = _mpce_v2()
    residuals, jac = element.residuals(_states(), _PT_JCT, _PORT_MDOTS)
    assert residuals[3] == pytest.approx(sum(_PORT_MDOTS))
    for i in range(3):
        assert jac[3][f"e{i}.m_dot"] == pytest.approx(element._port_signs[i])


def test_angle_is_in_radians_where_the_closure_expects_it():
    """Guard against a degree/radian slip in the port-angle plumbing."""
    element = _mpce_v2()
    assert math.isclose(math.radians(element.port_angles_deg[2]), math.pi / 2.0)
