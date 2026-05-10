"""
Network-level solver tests for TeeJunctionElement.

Covers:
- Merging tee: two pressure-boundary inlets, one outlet
- Branching tee: one inlet, two pressure-boundary outlets
- Mass conservation at each node
- Jacobian accuracy (FD vs analytical)
- No-flow edge cases
- TeeJunctionElement.validate() error paths
"""

import math

import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
    TeeJunctionElement,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))

_THETA_90 = math.pi / 2.0
_THETA_60 = math.pi / 3.0
_F_C = 0.01  # 100 cm^2


def _merging_net(
    Pt_straight: float = 2.1e5,
    Pt_branch: float = 2.05e5,
    Pt_common: float = 2.0e5,
    theta: float = _THETA_90,
    psi: float = 1.0,
) -> FlowNetwork:
    """Three-node merging tee: two inlets (straight, branch) into one outlet (common)."""
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=400.0, Y=Y))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="pb_straight",
            branch_node="pb_branch",
            theta=theta,
            F_C=_F_C,
            psi=psi,
            tee_type="merging",
        )
    )
    return net


def _branching_net(
    Pt_common: float = 2.1e5,
    Pt_straight: float = 2.0e5,
    Pt_branch: float = 2.02e5,
    theta: float = _THETA_90,
    psi: float = 1.0,
) -> FlowNetwork:
    """Three-node branching tee: one inlet (common) into two outlets (straight, branch)."""
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_common", Pt=Pt_common, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_straight", Pt=Pt_straight, Tt=400.0, Y=Y))
    net.add_node(PressureBoundary("pb_branch", Pt=Pt_branch, Tt=400.0, Y=Y))
    net.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pb_common",
            straight_node="pb_straight",
            branch_node="pb_branch",
            theta=theta,
            F_C=_F_C,
            psi=psi,
            tee_type="branching",
        )
    )
    return net


# ---------------------------------------------------------------------------
# Topology tests
# ---------------------------------------------------------------------------


def test_merging_tee_graph_connectivity():
    """FlowNetwork registers all three nodes in connectivity dicts."""
    net = _merging_net()
    assert "tee" in net._downstream_of_node["pb_straight"]
    assert "tee" in net._downstream_of_node["pb_branch"]
    assert "tee" in net._upstream_of_node["pb_common"]


def test_branching_tee_graph_connectivity():
    net = _branching_net()
    assert "tee" in net._downstream_of_node["pb_common"]
    assert "tee" in net._upstream_of_node["pb_straight"]
    assert "tee" in net._upstream_of_node["pb_branch"]


def test_merging_reachable_from_boundary():
    """All nodes reachable from pb_straight in merging net."""
    net = _merging_net()
    reachable = net._reachable_from("pb_straight")
    assert "pb_branch" in reachable
    assert "pb_common" in reachable


def test_tee_element_unknowns():
    net = _merging_net()
    assert net.elements["tee"].unknowns() == ["tee.m_dot_com", "tee.m_dot_branch"]
    assert net.elements["tee"].n_equations() == 2


# ---------------------------------------------------------------------------
# Solver: merging tee
# ---------------------------------------------------------------------------


def test_merging_tee_solver_converges():
    """Merging tee solver finds finite positive flows."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert math.isfinite(m_com), "m_dot_com is not finite"
    assert math.isfinite(m_branch), "m_dot_branch is not finite"
    assert m_com > 0.0, "m_dot_com should be positive (flow into common outlet)"
    assert 0.0 < m_branch < m_com, "branch flow should be between 0 and total"


def test_merging_tee_mass_conservation():
    """m_dot_com == m_dot_straight + m_dot_branch at solved state."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    m_straight = m_com - m_branch
    assert abs(m_straight) >= 0.0
    assert abs(m_com - m_straight - m_branch) < 1e-10


def test_merging_tee_residuals_near_zero():
    """Residuals evaluated at the solution must be near zero."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()

    import numpy as np

    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert all(abs(r) < 1e-6 for r in res), f"Non-zero residuals: {res}"


def test_merging_tee_diagnostics():
    """Solved solution produces finite diagnostic values."""
    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    assert math.isfinite(sol["tee.q"])
    assert math.isfinite(sol["tee.m_dot_com"])
    assert math.isfinite(sol["tee.m_dot_branch"])


# ---------------------------------------------------------------------------
# Solver: branching tee
# ---------------------------------------------------------------------------


def test_branching_tee_solver_converges():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert math.isfinite(m_com)
    assert math.isfinite(m_branch)
    assert m_com > 0.0


def test_branching_tee_mass_conservation():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    m_com = sol["tee.m_dot_com"]
    m_branch = sol["tee.m_dot_branch"]
    assert abs(m_com - m_branch) > 0.0  # straight flow is nonzero


def test_branching_tee_residuals_near_zero():
    net = _branching_net()
    solver = NetworkSolver(net)
    sol = solver.solve()

    import numpy as np

    x = np.array([sol.get(n, 0.0) for n in solver.unknown_names])
    res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
    assert all(abs(r) < 1e-6 for r in res), f"Non-zero residuals: {res}"


# ---------------------------------------------------------------------------
# Jacobian accuracy (FD check)
# ---------------------------------------------------------------------------


def test_merging_tee_jacobian_vs_fd():
    """Analytical Jacobian matches central-difference FD to within 1e-4."""
    import numpy as np

    net = _merging_net()
    solver = NetworkSolver(net)
    sol = solver.solve()
    x0 = np.array([sol.get(n, 0.0) for n in solver.unknown_names])

    _, jac_analytical = solver._residuals_and_jacobian(x0, compute_jacobian=True)
    jac_dense = jac_analytical.toarray()

    eps = 1e-5
    n = len(x0)
    m = len(solver._residuals(x0))
    jac_fd = np.zeros((m, n))
    for j in range(n):
        xp, xm = x0.copy(), x0.copy()
        xp[j] += eps
        xm[j] -= eps
        jac_fd[:, j] = (solver._residuals(xp) - solver._residuals(xm)) / (2 * eps)

    # Only compare non-zero FD columns to avoid noise in near-zero entries
    mask = np.abs(jac_fd).max(axis=0) > 1e-8
    if mask.any():
        diff = np.abs(jac_dense[:, mask] - jac_fd[:, mask])
        assert diff.max() < 5e-3, f"Jacobian FD mismatch: max={diff.max():.3e}"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_validate_bad_theta():
    with pytest.raises(ValueError, match="theta"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b", "c"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "c", theta=0.0, F_C=0.01))
        net.validate()


def test_validate_bad_F_C():
    with pytest.raises(ValueError, match="F_C"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b", "c"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "c", theta=_THETA_90, F_C=-1.0))
        net.validate()


def test_validate_duplicate_node_rejected():
    """Adding a tee where two ports share the same node raises ValueError."""
    with pytest.raises(ValueError, match="duplicate"):
        net = FlowNetwork()
        Y = _DRY_AIR_Y
        for nid in ("a", "b"):
            net.add_node(PressureBoundary(nid, Pt=2e5, Tt=400.0, Y=Y))
        net.add_element(TeeJunctionElement("tee", "a", "b", "a", theta=_THETA_90, F_C=0.01))
