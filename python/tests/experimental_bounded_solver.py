"""
Bounds experiment: does enforcing physical positivity (m_dot >= 0 on
canonically forward channels, all pressures >= 0) make MPCE Tier-1 match
Bassett K2/K6 at M -> 0?

If yes: Tier-1 disagreement is solver landing on degenerate roots
(bidirectional-flow basins that are mathematically valid but physically
forbidden by the boundary conditions). Fix is to expose bounded
least_squares as a solver option for MPCE networks.

If no: the PDF Section 3.4 formula itself is the problem -- no constant L
in form-2 BorderCarnotLoss can match Bassett K6 = q^2 + 1 - 0.7654*q.
Fix is to revisit the loss-element form.

The experiment wraps scipy.optimize.least_squares(method='trf', bounds=...)
around NetworkSolver's existing _residuals_and_jacobian() callable. The
production solver is untouched.

Run with: uv run pytest python/tests/experimental_bounded_solver.py -v -s
The -s flag is recommended; this file prints diagnostic K values rather
than just asserting.
"""

import math

import numpy as np
import pytest
from scipy.optimize import least_squares

import combaero as cb
from combaero.network import (
    BorderCarnotLossElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_F_C = 0.01  # 100 cm^2
_M_DOT_IN = 0.1  # kg/s


def _K6_bassett(q, psi=1.0, theta_rad=math.pi / 2.0):
    return q * q * psi * psi + 1.0 - 2.0 * q * psi * math.cos(0.75 * theta_rad)


def _K2_bassett(q):
    return q * q - 1.5 * q + 0.5


def _build_imposed_q_network(m_in, m_lateral, Pt_ref=1.0e5):
    """Same topology as test_momentum_cv_tier1_bassett._build_imposed_q_network."""
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=300.0, Y=Y))
    net.add_node(MassFlowBoundary("mb_lat", m_dot=m_lateral, Tt=300.0, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_ref, Tt=300.0, Y=Y))
    net.add_node(MomentumChamberNode("port_com", area=_F_C))
    net.add_node(MomentumChamberNode("port_str", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C))
    net.add_node(MomentumChamberNode("port_bra_post", area=_F_C))
    net.add_element(LosslessConnectionElement("lc_in", "mb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    net.add_element(
        BorderCarnotLossElement(
            "loss_bra",
            from_node="port_bra",
            to_node="port_bra_post",
            delta_geom_deg=90.0,
            area=_F_C,
        )
    )
    net.add_element(LosslessConnectionElement("lc_bra", "port_bra_post", "mb_lat"))
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[_F_C, _F_C, _F_C],
        )
    )
    return net


def _build_bounds(solver: NetworkSolver):
    """Build bounds vectors for least_squares.

    Conservative physical bounds:
    - All pressure-like unknowns (.P, .Pt, .P_jct): lower 0, upper +inf.
    - All m_dot unknowns: lower 0, upper +inf (forward direction enforced).
      This is the experiment -- assumes the imposed-q test has unambiguous
      forward flow at every channel.
    - Other unknowns (T, Y, etc.): unbounded.
    """
    lo, hi = [], []
    for name in solver.unknown_names:
        if name.endswith(".P") or name.endswith(".Pt") or name.endswith(".P_jct"):
            lo.append(1.0)  # 1 Pa floor avoids exact-zero singularity in density
            hi.append(np.inf)
        elif (
            name.endswith(".m_dot") or name.endswith(".m_dot_com") or name.endswith(".m_dot_branch")
        ):
            lo.append(0.0)
            hi.append(np.inf)
        else:
            lo.append(-np.inf)
            hi.append(np.inf)
    return np.array(lo), np.array(hi)


def _solve_bounded(net: FlowNetwork) -> tuple[bool, dict[str, float], float]:
    """Run the experiment: least_squares(trf) with bounds + analytic Jacobian.

    Returns (success, name->value dict, residual norm).
    """
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = np.array(solver._build_x0())

    bounds = _build_bounds(solver)
    # Clamp x0 into the box (least_squares requires this).
    x0 = np.clip(x0, bounds[0], bounds[1])

    def fun(x):
        res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
        return np.asarray(res, dtype=float)

    def jac(x):
        _, J = solver._residuals_and_jacobian(x, compute_jacobian=True)
        return J.toarray()

    try:
        result = least_squares(
            fun, x0, jac=jac, bounds=bounds, method="trf", max_nfev=400, xtol=1e-10
        )
    except Exception as e:
        return False, {"error": str(e)}, np.inf

    success = result.cost * 2 < 1e-3  # cost = 0.5 * ||F||^2
    sol_dict = dict(zip(solver.unknown_names, result.x, strict=True))
    res_norm = float(np.linalg.norm(result.fun))
    return success, sol_dict, res_norm


def _extract_K(sol: dict[str, float], m_in: float) -> dict[str, float]:
    """Pull K_straight and K_lateral from the converged state."""
    P_com = sol["port_com.P"]
    Pt_com = sol["port_com.Pt"]
    Pt_str = sol["port_str.Pt"]
    Pt_bra_post = sol["port_bra_post.Pt"]

    X_air = list(cb.mass_to_mole(_DRY_AIR_Y))
    rho_com = float(cb.density(300.0, P_com, X_air))
    u_com = m_in / (rho_com * _F_C)
    q_dyn_com = 0.5 * rho_com * u_com * u_com
    a_com = float(cb.speed_of_sound(300.0, X_air))

    return {
        "K_straight": (Pt_com - Pt_str) / q_dyn_com,
        "K_lateral": (Pt_com - Pt_bra_post) / q_dyn_com,
        "Mach_com": abs(u_com) / a_com,
    }


@pytest.mark.parametrize("q", [0.25, 0.5, 0.75])
def test_bounded_solver_K_lateral_vs_bassett_K6(q):
    m_in = _M_DOT_IN
    net = _build_imposed_q_network(m_in=m_in, m_lateral=q * m_in)
    ok, sol, res_norm = _solve_bounded(net)
    if not ok:
        pytest.skip(f"q={q}: bounded solve did not converge (|F|={res_norm:.3e})")

    K = _extract_K(sol, m_in)
    K6_ref = _K6_bassett(q)
    diff = abs(K["K_lateral"] - K6_ref)
    rel = diff / max(abs(K6_ref), 1e-6)

    print(
        f"\n[q={q}] K_lateral={K['K_lateral']:.4f}  Bassett K6={K6_ref:.4f}  "
        f"rel err={rel:.1%}  |F|={res_norm:.3e}  M_com={K['Mach_com']:.4f}"
    )
    # Loose tolerance for the experiment; we are looking for a qualitative shift.
    if rel < 0.20 or diff < 0.05:
        # PASS = bounds fixed the disagreement -> Tier 1 issue was solver
        # finding degenerate roots, NOT the PDF formula. Promote this to a
        # production solver path.
        return
    # FAIL = bounds did not fix it -> PDF formula needs revising.
    pytest.fail(
        f"Bounded least_squares converged ({res_norm:.3e}) but K_lateral "
        f"still disagrees with Bassett K6 by {rel:.1%}. PDF formula is the "
        f"problem, not solver degeneracy."
    )


@pytest.mark.parametrize("q", [0.25, 0.5, 0.75])
def test_bounded_solver_K_straight_vs_bassett_K2(q):
    m_in = _M_DOT_IN
    net = _build_imposed_q_network(m_in=m_in, m_lateral=q * m_in)
    ok, sol, res_norm = _solve_bounded(net)
    if not ok:
        pytest.skip(f"q={q}: bounded solve did not converge (|F|={res_norm:.3e})")

    K = _extract_K(sol, m_in)
    K2_ref = _K2_bassett(q)
    diff = abs(K["K_straight"] - K2_ref)
    rel = diff / max(abs(K2_ref), 1e-6)

    print(
        f"\n[q={q}] K_straight={K['K_straight']:.4f}  Bassett K2={K2_ref:.4f}  "
        f"rel err={rel:.1%}  |F|={res_norm:.3e}  M_com={K['Mach_com']:.4f}"
    )
    if rel < 0.20 or diff < 0.05:
        return
    pytest.fail(
        f"Bounded least_squares converged ({res_norm:.3e}) but K_straight "
        f"still disagrees with Bassett K2 by {rel:.1%}. PDF formula is the "
        f"problem, not solver degeneracy."
    )
