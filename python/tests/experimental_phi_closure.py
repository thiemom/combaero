"""
Phi flow-fraction closure prototype (exp/phi-closure-prototype branch).

Reparametrizes MPCE junction leg mass flows as fractions of the junction
throughflow: for a junction with common port c and K non-common legs,

    m_leg_k = phi_k * m_common          (element frame, k = 0..K-1)
    phi_{K-1} = 1 - sum_{k<K-1} phi_k

eliminating K-1 leg m_dots from the unknown vector and dropping the
junction sum-mass row (identically satisfied). With trf bounds
phi in [0, 1] and m_common >= 0 the feasible box IS the physical
orthant: the 2^N sign basins and the strict/soft barrier seam cannot
be entered at all (the soft penalty max(0, -e_i * mdot_i)^2 is zero
everywhere inside the box).

The production solver is untouched: the transform wraps
NetworkSolver._residuals_and_jacobian exactly like
experimental_bounded_solver.py wraps it, and scipy least_squares(trf)
solves the reduced system F_red(y) = drop_rows(F(T(y))).

Run smoke tests explicitly (file is not collected by default):
    uv run pytest python/tests/experimental_phi_closure.py -v -s
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import MultiPortChamberElement


class PhiTransformError(RuntimeError):
    """Raised when the network topology is outside the prototype's scope."""


class _DeadlineExceededError(RuntimeError):
    pass


class PhiTransform:
    """Variable transform y (phi space) <-> x (production unknown space).

    Built from a resolved NetworkSolver whose _build_x0() has populated
    unknown_names / _name_to_index, plus the Jacobian at x0 (used to
    locate each junction's sum-mass row numerically: it is the unique
    row whose entries are exactly the +-1 port signs on the leg m_dot
    columns and nothing else -- independent of solver assembly order).
    """

    def __init__(self, solver: NetworkSolver, x0: np.ndarray):
        self.solver = solver
        names = solver.unknown_names
        n2i = solver._name_to_index
        self.n_x = len(names)

        # --- Collect junctions and their leg layout -------------------
        # derived_by: x index of a derived leg m_dot -> (jct_id, k)
        self.junctions: dict[str, dict[str, Any]] = {}
        derived_x: dict[int, tuple[str, int]] = {}
        for eid, elem in solver.network.elements.items():
            if not isinstance(elem, MultiPortChamberElement):
                continue
            if len(elem.inlet_nodes) == 1:
                common_idx = 0
            elif len(elem.outlet_nodes) == 1:
                common_idx = len(elem.inlet_nodes)
            else:
                raise PhiTransformError(
                    f"junction '{eid}': >1 inlets and >1 outlets is outside "
                    "the prototype scope (no unique common port)."
                )
            leg_eids = list(elem._port_element_ids)
            if any(not le for le in leg_eids):
                raise PhiTransformError(f"junction '{eid}': unresolved port elements.")
            common_eid = leg_eids[common_idx]
            noncommon = [le for i, le in enumerate(leg_eids) if i != common_idx]
            leg_x = [n2i[f"{le}.m_dot"] for le in leg_eids]
            for k, le in enumerate(noncommon):
                xi = n2i[f"{le}.m_dot"]
                if xi in derived_x:
                    raise PhiTransformError(
                        f"element '{le}' is a non-common leg of two junctions; "
                        "conflict is outside the prototype scope."
                    )
                derived_x[xi] = (eid, k)
            self.junctions[eid] = {
                "common_eid": common_eid,
                "common_x": n2i[f"{common_eid}.m_dot"],
                "noncommon_eids": noncommon,
                "leg_x": leg_x,
                "signs": [float(s) for s in elem._port_signs],
                "K": len(noncommon),
            }
        if not self.junctions:
            raise PhiTransformError("network has no MultiPortChamberElement.")
        self.derived_x = derived_x

        # --- y layout: kept x entries, then phi slots ------------------
        self.kept_x: list[int] = [i for i in range(self.n_x) if i not in derived_x]
        self.x_to_y = {xi: yi for yi, xi in enumerate(self.kept_x)}
        self.y_names: list[str] = [names[xi] for xi in self.kept_x]
        self.phi_slot: dict[tuple[str, int], int] = {}
        for jid in sorted(self.junctions):
            K = self.junctions[jid]["K"]
            for k in range(K - 1):
                self.phi_slot[(jid, k)] = len(self.y_names)
                self.y_names.append(f"{jid}.phi[{k}]")
        self.n_y = len(self.y_names)

        # --- Locate each junction's sum-mass row from J(x0) ------------
        res0, J0 = solver._residuals_and_jacobian(x0, compute_jacobian=True)
        J0 = J0.toarray()
        self.n_rows = J0.shape[0]
        drop_rows: set[int] = set()
        for jid, info in self.junctions.items():
            cols = info["leg_x"]
            signs = info["signs"]
            matches = []
            for r in range(self.n_rows):
                row = J0[r]
                nz = np.nonzero(np.abs(row) > 1e-12)[0]
                if set(nz.tolist()) != set(cols):
                    continue
                # Values must be exactly the port signs (duplicated legs
                # would sum; none of the benchmark nets have those).
                if all(abs(row[c] - s) < 1e-9 for c, s in zip(cols, signs, strict=True)):
                    matches.append(r)
            if len(matches) != 1:
                raise PhiTransformError(
                    f"junction '{jid}': expected exactly 1 sum-mass row, "
                    f"found {len(matches)} (rows {matches})."
                )
            drop_rows.add(matches[0])
        self.keep_rows = np.array([r for r in range(self.n_rows) if r not in drop_rows], dtype=int)

    # -- forward map -----------------------------------------------------
    def _phi(self, y: np.ndarray, jid: str, k: int) -> tuple[float, dict[int, float]]:
        """phi_k value and gradient dict {y_idx: d phi_k / d y}."""
        K = self.junctions[jid]["K"]
        if k < K - 1:
            slot = self.phi_slot[(jid, k)]
            return float(y[slot]), {slot: 1.0}
        val = 1.0
        grad: dict[int, float] = {}
        for j in range(K - 1):
            slot = self.phi_slot[(jid, j)]
            val -= float(y[slot])
            grad[slot] = -1.0
        return val, grad

    def _m_elem(
        self,
        y: np.ndarray,
        xi: int,
        cache: dict[int, tuple[float, dict[int, float]]],
    ) -> tuple[float, dict[int, float]]:
        """Element m_dot value + gradient (recursive through chaining)."""
        if xi in cache:
            return cache[xi]
        if xi not in self.derived_x:
            yi = self.x_to_y[xi]
            out = (float(y[yi]), {yi: 1.0})
        else:
            jid, k = self.derived_x[xi]
            phi, dphi = self._phi(y, jid, k)
            m_c, dm_c = self._m_elem(y, self.junctions[jid]["common_x"], cache)
            grad = {yi: phi * g for yi, g in dm_c.items()}
            for yi, g in dphi.items():
                grad[yi] = grad.get(yi, 0.0) + m_c * g
            out = (phi * m_c, grad)
        cache[xi] = out
        return out

    def x_from_y(self, y: np.ndarray) -> np.ndarray:
        x = np.empty(self.n_x)
        cache: dict[int, tuple[float, dict[int, float]]] = {}
        for xi in range(self.n_x):
            if xi in self.derived_x:
                x[xi] = self._m_elem(y, xi, cache)[0]
            else:
                x[xi] = y[self.x_to_y[xi]]
        return x

    def dx_dy(self, y: np.ndarray) -> np.ndarray:
        D = np.zeros((self.n_x, self.n_y))
        cache: dict[int, tuple[float, dict[int, float]]] = {}
        for xi in range(self.n_x):
            if xi in self.derived_x:
                _, grad = self._m_elem(y, xi, cache)
                for yi, g in grad.items():
                    D[xi, yi] = g
            else:
                D[xi, self.x_to_y[xi]] = 1.0
        return D

    def y_from_x(self, x: np.ndarray) -> np.ndarray:
        """Initial y: kept entries copied; phi from leg/common flow ratios."""
        y = np.empty(self.n_y)
        for xi in self.kept_x:
            y[self.x_to_y[xi]] = x[xi]
        for jid in sorted(self.junctions):
            info = self.junctions[jid]
            K = info["K"]
            m_c = abs(float(x[info["common_x"]]))
            n2i = self.solver._name_to_index
            for k in range(K - 1):
                if m_c > 1e-12:
                    xi = n2i[f"{info['noncommon_eids'][k]}.m_dot"]
                    phi0 = float(np.clip(abs(float(x[xi])) / m_c, 0.05, 0.95))
                else:
                    phi0 = 1.0 / K
                y[self.phi_slot[(jid, k)]] = phi0
            # Keep the implied last fraction strictly positive.
            slots = [self.phi_slot[(jid, k)] for k in range(K - 1)]
            if slots:
                total = sum(y[s] for s in slots)
                if total > 0.95:
                    for s in slots:
                        y[s] *= 0.95 / total
        return y

    def bounds(self, bound_all_mdots: bool = False) -> tuple[np.ndarray, np.ndarray]:
        lo = np.full(self.n_y, -np.inf)
        hi = np.full(self.n_y, np.inf)
        common_x = {info["common_x"] for info in self.junctions.values()}
        for yi, name in enumerate(self.y_names):
            if ".phi[" in name:
                lo[yi], hi[yi] = 0.0, 1.0
            elif name.endswith((".P", ".Pt", ".P_jct")):
                lo[yi] = 1.0  # 1 Pa floor avoids zero-density singularity
            elif name.endswith(".m_dot"):
                xi = self.kept_x[yi]
                if xi in common_x or bound_all_mdots:
                    lo[yi] = 0.0
        return lo, hi


def solve_phi(
    net: FlowNetwork,
    timeout: float = 60.0,
    max_nfev: int = 2000,
    bound_all_mdots: bool = False,
) -> dict[str, Any]:
    """Cold phi-closure solve. Returns a solver-style result dict."""
    t0 = time.perf_counter()
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = np.asarray(solver._build_x0(), dtype=float)
    tf = PhiTransform(solver, x0)
    lo, hi = tf.bounds(bound_all_mdots=bound_all_mdots)
    y0 = np.clip(tf.y_from_x(x0), lo + 1e-9, hi - 1e-9)

    # Row weighting derived from the Jacobian at y0: each row is
    # normalized by its characteristic magnitude ||J_red[r, :] * D_y||_inf
    # (residual change over a characteristic variable change). Without
    # this the raw least-squares landscape mixes Pa- and kg/s-magnitude
    # rows (cond ~ 8e8 observed) and trf crawls at maxF ~ 5e-2 -- the
    # UNREDUCED bounded trf shows the identical crawl, so this is a
    # weighting problem, not a transform problem. NOTE: the production
    # _build_residual_scales vector is NOT suitable here -- it scales
    # MPCE impulse rows (Pa magnitude) by ref_mdot, which a least-
    # squares objective cannot tolerate (breaks the three-port net
    # outright), while hybr root-finding is insensitive to it.
    D_y = np.abs(y0).copy()
    D_y[D_y < 1e-12] = 1.0

    deadline = t0 + timeout
    n_evals = [0]

    def raw_jac_red(y: np.ndarray) -> np.ndarray:
        _, J = solver._residuals_and_jacobian(tf.x_from_y(y), compute_jacobian=True)
        return (J.toarray() @ tf.dx_dy(y))[tf.keep_rows]

    J0 = raw_jac_red(y0)
    row_scale = np.max(np.abs(J0) * D_y[None, :], axis=1)
    row_scale[row_scale < 1e-12] = 1.0
    weights = 1.0 / row_scale

    def fun(y: np.ndarray) -> np.ndarray:
        if time.perf_counter() > deadline:
            raise _DeadlineExceededError
        n_evals[0] += 1
        try:
            x = tf.x_from_y(y)
            res, _ = solver._residuals_and_jacobian(x, compute_jacobian=False)
            return np.asarray(res, dtype=float)[tf.keep_rows] * weights
        except _DeadlineExceededError:
            raise
        except Exception:
            # Unphysical intermediate state: large flat penalty.
            return np.full(len(tf.keep_rows), 1.0e6)

    def jac(y: np.ndarray) -> np.ndarray:
        if time.perf_counter() > deadline:
            raise _DeadlineExceededError
        try:
            return weights[:, None] * raw_jac_red(y)
        except _DeadlineExceededError:
            raise
        except Exception:
            return np.eye(len(tf.keep_rows), tf.n_y)

    message = ""
    try:
        result = least_squares(
            fun,
            y0,
            jac=jac,
            bounds=(lo, hi),
            method="trf",
            x_scale="jac",
            max_nfev=max_nfev,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-10,
        )
        y_sol = result.x
        message = str(result.message)
    except _DeadlineExceededError:
        return {
            "__success__": False,
            "__message__": f"phi solve timed out after {timeout:.1f} s.",
            "__final_norm__": float("inf"),
            "__wall__": time.perf_counter() - t0,
            "__nfev__": n_evals[0],
            "__phi__": {},
            "__phi_pinned__": {},
        }

    x_sol = tf.x_from_y(y_sol)
    # Correctness gate on the FULL system (including the dropped mass
    # rows): the transform must satisfy the original problem, not a
    # different one.
    try:
        full_res = np.asarray(solver._residuals(x_sol), dtype=float)
        # Same gate as production _solve_impl: 2-norm of the raw
        # residual (including the dropped mass rows) below 1e-3.
        final_norm = float(np.linalg.norm(full_res))
    except Exception:
        final_norm = float("inf")
    success = final_norm < 1e-3

    sol = dict(zip(solver.unknown_names, x_sol.tolist(), strict=True))
    if success:
        bad = [
            eid
            for eid, elem in net.elements.items()
            if hasattr(elem, "verify_solution_consistent")
            and not elem.verify_solution_consistent(sol)
        ]
        if bad:
            success = False
            message = f"wrong-direction artifact root at junction(s) {bad}"

    phi_vals = {name: float(y_sol[i]) for i, name in enumerate(tf.y_names) if ".phi[" in name}
    pinned = {name: v for name, v in phi_vals.items() if v < 1e-6 or v > 1.0 - 1e-6}
    sol.update(
        {
            "__success__": success,
            "__message__": message,
            "__final_norm__": final_norm,
            "__wall__": time.perf_counter() - t0,
            "__nfev__": n_evals[0],
            "__phi__": phi_vals,
            "__phi_pinned__": pinned,
        }
    )
    return sol


# ---------------------------------------------------------------------------
# Smoke tests (run explicitly: uv run pytest python/tests/experimental_phi_closure.py)
# ---------------------------------------------------------------------------


def test_phi_three_port_reaches_mass_conserving_root():
    """The doomed-primary exemplar: production cold solve stalls in both
    hybr and LM phases; phi-closure should reach the root directly."""
    from test_multi_port_chamber import _three_port_net

    net = _three_port_net()
    sol = solve_phi(net, timeout=60.0)
    print(
        f"\n[three_port] success={sol['__success__']} |F|={sol['__final_norm__']:.3e} "
        f"wall={sol['__wall__']:.1f}s nfev={sol['__nfev__']} phi={sol['__phi__']}"
    )
    assert sol["__success__"], sol["__message__"]
    m_in = sol["ch_in.m_dot"]
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    assert abs(-m_in + m_str + m_bra) < 1e-6
    assert m_in > 0.0


def test_phi_matches_cold_root_on_branch_tee_orifice_net():
    """Same root as the production cold start on a healthy branch net."""
    from test_network_compressible import _mfb_branch_tee_orifice_network

    net_cold = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net_cold)
    sol_cold = solver.solve(timeout=60.0, auto_retry=False)
    assert sol_cold["__success__"], sol_cold.get("__message__")

    net_phi = _mfb_branch_tee_orifice_network(0.2, "compressible")
    sol_phi = solve_phi(net_phi, timeout=60.0)
    print(
        f"\n[mfb_tee] success={sol_phi['__success__']} |F|={sol_phi['__final_norm__']:.3e} "
        f"wall={sol_phi['__wall__']:.1f}s nfev={sol_phi['__nfev__']} phi={sol_phi['__phi__']}"
    )
    assert sol_phi["__success__"], sol_phi["__message__"]
    for key in ("ch_str.m_dot", "ch_bra.m_dot"):
        assert abs(sol_phi[key] - sol_cold[key]) / max(abs(sol_cold[key]), 1e-6) < 1e-3


def test_phi_reduced_jacobian_matches_fd():
    """Analytic J_red = J_x @ dx/dy against central differences."""
    from test_multi_port_chamber import _three_port_net

    net = _three_port_net()
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = np.asarray(solver._build_x0(), dtype=float)
    tf = PhiTransform(solver, x0)
    y0 = tf.y_from_x(x0)

    def f_red(y: np.ndarray) -> np.ndarray:
        res, _ = solver._residuals_and_jacobian(tf.x_from_y(y), compute_jacobian=False)
        return np.asarray(res, dtype=float)[tf.keep_rows]

    _, J = solver._residuals_and_jacobian(tf.x_from_y(y0), compute_jacobian=True)
    J_red = (J.toarray() @ tf.dx_dy(y0))[tf.keep_rows]

    J_fd = np.zeros_like(J_red)
    for j in range(tf.n_y):
        eps = max(abs(y0[j]) * 1e-7, 1e-9)
        yp, ym = y0.copy(), y0.copy()
        yp[j] += eps
        ym[j] -= eps
        J_fd[:, j] = (f_red(yp) - f_red(ym)) / (2.0 * eps)

    scale = np.maximum(np.abs(J_fd), np.abs(J_red)) + 1.0
    max_rel = float(np.max(np.abs(J_red - J_fd) / scale))
    print(f"\n[jac_fd] max rel err = {max_rel:.3e}")
    assert max_rel < 5e-3
