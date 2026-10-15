"""TEMPORARY Linux diagnostic for PR #579's Validation failure. Not for merge."""

import warnings

import numpy as np
from scipy.optimize._numdiff import approx_derivative
from test_constant_k_tee import _merge_net

import combaero as cb
from combaero import _solver_tools as st
from combaero.network import NetworkSolver


def test_diag() -> None:
    lines = []
    net = _merge_net(0.5, 1.2)
    s = NetworkSolver(net)
    net.resolve_all_topology()
    net.validate()
    x0 = np.asarray(s._build_x0(), dtype=float)
    r0, J = s._residuals_and_jacobian(x0)
    J = J.toarray()
    fd = approx_derivative(
        lambda x: s._residuals(x),
        x0,
        method="2-point",
        abs_step=np.maximum(np.abs(x0) * 1e-7, 1e-10),
    )
    sc = np.maximum(np.abs(fd), np.abs(J)) + 1.0
    e = np.abs(J - fd) / sc
    for k in np.argsort(-e.ravel())[:6]:
        i, j = np.unravel_index(k, e.shape)
        lines.append(
            f"row {i} col {s.unknown_names[j]} an {J[i, j]:.9e} fd {fd[i, j]:.9e} err {e[i, j]:.2e}"
        )
    lines.append("r0 " + " ".join(f"{v:.6e}" for v in r0))
    X = list(cb.species.dry_air())
    Y = list(cb.mole_to_mass(X))
    for T0 in (300.0, 318.0, 320.0):
        for m in (0.12, 0.18, 0.3):
            a = st.momentum_chamber_residual_and_jacobian(1.9e5, 1.9e5, m, T0, Y, 1.9635e-3)
            b = st.momentum_chamber_residual_and_jacobian(1.9e5, 1.9e5, m, T0 + 1e-4, Y, 1.9635e-3)
            c = st.momentum_chamber_residual_and_jacobian(1.9e5, 1.9e5, m, T0 - 1e-4, Y, 1.9635e-3)
            lines.append(
                f"T0 {T0} m {m}: Ts {a.T_static:.12f} M {a.M:.6f} res {a.residual:.9e} "
                f"dT an {a.d_res_dT:.9e} fd {(b.residual - c.residual) / 2e-4:.9e} "
                f"Ts+ {b.T_static:.12f} Ts- {c.T_static:.12f}"
            )
    warnings.warn("CI-DIAG-579\n" + "\n".join(lines), UserWarning, stacklevel=1)
