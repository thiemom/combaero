"""The LM fallback on a Jacobian that is singular but for underflow-level entries (#1020).

MINPACK lmder treats a Jacobian as singular only on an exact zero pivot. On
stress effusion 34 (an EffusionLiner, 5 segments) a far iterate drives every
plate's McGreehan-Schotsch crossflow Cd = exp(-Rv^n) to ~1e-308, and that
plate row is the only one fixing chamber ln__b1's pressure level: columns
ln__b1.P and ln__b1.Pt cancel exactly elsewhere. The pivot came out 1e-308,
the step overflowed to NaN, lmder's trust radius became NaN, and every later
LM iterate was NaN (129 of 196 evaluations in an 8 s run; never converged).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import combaero.network.solver as solver_mod
from combaero.network import NetworkSolver

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402


def test_lm_fallback_never_proposes_a_non_finite_iterate(monkeypatch: pytest.MonkeyPatch) -> None:
    """Effusion 34 cold, as the stress sweep solves it minus the retries: the
    LM phase sees only finite iterates and converges (before: NaN from LM
    evaluation 16 on, outcome no_progress at |F| 8.78)."""
    lm_x: list[np.ndarray] = []
    real_root = solver_mod.root

    def traced_root(fun: Any, x0: np.ndarray, method: str = "hybr", **kw: Any) -> Any:
        if method != "lm":
            return real_root(fun, x0, method=method, **kw)

        def f(x: np.ndarray) -> Any:
            lm_x.append(np.array(x, copy=True))
            return fun(x)

        return real_root(f, x0, method=method, **kw)

    monkeypatch.setattr(solver_mod, "root", traced_root)
    net, _ = FAMILIES["effusion"](34)
    solver = NetworkSolver(net)
    solver.homotopy_retry = False
    solver.bc_ramp_retry = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = solver.solve(timeout=60.0)

    assert lm_x, "the LM fallback did not run: the case no longer exercises it"
    n_bad = sum(not np.all(np.isfinite(x)) for x in lm_x)
    assert n_bad == 0, f"{n_bad} of {len(lm_x)} LM iterates non-finite"
    assert res["__converged__"], res.get("__message__")
    assert res["__final_norm__"] < 1e-3
