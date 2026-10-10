"""The C++ junction kernels against the Python, live and through their snapshots.

Production evaluates the Mynard closure in C++ (include/mynard_junction.h,
behind ``_core.mpce_residuals_and_jacobian``). Fidelity to Mynard's figures is
proven on the Python kernel (test_junction_mynard_fidelity.py), so it carries
over to production only while the two give the same numbers.

The C++ ctests hold the C++ to golden headers generated from the Python. A
snapshot cannot notice the Python moving, though, and one did: the
whole-element header predated the compressible reference head (#361). It
stored no gamma, so the C++ test ran the incompressible fallback and checked a
path production never takes, while the Python had drifted 1e-4 from it. Two
guards close that:

1. **Live.** The C++ K, read back through the binding, equals the Python K at
   every Mynard sample geometry, for each configuration knob the C++ takes.
   Measured 2026-10-10: 632 cases, max difference 1.1e-14.
2. **Fresh.** Both golden headers regenerate, in memory, to what is committed.

The recovery switch (``dividing_streamline_recovery``) is Python-only: the C++
carries production's recovery as a constant, so the live check compares the
Python at that default. Mynard's own configuration (recovery off) is a
Python-only result, and the term is additive and confined to collinear
continuing collectors.
"""

from __future__ import annotations

import importlib.util
import math
import re
from pathlib import Path

import numpy as np
import pytest

from combaero import _core
from combaero.network._mynard2010 import junction_loss_coefficient
from validation.junction import mynard_fidelity as mf

_DATA = Path(__file__).resolve().parents[2] / "validation" / "junction" / "data"
_RHO, _P = 1.2, 1.0e5


def _geometries() -> list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    spec_of = {s.key: s for s in mf.PANELS}
    out = []
    for p in mf.samples():
        spec = spec_of[p.panel]
        kw: dict = {"lam": 0.5, "phi_deg": None, "area": 1.0}
        kw[{"lam": "lam", "phi": "phi_deg", "area": "area"}.get(spec.sweep, "lam")] = (
            p.x if spec.sweep != "Re" else 0.5
        )
        Q, A, th, _ = mf._junction(spec.flow_type, **kw)
        out.append((p.panel, Q / A, A, th))
    return out


_GEOMETRIES = _geometries()


def _common(U: np.ndarray, A: np.ndarray) -> int:
    sup = U * A > 0.0
    return int(np.flatnonzero(sup)[0]) if sup.sum() == 1 else int(np.flatnonzero(~sup)[0])


def _cpp_K(
    U: np.ndarray, A: np.ndarray, theta: np.ndarray, alpha: float, eta: float
) -> tuple[np.ndarray, int]:
    geom = _core.MpceGeometry()
    geom.area = [float(a) for a in A]
    geom.theta_rad = [float(t) for t in theta]
    geom.port_sign = [-1.0, -1.0, -1.0]  # outer flow = junction flow, positive in
    geom.joining_etransfer_alpha = alpha
    geom.eta_scale = eta
    geom.gamma = [0.0, 0.0, 0.0]  # K does not depend on the reference head
    r = _core.mpce_residuals_and_jacobian(
        [_P] * 3, [_P] * 3, [_RHO] * 3, [_RHO / _P] * 3, [float(v) for v in U * _RHO * A], _P, geom
    )
    assert r.valid
    return np.asarray(r.k_per_port), int(r.common_port)


def _to_common_at_pi(U: np.ndarray, A: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """The element re-points its common port along the main duct at pi, so a
    geometry is handed over already rotated that way. Mynard's closure works
    on relative angles; the rotation is checked to change nothing."""
    return theta + (math.pi - theta[_common(U, A)])


def test_rotation_leaves_the_python_kernel_unchanged():
    worst = 0.0
    for _, U, A, th in _GEOMETRIES:
        a = np.atleast_1d(junction_loss_coefficient(U, A, th, 0.0, 1.0).K)
        b = np.atleast_1d(junction_loss_coefficient(U, A, _to_common_at_pi(U, A, th), 0.0, 1.0).K)
        worst = max(worst, float(np.abs(a - b).max()))
    assert worst < 1e-13


@pytest.mark.parametrize("alpha", [0.0, 0.2])
@pytest.mark.parametrize("eta", [0.0, 1.0])
def test_cpp_kernel_equals_python_at_every_mynard_sample(alpha, eta):
    for panel, U, A, th in _GEOMETRIES:
        rot = _to_common_at_pi(U, A, th)
        k_cpp, common = _cpp_K(U, A, rot, alpha, eta)
        assert common == _common(U, A), panel
        k_py = np.atleast_1d(junction_loss_coefficient(U, A, rot, alpha, eta).K)
        others = [i for i in range(3) if i != common]
        np.testing.assert_allclose(k_cpp[others], k_py, rtol=0.0, atol=1e-12, err_msg=panel)


def test_the_live_comparison_sees_a_difference():
    """Falsify the check: the C++ with eta on against the Python with it off
    must disagree somewhere, by far more than the tolerance."""
    worst = 0.0
    for _, U, A, th in _GEOMETRIES:
        rot = _to_common_at_pi(U, A, th)
        k_cpp, common = _cpp_K(U, A, rot, 0.0, 1.0)
        k_py = np.atleast_1d(junction_loss_coefficient(U, A, rot, 0.0, 0.0).K)
        others = [i for i in range(3) if i != common]
        worst = max(worst, float(np.abs(k_cpp[others] - k_py).max()))
    assert worst > 0.1


# ---------------------------------------------------------------------------
# Snapshot freshness
# ---------------------------------------------------------------------------

_NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:e[-+]?\d+)?|nan|inf", re.IGNORECASE)


def _render(name: str) -> str:
    spec = importlib.util.spec_from_file_location(name, _DATA / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.render()


def _split(text: str) -> tuple[str, np.ndarray]:
    numbers = np.array([float(x) for x in _NUMBER.findall(text)])
    return _NUMBER.sub("#", text), numbers


@pytest.mark.parametrize(
    ("generator", "header"),
    [
        ("generate_mynard_reference", "mynard_reference_data.h"),
        ("generate_mpce_reference", "mpce_reference_data.h"),
    ],
)
def test_golden_header_is_current(generator, header):
    """The committed header is what the Python produces today. Compared
    number by number, not byte by byte: the Jacobians are central differences,
    and in the mass-flow columns (step ~1e-6 kg/s) a last-bit difference in
    the residual becomes ~1e-6 absolute, which can differ between platforms.
    The tolerance (1e-5 absolute, 1e-7 relative) still fails the drift this
    test was written after by four orders: 0.3 Pa on a 3085 Pa residual."""
    fresh_text, fresh = _split(_render(generator))
    committed_text, committed = _split((_DATA / header).read_text())
    assert fresh_text == committed_text, (
        f"{header}: structure changed -- regenerate it with {generator}.py"
    )
    np.testing.assert_allclose(
        fresh,
        committed,
        rtol=1e-7,
        atol=1e-5,
        err_msg=f"{header} is stale -- regenerate it with {generator}.py",
    )
