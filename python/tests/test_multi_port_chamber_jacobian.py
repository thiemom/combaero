"""FD guardrail for the MPCE-v1 analytic Jacobian.

``multi_port_chamber_residuals_and_jacobian`` returns residuals and every
partial analytically, but nothing pinned those partials against the residual
they claim to differentiate. The whole-row gap is where a wrong derivative
hides: it costs iterations and robustness rather than correctness of a
converged root, so a suite can stay green while a partial is simply wrong
(exactly the failure found in MPCEv2, issue #271).

Any change to the impulse residual -- notably the angle projection debated in
issue #272 -- must keep these within tolerance.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero import _solver_tools

_Y = list(cb.mole_to_mass(cb.species.dry_air()))

# Separating tee: port 0 common inflow (axial), 1 straight, 2 lateral.
_BASE = {
    "P_jct": 100_000.0,
    "P": [100_000.0, 99_800.0, 99_500.0],
    "mdot": [-0.10, 0.06, 0.04],
    "T": [300.0, 300.5, 301.0],
    "Y": [_Y, _Y, _Y],
    "A": [0.01, 0.01, 0.008],
    "theta_rad": [0.0, 0.0, math.pi / 2.0],
}

# eps chosen per variable so the central difference is well conditioned in the
# variable's own units: Pa, K, kg/s.
_PERTURBATION = {"P": 1.0, "T": 1e-3, "mdot": 1e-6}


def _residuals(**overrides: object) -> list[float]:
    kwargs = {**_BASE, **overrides}
    result = _solver_tools.multi_port_chamber_residuals_and_jacobian(**kwargs)
    return list(result.impulse_residuals)


def _central_difference(key: str, index: int, row: int) -> float:
    eps = _PERTURBATION[key]
    up, down = list(_BASE[key]), list(_BASE[key])
    up[index] += eps
    down[index] -= eps
    return (_residuals(**{key: up})[row] - _residuals(**{key: down})[row]) / (2.0 * eps)


@pytest.fixture(scope="module")
def analytic():
    return _solver_tools.multi_port_chamber_residuals_and_jacobian(**_BASE)


@pytest.mark.parametrize("port", [0, 1, 2])
@pytest.mark.parametrize("name, key", [("dR_dP", "P"), ("dR_dT", "T"), ("dR_dmdot", "mdot")])
def test_local_port_partials_match_finite_difference(analytic, port: int, name: str, key: str):
    """Each port's impulse row differentiated wrt its own state."""
    expected = _central_difference(key, port, port)
    actual = getattr(analytic.port_jac[port], name)
    assert actual == pytest.approx(expected, rel=1e-5, abs=1e-9), (
        f"port {port} {name}: analytic={actual:.6e} fd={expected:.6e}"
    )


@pytest.mark.parametrize("port", [1, 2])
@pytest.mark.parametrize(
    "name, key",
    [
        ("cross_dR_dP_axial", "P"),
        ("cross_dR_dT_axial", "T"),
        ("cross_dR_dmdot_axial", "mdot"),
    ],
)
def test_cross_coupling_partials_match_finite_difference(analytic, port: int, name: str, key: str):
    """Non-axial rows also depend on the axial reference port's state.

    This is the coupling that carries Bassett's -2*q*cos((3/4)*theta) term into
    K6; a wrong partial here degrades the lateral silently.
    """
    expected = _central_difference(key, 0, port)
    actual = getattr(analytic, name)[port]
    assert actual == pytest.approx(expected, rel=1e-5, abs=1e-9), (
        f"port {port} {name}: analytic={actual:.6e} fd={expected:.6e}"
    )


def test_mass_residual_is_the_port_sum(analytic):
    assert analytic.mass_residual == pytest.approx(sum(_BASE["mdot"]))
