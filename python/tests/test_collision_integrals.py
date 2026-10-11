"""Collision integral Omega*(2,2): exact Lennard-Jones nodes, tabulated slopes,
C1 cubic Hermite interpolation, analytic transport derivatives (#485).

The delta* = 0 column comes from quadrature of the classical collision
integrals (thermo_data_generator/collision_integrals.py). It is checked here
against Kim & Monroe (2014, J. Comput. Phys. 273:358), an independent
high-accuracy computation, at the nodes and between them. The former
Monchick-Mason column was off by up to 0.6% at the nodes, and linear
interpolation in raw T* added up to 2.7% between them with 2-17% kinks in
dmu/dT.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import combaero as cb
from combaero import _core

# Kim & Monroe (2014) Omega*(2,2) fit, 0.3 <= T* <= 400, with the scaling of
# their coefficient table (as implemented in the chapensk package).
_KM_A = -0.92032979
_KM_B = [2.3508044, 0.50110649, -0.47193769, 0.15806367, -0.026367184, 0.0018120118]
_KM_C = [1.6330213, -0.69795156, 0.16096572, -0.022109440, 0.0017031434, -0.000056699986]

NODES = [
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    1.2,
    1.4,
    1.6,
    1.8,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    10.0,
    12.0,
    14.0,
    16.0,
    18.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
    50.0,
    75.0,
    100.0,
    150.0,
    200.0,
    300.0,
    400.0,
]


def _km(T: float) -> float:
    return _KM_A + sum(
        _KM_B[k] / T ** (k + 1) + _KM_C[k] * math.log(T) ** (k + 1) for k in range(6)
    )


def _mix(d: dict[str, float]) -> list[float]:
    X = [0.0] * cb.num_species()
    for k, v in d.items():
        X[cb.species_index_from_name(k)] = v
    s = sum(X)
    return [x / s for x in X]


@pytest.mark.parametrize("where", ["nodes", "between"])
def test_lennard_jones_column_matches_kim_monroe(where: str) -> None:
    """<= 1e-4 everywhere in 0.3-400, between nodes too (the interpolant)."""
    if where == "nodes":
        pts = NODES
    else:
        pts = [
            math.exp((1 - f) * math.log(a) + f * math.log(b))
            for a, b in zip(NODES, NODES[1:], strict=False)
            for f in (0.25, 0.5, 0.75)
        ]
    dev = max(abs(_core.omega22_and_derivative(t, 0.0).omega / _km(t) - 1) for t in pts)
    assert dev < 1e-4


@pytest.mark.parametrize("delta", [0.0, 0.25, 0.6, 1.0, 1.7, 2.5])
def test_the_interpolant_is_c1_at_every_node(delta: float) -> None:
    """The slope is the same on both sides of a node (it jumped 2-17%)."""
    for t in [0.2] + NODES[:-1]:
        lo = _core.omega22_and_derivative(t * (1 - 1e-9), delta).domega_dTstar
        hi = _core.omega22_and_derivative(t * (1 + 1e-9), delta).domega_dTstar
        assert hi == pytest.approx(lo, rel=1e-6)


@pytest.mark.parametrize("delta", [0.0, 0.25, 1.0, 2.5])
def test_the_collision_integral_derivative_is_its_slope(delta: float) -> None:
    for t in (0.15, 0.45, 1.3, 4.7, 23.0, 87.0, 260.0, 600.0):
        r = _core.omega22_and_derivative(t, delta)
        h = 1e-6 * t
        fd = (
            _core.omega22_and_derivative(t + h, delta).omega
            - _core.omega22_and_derivative(t - h, delta).omega
        ) / (2 * h)
        assert r.domega_dTstar == pytest.approx(fd, rel=1e-6)


def test_polar_columns_keep_the_tables_ratio_to_the_lennard_jones_column() -> None:
    """delta* > 0 is Monchick-Mason scaled by the delta* = 0 column's
    exact/tabulated ratio, so each polar column keeps the table's ratio to its
    own LJ column. (Spot checks against the table: T* = 10 and 100.)"""
    mm = {  # (T*, delta*) -> Monchick-Mason value; delta* = 0 first
        10.0: {0.0: 0.82435, 1.0: 0.8356, 2.5: 0.8901},
        100.0: {0.0: 0.5887, 1.0: 0.5903, 2.5: 0.5885},
    }
    for t, row in mm.items():
        o0 = _core.omega22_and_derivative(t, 0.0).omega
        for d in (1.0, 2.5):
            assert _core.omega22_and_derivative(t, d).omega / o0 == pytest.approx(
                row[d] / row[0.0], rel=1e-9
            )


CASES = {
    "air": None,
    "H2O/N2": {"H2O": 0.3, "N2": 0.7},
    "NH3/air": {"NH3": 0.2, "N2": 0.62, "O2": 0.18},
    "H2": {"H2": 1.0},
    "CO2": {"CO2": 1.0},
    "products": {"N2": 0.71, "H2O": 0.18, "CO2": 0.09, "O2": 0.02},
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_transport_derivatives_are_analytic_and_exact(case: str) -> None:
    """d mu/dT and d k/dT through the interpolation, Eucken factors and both
    mixing rules. Away from 1000 K: the NASA polynomials' range join puts a
    kink in cp, so k's derivative is one-sided exactly there."""
    X = list(cb.species.dry_air()) if CASES[case] is None else _mix(CASES[case])
    for T in np.concatenate([np.linspace(230.0, 990.0, 20), np.linspace(1010.0, 3000.0, 20)]):
        r = cb.transport_and_dT(T, 1e5, X)
        h = 1e-6 * T
        fmu = (cb.viscosity(T + h, 1e5, X) - cb.viscosity(T - h, 1e5, X)) / (2 * h)
        fk = (cb.thermal_conductivity(T + h, 1e5, X) - cb.thermal_conductivity(T - h, 1e5, X)) / (
            2 * h
        )
        assert r.dmu_dT == pytest.approx(fmu, rel=1e-6)
        assert r.dk_dT == pytest.approx(fk, rel=1e-6)
        assert r.mu == pytest.approx(cb.viscosity(T, 1e5, X), rel=1e-14)


def test_air_viscosity_is_smooth_in_temperature() -> None:
    """Second differences of mu(T): spike ratio ~1700 with the old kinks."""
    T = np.linspace(250.0, 300.0, 2001)
    mu = np.array([cb.viscosity(t, 1.5e5, cb.species.dry_air()) for t in T])
    d2 = np.abs(np.diff(mu, 2))
    assert d2.max() / np.median(d2) < 2.0
