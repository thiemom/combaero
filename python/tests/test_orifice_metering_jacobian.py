"""The metering orifice's own flow column carries dCd/dm (through Re_D).

ReaderHarrisGallagher (the default) and Stolz read the pipe Reynolds
number Re_D = 4 |m| / (pi D mu), so the row m - Cd(Re_D(m)) K has slope
1 - K dCd/dm. The Jacobian gave 1: 1.0% off on a 2 cm orifice at 1600 K,
found by the FD check behind #636's tests.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


@pytest.mark.parametrize("Re", [500.0, 5e3, 3e4, 1e5, 1e6])
def test_the_rhg_derivative_is_exact_and_its_value_unchanged(Re: float) -> None:
    g = cb.OrificeGeometry()
    g.d, g.D = 0.02, 0.05

    def st(r: float) -> cb.OrificeState:
        s = cb.OrificeState()
        s.Re_D = r
        return s

    cd, dcd = cb.Cd_sharp_thin_plate_and_dRe(g, st(Re))
    assert cd == cb.Cd_sharp_thin_plate(g, st(Re))
    h = 1e-6 * Re
    fd = (cb.Cd_sharp_thin_plate(g, st(Re + h)) - cb.Cd_sharp_thin_plate(g, st(Re - h))) / (2 * h)
    assert dcd == pytest.approx(fd, rel=1e-6)


def _net(regime: str, correlation: str) -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.6e5, Tt=1500.0))
    g.add_node(PlenumNode("hp"))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    for eid, a, b, d in (("o1", "hi", "hp", 0.02), ("o2", "hp", "ho", 0.015)):
        g.add_element(
            OrificeElement(eid, a, b, Cd=0.8, diameter=d, regime=regime, correlation=correlation)
        )
    return g


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("correlation", ["ReaderHarrisGallagher", "Stolz"])
def test_the_orifice_rows_flow_column_is_exact(regime: str, correlation: str) -> None:
    s = NetworkSolver(_net(regime, correlation))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    rows = s._row_kinds(x)
    for eid in ("o1", "o2"):
        i = next(k for k, (owner, _) in enumerate(rows) if f"'{eid}'" in owner)
        j = s.unknown_names.index(f"{eid}.m_dot")
        h = 1e-6 * abs(x[j])
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
        assert J[i, j] == pytest.approx(fd, rel=1e-7), (eid, correlation, regime)
        # The column really carries dCd/dm. Stolz's Re term (0.0029 b^2.5
        # (1e6/Re)^0.75) is the weaker one: |J - 1| ~ 7e-5 here; a missing
        # derivative gives exactly 0.
        assert abs(J[i, j] - 1.0) > 1e-5
