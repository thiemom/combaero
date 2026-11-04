"""The orifice row's Jacobian against central differences, every column (#964).

Three gaps, all pre-existing:

* The Reynolds-dependent Cd (Reader-Harris/Gallagher, the Idelchik
  and McGreehan-Schotsch holes, Lichtarowicz) reads Re ~ |m| / mu(T, Y). The
  row carried mu's temperature column but not its composition: with an RHG
  orifice behind a combustor the row's fuel-flow column was 7e-4 off.
* The kernels' composition (Y) columns were forward differences, up to 3e-3
  of the row off a central one. They are analytic now: the density for the
  incompressible kernel, the nozzle's isentrope (and the approach static
  temperature on it), Jobson's and ISO's gamma and R by DualN, the approach
  pipe's static pressure subsonic and past choke.
* Idelchik's beveled and rounded holes returned their tables' GEOMETRY slope
  as dCd/dRe (0.83 and 1.53 per unit Re for a Cd constant in Re), so the
  self-consistent row's denominator 1 - K dCd/dm was of order -1e5 and every
  column but the flow's was ~0.

The error is scaled per column, |J_ij - FD_ij| s_j / max_k |J_ik| s_k with s
the column's scale (the pressure, the temperature, 1 for a mass fraction, the
flow), so columns of different units compare.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    ChannelElement,
    CombustorNode,
    LosslessConnectionElement,
    MassFlowBoundary,
    NetworkMixtureState,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

NS = cb.num_species()
TOL = 1e-5


def _burnt_Y() -> list[float]:
    X = [0.0] * NS
    for name, x in (("N2", 0.72), ("O2", 0.10), ("CO2", 0.06), ("H2O", 0.11), ("AR", 0.009)):
        X[cb.species_index_from_name(name)] = x
    X[cb.species_index_from_name("CH4")] = 0.001
    s = sum(X)
    return list(cb.mole_to_mass([x / s for x in X]))


_CORR_KW = {
    "fixed": {"Cd": 0.62},
    "ReaderHarrisGallagher": {},
    "IdelchikThick": {"plate_thickness": 0.004},
    "IdelchikBeveled": {"bevel_depth": 0.001},
    "IdelchikRounded": {"edge_radius": 0.001},
    "McGreehanSchotsch": {"plate_thickness": 0.02},
}


def _orifice(corr: str, regime: str, model: str | None, beta: float) -> OrificeElement:
    """An orifice between two plena, in a pipe of diameter d / beta for
    beta > 0, else a hole in a wall."""
    d = 0.01
    n = FlowNetwork()
    n.add_node(PressureBoundary("src", Pt=10e5, Tt=900.0, Y=_burnt_Y()))
    n.add_node(PlenumNode("a"))
    n.add_node(PlenumNode("b"))
    n.add_node(PressureBoundary("out", Pt=5e5, Tt=900.0, Y=_burnt_Y()))
    in_pipe = beta > 0.0
    n.add_element(
        ChannelElement("pu", "src", "a", length=0.5, diameter=d / beta if in_pipe else 0.05)
    )
    ox = OrificeElement(
        "ox",
        "a",
        "b",
        diameter=d,
        correlation=corr,
        regime=regime,
        in_pipe=in_pipe,
        compressible_model=model,
        **_CORR_KW[corr],
    )
    n.add_element(ox)
    n.add_element(LosslessConnectionElement("ld", "b", "out"))
    n.resolve_all_topology()
    return ox


def _residual(ox: OrificeElement, x: dict[str, float]) -> tuple[float, dict[str, float]]:
    Y = [x[f"a.Y[{k}]"] for k in range(NS)]
    a = NetworkMixtureState(P=x["a.P"], Pt=x["a.Pt"], Tt=x["a.Tt"], m_dot=x["ox.m_dot"], Y=Y)
    b = NetworkMixtureState(P=x["b.P"], Pt=x["b.Pt"], Tt=900.0, m_dot=x["ox.m_dot"], Y=_burnt_Y())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        F, J = ox.residuals(a, b)
    return float(F[0]), dict(J[0])


def _state(ox: OrificeElement, PR: float) -> dict[str, float]:
    """The row at its own consistent flow, Pt = 10 bar, Tt = 900 K."""
    x = {"a.P": 10e5, "a.Pt": 10e5, "a.Tt": 900.0, "b.P": PR * 10e5, "b.Pt": PR * 10e5}
    x["ox.m_dot"] = 0.05
    for k, y in enumerate(_burnt_Y()):
        x[f"a.Y[{k}]"] = y
    for _ in range(6):
        x["ox.m_dot"] -= _residual(ox, x)[0]
    return x


def _scale(key: str, x: dict[str, float]) -> float:
    return 1.0 if ".Y[" in key else abs(x[key])


def _fd(ox: OrificeElement, x: dict[str, float], key: str) -> float:
    v = x[key]
    h = 1e-6 if ".Y[" in key else 1e-6 * abs(v)
    if ".Y[" in key and v < 2 * h:
        f = [_residual(ox, {**x, key: v + i * h})[0] for i in range(3)]
        return (-3 * f[0] + 4 * f[1] - f[2]) / (2 * h)
    return (_residual(ox, {**x, key: v + h})[0] - _residual(ox, {**x, key: v - h})[0]) / (2 * h)


def _worst(ox: OrificeElement, x: dict[str, float]) -> tuple[float, str]:
    _, J = _residual(ox, x)
    keys = list(x)
    fd = {k: _fd(ox, x, k) for k in keys}
    smax = max(max(abs(J.get(k, 0.0)), abs(fd[k])) * _scale(k, x) for k in keys)
    errs = {k: abs(J.get(k, 0.0) - fd[k]) * _scale(k, x) / smax for k in keys}
    k = max(errs, key=errs.get)
    return errs[k], k


CASES = [
    # (correlation, regime, model, beta, P_down / Pt)
    ("ReaderHarrisGallagher", "incompressible", None, 0.5, 0.7),
    ("ReaderHarrisGallagher", "compressible", "nozzle", 0.5, 0.7),
    ("ReaderHarrisGallagher", "compressible", "nozzle", 0.5, 0.4),
    ("ReaderHarrisGallagher", "compressible", "jobson1955", 0.5, 0.4),
    ("ReaderHarrisGallagher", "compressible", "iso5167", 0.5, 0.95),
    # A wide bore: the approach pipe beyond its sonic limit at the root.
    ("ReaderHarrisGallagher", "compressible", "nozzle", 0.99, 0.4),
    ("fixed", "compressible", "nozzle", 0.0, 0.7),
    ("fixed", "compressible", "nozzle", 0.95, 0.4),
    ("IdelchikThick", "compressible", "nozzle", 0.0, 0.7),
    ("IdelchikThick", "incompressible", None, 0.5, 0.7),
    ("McGreehanSchotsch", "compressible", "nozzle", 0.0, 0.7),
    ("IdelchikBeveled", "compressible", "nozzle", 0.0, 0.7),
    ("IdelchikRounded", "incompressible", None, 0.0, 0.7),
]


@pytest.mark.parametrize(("corr", "regime", "model", "beta", "PR"), CASES)
def test_every_orifice_column_matches_a_central_difference(
    corr: str, regime: str, model: str | None, beta: float, PR: float
) -> None:
    ox = _orifice(corr, regime, model, beta)
    x = _state(ox, PR)
    for fac in (1.0, 1.3):
        err, key = _worst(ox, {**x, "ox.m_dot": fac * x["ox.m_dot"]})
        assert err < TOL, (key, err, fac)


def test_the_wide_bore_case_is_past_its_approach_pipes_sonic_limit() -> None:
    """The beta 0.99 case above really reads the choked approach branch."""
    ox = _orifice("ReaderHarrisGallagher", "compressible", "nozzle", 0.99)
    x = _state(ox, 0.4)
    A = math.pi * ox.upstream_diameter**2 / 4
    Y = [x[f"a.Y[{k}]"] for k in range(NS)]
    st = _solver_tools.approach_static_pressure(10e5, 900.0, abs(x["ox.m_dot"]) / A, Y, False)
    assert st.choked


def test_beveled_and_rounded_holes_have_no_reynolds_slope() -> None:
    for corr, hole in (
        ("Idelchik1966Beveled", cb.DischargeHoleGeometry(d=0.01, L=0.0, r=0.0)),
        ("Idelchik1966Rounded", cb.DischargeHoleGeometry(d=0.01, L=0.0, r=0.001)),
    ):
        if corr == "Idelchik1966Beveled":
            hole.bevel = 0.001
        sel = getattr(cb.DischargeCdCorrelation, corr)
        _, dRe, _ = cb.discharge_cd_and_derivatives(sel, hole, cb.DischargeHoleState(Re=1e5))
        assert dRe == 0.0, corr


def _combustor_net(corr: str) -> FlowNetwork:
    """#788's network: air and fuel into a combustor, a channel choked into a
    junction, an in-pipe exit orifice. The fuel flow moves the orifice's
    upstream composition, hence its Cd's viscosity."""
    X_AIR = cb.species.dry_air()
    Y_AIR = list(cb.mole_to_mass(X_AIR))
    X = [0.0] * NS
    X[cb.species_index_from_name("CH4")] = 1.0
    Y_CH4 = list(cb.mole_to_mass(X))
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=32.0, Tt=897.0, Y=Y_AIR))
    n.add_node(MassFlowBoundary("fuel", m_dot=1.0, Tt=300.0, Y=Y_CH4))
    n.add_node(PlenumNode("pl"))
    n.add_node(CombustorNode("cmb"))
    n.add_node(PlenumNode("jct"))
    n.add_node(PressureBoundary("out", Pt=20e5, Tt=300.0, Y=Y_AIR))
    n.add_element(LosslessConnectionElement("ea", "air", "pl"))
    n.add_element(OrificeElement("of", "fuel", "pl", Cd=0.6, diameter=0.02, correlation="fixed"))
    n.add_element(LosslessConnectionElement("ep", "pl", "cmb"))
    n.add_element(ChannelElement("ch", "cmb", "jct", length=1.0, diameter=0.1))
    n.add_element(
        OrificeElement(
            "ox",
            "jct",
            "out",
            diameter=0.099,
            correlation=corr,
            in_pipe=True,
            compressible_model="nozzle",
        )
    )
    return n


def test_the_rhg_orifice_rows_fuel_column_behind_a_combustor() -> None:
    s = NetworkSolver(_combustor_net("ReaderHarrisGallagher"))
    s.homotopy_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30.0)
    assert r["__success__"]
    assert r["ch.mach_out"] == pytest.approx(1.0, abs=1e-9)
    x = np.array(r["__x_solution__"])
    F, J = s._residuals_and_jacobian(x)
    J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    i = next(k for k, row in enumerate(s._row_kinds(x)) if "'ox'" in str(row))
    for name in ("of.m_dot", "ea.m_dot"):
        j = s.unknown_names.index(name)
        # A wider step than the default: the row's own consistent-flow solve
        # (to 1e-13) reads as noise at 1e-6.
        h = 1e-4 * abs(x[j])
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        fd = (s._residuals_and_jacobian(xp)[0][i] - s._residuals_and_jacobian(xm)[0][i]) / (2 * h)
        assert abs(J[i, j] - fd) / np.max(np.abs(J[i])) < 1e-5, name
