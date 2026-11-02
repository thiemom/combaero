"""#932: a ThermalWall whose layers' k depends on T. The heat moves with each
layer's t/k, and t/k with the layer's mean temperature at the k(T) fixed
point; the Jacobian now carries that path, with dk/dT analytic from the
material data (the #882 envelope included)."""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
MATERIALS = [m for m in cb.list_materials() if m != "al6061"]


def _grid(name: str) -> list[float]:
    info = cb.get_material_info(name)
    lo, hi = info["T_envelope_min"], info["T_envelope_max"]
    # Below, across and above the data and the envelope, at odd offsets
    # (off the tables' printed points).
    return list(np.linspace(0.5 * lo, hi + 150.0, 157) + 0.3137)


@pytest.mark.parametrize("extrapolate", [False, True])
@pytest.mark.parametrize("name", MATERIALS)
def test_the_slope_is_the_conductivity_s(name: str, extrapolate: bool) -> None:
    """k is get_material_conductivity's bit for bit; dk/dT is its derivative
    wherever it is smooth (central differences agree from both sides)."""
    for T in _grid(name):
        try:
            k_ref = cb.get_material_conductivity(name, T, extrapolate=extrapolate)
        except RuntimeError:
            with pytest.raises(RuntimeError):
                cb.get_material_conductivity_and_slope(name, T, extrapolate=extrapolate)
            continue
        k, dk = cb.get_material_conductivity_and_slope(name, T, extrapolate=extrapolate)
        assert k == k_ref
        h = 1e-4

        def kk(t: float) -> float:
            return cb.get_material_conductivity(name, t, extrapolate=extrapolate)

        left = (kk(T) - kk(T - h)) / h
        right = (kk(T + h) - kk(T)) / h
        if not math.isclose(left, right, rel_tol=1e-4, abs_tol=1e-6):
            continue  # a printed point or the fade's end: only C0 there
        assert dk == pytest.approx(0.5 * (left + right), rel=1e-4, abs=1e-6), (name, T)


def test_the_envelope_slope_fades_to_zero() -> None:
    """At the data's end the slope is the end slope; it ramps to 0 over the
    100 K fade and is 0 beyond, and the linear continuation keeps it."""
    info = cb.get_material_info("hastelloy_x")
    T_max = info["T_max"]
    _, end = cb.get_material_conductivity_and_slope("hastelloy_x", T_max - 1e-3)
    _, half = cb.get_material_conductivity_and_slope("hastelloy_x", T_max + 50.0)
    _, flat = cb.get_material_conductivity_and_slope("hastelloy_x", T_max + 150.0)
    _, lin = cb.get_material_conductivity_and_slope("hastelloy_x", T_max + 150.0, extrapolate=True)
    assert end > 0.0
    assert half == pytest.approx(0.5 * end, rel=1e-6)
    assert flat == 0.0
    assert lin == pytest.approx(end, rel=1e-6)


def _fixed_point(p: list[float], t: list[float], kfun: list) -> list[float]:
    """t/k_i(Tm_i) iterated to convergence (the wall's own loop)."""
    T_hot, T_cold, h_hot, h_cold = p
    r = [ti / kf(0.5 * (T_hot + T_cold)) for ti, kf in zip(t, kfun, strict=True)]
    for _ in range(200):
        prof = cb.wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, r, 1e-4)[0]
        new = [
            ti / kf(0.5 * (prof[i] + prof[i + 1]))
            for i, (ti, kf) in enumerate(zip(t, kfun, strict=True))
        ]
        if max(abs(a - b) for a, b in zip(new, r, strict=True)) < 1e-18:
            return new
        r = new
    return r


def test_the_fixed_point_derivative_is_exact() -> None:
    """Steep synthetic k(T) in a three-layer stack: dR_wall/dp against
    central differences of the converged fixed point."""
    kfun = [
        lambda T: 1.0 + 2e-3 * T,
        lambda T: 30.0 - 0.01 * T,
        lambda T: 0.5 + 1e-6 * T * T,
    ]
    dkfun = [lambda T: 2e-3, lambda T: -0.01, lambda T: 2e-6 * T]
    t = [3e-4, 2e-3, 1e-3]
    p = [1600.0, 500.0, 900.0, 300.0]
    r = _fixed_point(p, t, kfun)
    prof = cb.wall_temperature_profile(*p, r, 1e-4)[0]
    dr = []
    for i in range(3):
        Tm = 0.5 * (prof[i] + prof[i + 1])
        dr.append(-t[i] * dkfun[i](Tm) / kfun[i](Tm) ** 2)
    got = _solver_tools.wall_resistance_fixed_point_derivatives(*p, r, 1e-4, dr)
    for j in range(4):
        e = 1e-6 * p[j]
        hi, lo = list(p), list(p)
        hi[j] += e
        lo[j] -= e
        fd = (sum(_fixed_point(hi, t, kfun)) - sum(_fixed_point(lo, t, kfun))) / (2 * e)
        assert got[j] == pytest.approx(fd, rel=1e-6), j
    # A constant-k stack has none.
    assert (
        list(_solver_tools.wall_resistance_fixed_point_derivatives(*p, r, 1e-4, [0.0] * 3))
        == [0.0] * 4
    )


@pytest.mark.parametrize(
    "wall",
    [
        {"contact_area": 0.05},
        {"area_a": 0.03, "area_b": 0.009},
    ],
)
def test_the_coupling_derivatives_include_k_of_T(wall: dict) -> None:
    """ss310, and YSZ on hastelloy X: every input
    derivative of compute_coupling against central differences (k held
    fixed misses by 1e-4..5e-3 here)."""
    for layers in (
        [WallLayer(0.002, 20.0, material="ss310")],
        [
            WallLayer(0.0003, 1.0, material="ysz"),
            WallLayer(0.002, 20.0, material="hastelloy_x"),
        ],
    ):
        w = ThermalWall("w", "a", "b", layers=layers, **wall)
        args = [800.0, 1300.0, 0.05, 150.0, 400.0, 0.05, 50.0, 20.0]
        r = w.compute_coupling(*args)
        for i, name in ((0, "dQ_dh_a"), (1, "dQ_dT_aw_a"), (3, "dQ_dh_b"), (4, "dQ_dT_aw_b")):
            e = 1e-5 * abs(args[i])
            hi, lo = list(args), list(args)
            hi[i] += e
            lo[i] -= e
            fd = (w.compute_coupling(*hi).Q - w.compute_coupling(*lo).Q) / (2 * e)
            assert getattr(r, name) == pytest.approx(fd, rel=1e-8), (layers[-1].material, name)


def _pair(layers: list[WallLayer], T_hot: float, T_cold: float, **wall: object) -> FlowNetwork:
    """Hot gas duct coupled to a small coolant duct through a wall."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=T_hot, Y=AIR))
    n.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=T_hot, Y=AIR))
    n.add_node(MassFlowBoundary("c_in", m_dot=0.02, Tt=T_cold, Y=AIR))
    n.add_node(PressureBoundary("c_out", Pt=1.5e5, Tt=T_cold, Y=AIR))
    for eid, a, b, d in (("hot", "h_in", "h_out", 0.04), ("cold", "c_in", "c_out", 0.01)):
        n.add_element(
            ChannelElement(
                eid, a, b, length=0.3, diameter=d, surface=ConvectiveSurface(model=SmoothModel())
            )
        )
    n.add_wall(ThermalWall("w", "hot", "cold", layers=layers, **wall))
    return n


@pytest.mark.parametrize(
    ("layers", "T_hot", "T_cold", "wall"),
    [
        # ss310 below its data (the envelope's low fade), one area / per side
        ([WallLayer(0.002, 20.0, material="ss310")], 900.0, 350.0, {}),
        (
            [WallLayer(0.002, 20.0, material="ss310")],
            900.0,
            350.0,
            {"area_a": 0.03, "area_b": 0.009},
        ),
        # YSZ + bond coat + inconel 718, per side
        (
            [
                WallLayer(0.0003, 1.0, material="ysz"),
                WallLayer(0.0001, 10.0, material="nicocraly"),
                WallLayer(0.002, 20.0, material="inconel718"),
            ],
            1500.0,
            350.0,
            {"area_a": 0.03, "area_b": 0.02},
        ),
        # hastelloy X past its data (the envelope's hot fade, and the
        # linear continuation)
        ([WallLayer(0.002, 20.0, material="hastelloy_x")], 1900.0, 1150.0, {}),
        (
            [WallLayer(0.002, 20.0, material="hastelloy_x", extrapolate=True)],
            1900.0,
            1150.0,
            {},
        ),
    ],
)
def test_a_k_of_T_wall_network_has_exact_rows(
    layers: list[WallLayer], T_hot: float, T_cold: float, wall: dict
) -> None:
    """Every row against central differences, to the constant-k level."""
    s = NetworkSolver(_pair(layers, T_hot, T_cold, **wall))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    FD = np.zeros_like(J)
    for j in range(len(x)):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        FD[:, j] = (F(x + e) - F(x - e)) / (2 * h)
    for i in range(J.shape[0]):
        scale = (np.abs(J[i]) * np.abs(x)).max()
        err = (np.abs(J[i] - FD[i]) * np.abs(x)).max()
        # Constant-k walls reach 4e-8 here (5e-7 at 1900 K, a pressure
        # column's difference noise); the k(T) gap was 1e-5..4e-4.
        assert err <= 1e-6 * scale, (i, s.unknown_names[int(np.argmax(np.abs(J[i] - FD[i])))])
