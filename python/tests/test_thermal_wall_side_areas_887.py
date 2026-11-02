"""#887: a ThermalWall with its two faces of different area -- each side's
convective resistance on its own face, the layers on their conduction area
(log mean, exact for a thick cylinder), and its Jacobian."""

from __future__ import annotations

import json
import math
import warnings

import numpy as np
import pytest

import combaero as cb
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
RI, RO, L = 0.01, 0.025, 0.3
A_I, A_O = 2 * math.pi * RI * L, 2 * math.pi * RO * L


def _cylinder(*layers: tuple[float, float], wall_area: str = "log_mean") -> ThermalWall:
    return ThermalWall(
        "w",
        "a",
        "b",
        layers=[WallLayer(thickness=t, conductivity=k) for t, k in layers],
        area_a=A_I,
        area_b=A_O,
        wall_area=wall_area,
    )


def test_a_thick_cylinder_has_the_analytic_heat() -> None:
    """Q = dT / (1/(h_i A_i) + ln(r_o/r_i)/(2 pi k L) + 1/(h_o A_o))."""
    ha, hb, Ta, Tb, k = 800.0, 150.0, 1200.0, 400.0, 20.0
    r = _cylinder((RO - RI, k)).compute_coupling(ha, Ta, 1.0, hb, Tb, 1.0)
    exact = (Ta - Tb) / (
        1 / (ha * A_I) + math.log(RO / RI) / (2 * math.pi * k * L) + 1 / (hb * A_O)
    )
    assert pytest.approx(exact, rel=1e-13) == r.Q
    assert r.T_hot == pytest.approx(Ta - r.Q / (ha * A_I), rel=1e-13)


def test_a_layered_cylinder_is_exact_layer_by_layer() -> None:
    """The area runs linearly with depth through the stack, as with radius."""
    ha, hb, Ta, Tb = 800.0, 150.0, 1200.0, 400.0
    t1, k1, k2 = 0.006, 20.0, 2.0
    r = _cylinder((t1, k1), (RO - RI - t1, k2)).compute_coupling(ha, Ta, 1.0, hb, Tb, 1.0)
    R = math.log((RI + t1) / RI) / (2 * math.pi * k1 * L) + math.log(RO / (RI + t1)) / (
        2 * math.pi * k2 * L
    )
    exact = (Ta - Tb) / (1 / (ha * A_I) + R + 1 / (hb * A_O))
    assert pytest.approx(exact, rel=1e-13) == r.Q


@pytest.mark.parametrize("wall_area", ["log_mean", "slab_a", "slab_b"])
def test_the_coupling_jacobian_is_exact(wall_area: str) -> None:
    w = _cylinder((0.006, 20.0), (0.009, 2.0), wall_area=wall_area)
    args = [800.0, 1200.0, 1.0, 150.0, 400.0, 1.0, 3.0, 0.8]
    r = w.compute_coupling(*args)
    for i, name in ((0, "dQ_dh_a"), (1, "dQ_dT_aw_a"), (3, "dQ_dh_b"), (4, "dQ_dT_aw_b")):
        e = 1e-6 * abs(args[i])
        hi, lo = list(args), list(args)
        hi[i] += e
        lo[i] -= e
        fd = (w.compute_coupling(*hi).Q - w.compute_coupling(*lo).Q) / (2 * e)
        assert getattr(r, name) == pytest.approx(fd, rel=1e-6), name


def test_equal_faces_are_the_one_area_wall() -> None:
    """area_a == area_b gives exactly the one-area wall."""
    args = (800.0, 1200.0, 0.05, 150.0, 400.0, 0.08, 3.0, 0.8)
    layers = [WallLayer(thickness=0.003, conductivity=20.0)]
    one = ThermalWall("w", "a", "b", layers=layers, contact_area=0.05).compute_coupling(*args)
    two = ThermalWall("w", "a", "b", layers=layers, area_a=0.05, area_b=0.05).compute_coupling(
        *args
    )
    for f in ("Q", "T_hot", "dQ_dh_a", "dQ_dh_b", "dQ_dT_aw_a", "dQ_dT_aw_b", "UA"):
        assert getattr(two, f) == pytest.approx(getattr(one, f), rel=1e-14), f


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"area_a": 0.1}, "both area_a and area_b"),
        ({"area_a": 0.1, "area_b": 0.2, "contact_area": 0.1}, "give one"),
        ({"area_a": 0.0, "area_b": 0.2}, "area_a must be > 0"),
        ({"area_a": 0.1, "area_b": float("nan")}, "area_b must be > 0"),
        ({"area_a": 0.1, "area_b": 0.2, "wall_area": "mean"}, "wall_area"),
    ],
)
def test_bad_faces_are_refused(kw: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        ThermalWall("w", "a", "b", layers=[WallLayer(0.003, 20.0)], **kw)


def _pair(material: str = "ss310", **wall: object) -> FlowNetwork:
    """Hot gas duct coupled to a small coolant duct through a wall."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=900.0, Y=AIR))
    n.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=900.0, Y=AIR))
    n.add_node(MassFlowBoundary("c_in", m_dot=0.02, Tt=350.0, Y=AIR))
    n.add_node(PressureBoundary("c_out", Pt=1.5e5, Tt=350.0, Y=AIR))
    for eid, a, b, d in (("hot", "h_in", "h_out", 0.04), ("cold", "c_in", "c_out", 0.01)):
        n.add_element(
            ChannelElement(
                eid, a, b, length=0.3, diameter=d, surface=ConvectiveSurface(model=SmoothModel())
            )
        )
    n.add_wall(
        ThermalWall("w", "hot", "cold", layers=[WallLayer(0.002, 20.0, material=material)], **wall)
    )
    return n


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve()


def test_a_network_with_a_per_side_wall_has_exact_rows() -> None:
    """A constant-k layer: with a k(T) material the wall's Q also moves with
    the wall temperatures through k, which the Jacobian leaves out for one
    area and per-side areas alike (8.6e-5 / 7.1e-5 scaled, measured), so it
    would not isolate what this change adds."""
    s, r = _solve(_pair(material="generic", area_a=0.03, area_b=0.009))
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    _, J = s._residuals_and_jacobian(x)
    J = J.toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    for i in range(J.shape[0]):
        scale = (np.abs(J[i]) * np.abs(x)).max()
        for j, name in enumerate(s.unknown_names):
            h = 1e-6 * max(abs(x[j]), 1e-3)
            e = np.zeros_like(x)
            e[j] = h
            fd = (F(x + e)[i] - F(x - e)[i]) / (2 * h)
            assert J[i, j] * abs(x[j]) == pytest.approx(
                fd * abs(x[j]), rel=1e-5, abs=1e-7 * scale
            ), (
                i,
                name,
            )


def test_a_smaller_cold_face_passes_less_heat() -> None:
    """The coolant side's smaller face adds resistance: less heat than the
    one-area wall on the hot face, more than one on the cold face."""
    q = {}
    for tag, wall in (
        ("hot face", {"contact_area": 0.03}),
        ("per side", {"area_a": 0.03, "area_b": 0.009, "wall_area": "slab_a"}),
        ("cold face", {"contact_area": 0.009}),
    ):
        _, r = _solve(_pair(**wall))
        assert r["__success__"], r["__message__"]
        q[tag] = abs(r["w.Q"])
    assert q["hot face"] > q["per side"] > q["cold face"]


def test_the_gui_builds_and_exports_a_per_side_wall() -> None:
    from gui.backend.export import network_to_gui
    from gui.backend.runner import NetworkRunner

    net = _pair(area_a=0.03, area_b=0.009, wall_area="slab_a")
    doc = network_to_gui(net)
    wall = next(e for e in doc["edges"] if e["id"] == "w")["data"]
    assert (wall["area_a"], wall["area_b"], wall["wall_area"]) == (0.03, 0.009, "slab_a")
    _, direct = _solve(_pair(area_a=0.03, area_b=0.009, wall_area="slab_a"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gui = NetworkRunner.from_dict(json.loads(json.dumps(doc))).solve(timeout=60.0)
    assert gui.success, gui.message
    assert gui._raw["w.Q"] == pytest.approx(direct["w.Q"], rel=1e-6)
