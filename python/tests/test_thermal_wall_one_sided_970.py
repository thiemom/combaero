"""#970: a ThermalWall with no stream on side B -- an adiabatic side, or a
heat-loss side (a fixed external temperature behind an h).

Adiabatic: the wall passes no heat and sits at side A's face balance (its
T_aw, or with gas radiation the face's radiation balance). HeatLoss: side
A's stream exchanges with T_ext through the series resistance, the
epsilon-NTU integral against a reservoir (C++ ``wall_exchange_and_jacobian``
with C_b = inf); the heat leaves the network, booked as ``{wall}.Q_loss``.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    FlowNetwork,
    GasRadiation,
    HeatLoss,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
M_DOT = 0.02  # kg/s


def _ducts(T_in: float = 900.0, m: float = M_DOT, **wall: object) -> FlowNetwork:
    """Two walled ducts in series through a plenum, each wall one-sided
    (``wall`` the walls' keywords); no walls when ``wall`` is empty."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("in", m_dot=m, Tt=T_in, Y=AIR))
    n.add_node(PlenumNode("mid"))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=T_in, Y=AIR))
    for eid, a, b in (("d1", "in", "mid"), ("d2", "mid", "out")):
        n.add_element(
            ChannelElement(
                eid, a, b, length=0.4, diameter=0.03, surface=ConvectiveSurface(model=SmoothModel())
            )
        )
        if wall:
            n.add_wall(ThermalWall(f"w_{eid}", eid, None, **wall))  # type: ignore[arg-type]
    return n


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return s, r


# --- the API ---------------------------------------------------------------


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"element_b": None}, "needs side_b"),
        ({"element_b": None, "side_b": "insulated"}, "needs side_b"),
        ({"element_b": "b", "side_b": "adiabatic"}, "drop side_b"),
        (
            {
                "element_b": None,
                "side_b": "adiabatic",
                "radiation_b": GasRadiation(beam_length=0.05),
            },
            "no radiation_b",
        ),
        (
            {"element_b": None, "side_b": "adiabatic", "film": "f", "film_side": "b"},
            "film_side must be 'a'",
        ),
        (
            {
                "element_b": None,
                "side_b": "adiabatic",
                "radiation_a": CasingRadiation(area_ratio=0.7),
            },
            "casing that no stream carries",
        ),
    ],
)
def test_a_side_without_a_stream_is_validated(kw: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        ThermalWall("w", "a", **kw)


@pytest.mark.parametrize(
    ("T_ext", "h"), [(0.0, 10.0), (-5.0, 10.0), (300.0, 0.0), (300.0, math.nan)]
)
def test_heat_loss_needs_a_temperature_and_an_h(T_ext: float, h: float) -> None:
    with pytest.raises(ValueError, match="HeatLoss"):
        HeatLoss(T_ext=T_ext, h=h)


def test_the_endpoints() -> None:
    two = ThermalWall("w", "a", "b")
    adia = ThermalWall("w", "a", None, side_b="adiabatic")
    loss = ThermalWall("w", "a", None, side_b=HeatLoss(300.0, 5.0))
    assert two.endpoints == two.heated_endpoints == ("a", "b")
    assert adia.endpoints == ("a",) and adia.heated_endpoints == ()
    assert loss.endpoints == loss.heated_endpoints == ("a",)
    assert adia.adiabatic and not loss.adiabatic and loss.heat_loss == HeatLoss(300.0, 5.0)


def test_an_unknown_element_a_is_still_refused() -> None:
    n = _ducts()
    with pytest.raises(ValueError, match="Unknown element_a"):
        n.add_wall(ThermalWall("w", "nope", None, side_b="adiabatic"))


# --- adiabatic ---------------------------------------------------------------


def test_an_adiabatic_wall_changes_nothing() -> None:
    """The network solves as without the walls, to the bit; each wall passes
    no heat and its whole stack sits at side A's T_aw."""
    _, bare = _solve(_ducts())
    layers = [WallLayer(0.0003, 1.0, material="ysz"), WallLayer(0.002, 20.0, material="ss310")]
    s, r = _solve(_ducts(side_b="adiabatic", layers=layers, contact_area=0.02))
    np.testing.assert_array_equal(r["__x_solution__"], bare["__x_solution__"])
    for w in ("w_d1", "w_d2"):
        assert r[f"{w}.Q"] == 0.0
        assert r[f"{w}.Q_loss"] == 0.0
        assert r[f"{w}.T_interface"] == [r[f"{w}.T_aw_a"]] * 3
        assert r[f"{w}.T_hot"] == r[f"{w}.T_aw_a"]
        assert r["__second_law_walls__"][w]["verdict"] == "pass"
    assert r["__ledger__"]["Q_loss"] == 0.0
    assert not any(k.endswith(".Q_wall_out") for k in r)
    # Not heated: the ducts keep their adiabatic momentum.
    assert not s.network.elements["d1"]._wall_heated


def test_an_adiabatic_wall_with_gas_radiation_sits_at_the_face_balance() -> None:
    """Gas radiation into the face (from a given 2000 K gas over a 900 K
    stream): the wall runs hot of T_aw by what its
    convection must hand back, h (T_s - T_aw) = R1 (Lefebvre Eq. 8.7), and
    the stream still gets no net heat."""
    rad = GasRadiation(
        beam_length=0.05, emissivity_wall=0.7, T_gas=2000.0, fuel_air_ratio=0.03, pressure=1e5
    )
    _, bare = _solve(_ducts())
    _, r = _solve(_ducts(side_b="adiabatic", radiation_a=rad, contact_area=0.03))
    np.testing.assert_array_equal(r["__x_solution__"], bare["__x_solution__"])
    T_s, T_aw, h = r["w_d1.T_surface_a"], r["w_d1.T_aw_a"], r["w_d1.h_a"]
    assert T_s > T_aw + 1.0
    assert r["w_d1.T_interface"][0] == T_s
    R1 = cb.liner_gas_radiation(2000.0, T_s, r["w_d1.eps_g_a"], 0.7).R1
    assert h * (T_s - T_aw) == pytest.approx(R1, rel=1e-9)
    assert r["w_d1.q_rad_a"] == pytest.approx(R1, rel=1e-9)
    assert r["w_d1.Q"] == 0.0


# --- heat loss ---------------------------------------------------------------


def test_heat_loss_is_the_exchange_against_a_reservoir() -> None:
    """Q is the C++ exchange with C_b = inf, Q = C_a (1 - exp(-UA/C_a)) dT,
    at the stream's capacity rate (cp the secant from its inlet to T_ext);
    it leaves the network: the ledger closes on Q_loss, and the wall
    generates entropy, Q (1/T_ext - mean 1/T_stream) > 0."""
    loss = HeatLoss(T_ext=300.0, h=40.0)
    s, r = _solve(_ducts(side_b=loss, layers=[WallLayer(0.003, 1.5)], contact_area=0.04))
    x = np.asarray(r["__x_solution__"])
    for w, feed in (("w_d1", "in"), ("w_d2", "mid")):
        Q = r[f"{w}.Q"]
        assert Q > 0.0 and r[f"{w}.Q_loss"] == Q
        st = s._get_node_state(s.network.nodes[feed], x)
        cp = cb.cp_secant_mass_and_jacobian(float(st.Tt), 300.0, cb.mass_to_mole(AIR))[0]
        C = M_DOT * cp
        UA = r[f"{w}.UA"]
        dT = r[f"{w}.T_aw_a"] - 300.0
        assert pytest.approx(C * (1.0 - math.exp(-UA / C)) * dT, rel=1e-12) == Q
        R = 1.0 / r[f"{w}.h_a"] + 0.003 / 1.5 + 1.0 / 40.0
        assert pytest.approx(0.04 / R, rel=1e-12) == UA
        assert r[f"{w}.NTU"] == pytest.approx(UA / C, rel=1e-12)
        assert r["__second_law_walls__"][w]["verdict"] == "pass"
        assert r["__second_law_walls__"][w]["sigma_wall"] > 0.0
    led = r["__ledger__"]
    assert led["Q_loss"] == pytest.approx(r["w_d1.Q"] + r["w_d2.Q"], rel=1e-12)
    assert abs(led["energy_closure_rel"]) < 1e-9
    # The outlet stream carries what is left.
    h_in = cb.h_mass(900.0, cb.mass_to_mole(AIR))
    h_mid = cb.h_mass(float(r["mid.Tt"]), cb.mass_to_mole(AIR))
    assert M_DOT * (h_in - h_mid) == pytest.approx(r["w_d1.Q"], rel=1e-9)


def test_heat_loss_never_takes_more_than_the_stream_holds() -> None:
    """A huge h and a thin wall: the stream at most reaches T_ext."""
    loss = HeatLoss(T_ext=300.0, h=1e7)
    _, r = _solve(_ducts(m=0.002, side_b=loss, contact_area=5.0))
    assert r["w_d1.effectiveness"] == pytest.approx(1.0, abs=1e-6)
    assert r["mid.Tt"] == pytest.approx(300.0, abs=0.5)
    assert r["mid.Tt"] >= 300.0 - 1e-9


def test_heat_loss_with_per_side_areas_is_the_series_resistance() -> None:
    """1/UA = 1/(h_a A_a) + t/(k A_m) + 1/(h_ext A_b), log-mean layer area."""
    A_a, A_b = 0.03, 0.045
    loss = HeatLoss(T_ext=320.0, h=25.0)
    _, r = _solve(_ducts(side_b=loss, layers=[WallLayer(0.004, 2.0)], area_a=A_a, area_b=A_b))
    A_m = (A_b - A_a) / math.log(A_b / A_a)
    R = 1.0 / (r["w_d1.h_a"] * A_a) + 0.004 / (2.0 * A_m) + 1.0 / (25.0 * A_b)
    assert r["w_d1.UA"] == pytest.approx(1.0 / R, rel=1e-12)


# --- the Jacobian ------------------------------------------------------------


@pytest.mark.parametrize(
    "wall",
    [
        {"side_b": HeatLoss(300.0, 40.0), "layers": [WallLayer(0.003, 1.5)], "contact_area": 0.04},
        {
            "side_b": HeatLoss(350.0, 15.0),
            "layers": [
                WallLayer(0.0003, 1.0, material="ysz"),
                WallLayer(0.002, 20.0, material="ss310"),
            ],
            "area_a": 0.03,
            "area_b": 0.04,
        },
        {
            "side_b": HeatLoss(700.0, 60.0),
            "layers": [WallLayer(0.0012, 26.0)],
            "contact_area": 0.03,
            "radiation_a": GasRadiation(beam_length=0.03, emissivity_wall=0.7),
        },
        {"side_b": "adiabatic", "layers": [WallLayer(0.002, 20.0, material="ss310")]},
    ],
)
def test_one_sided_rows_match_central_differences(wall: dict) -> None:
    """Every row against central differences, scored as |J_ij - FD_ij| |x_j|
    over the row's largest |J_ij x_j| (the project's metric)."""
    T_in = 1600.0 if "radiation_a" in wall else 900.0
    n = _ducts(T_in=T_in, **wall)
    if "radiation_a" in wall:
        # a fuel/air ratio from the composition: burnt gas, not air
        n.nodes["in"].Y = list(cb.mole_to_mass(_burnt()))
        n.nodes["out"].Y = list(n.nodes["in"].Y)
    s, r = _solve(n)
    x = np.array(r["__x_solution__"])
    if wall["side_b"] != "adiabatic":
        assert min(r[f"w_{d}.Q"] for d in ("d1", "d2")) > 0.0
    for x0 in (x, x * (1.0 + 0.01 * np.cos(np.arange(len(x))))):
        _, J = s._residuals_and_jacobian(x0)
        J = J.toarray()

        def F(v: np.ndarray) -> np.ndarray:
            return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

        FD = np.zeros_like(J)
        for j in range(len(x0)):
            h = 1e-6 * max(abs(x0[j]), 1e-3)
            e = np.zeros_like(x0)
            e[j] = h
            FD[:, j] = (F(x0 + e) - F(x0 - e)) / (2 * h)
        for i in range(J.shape[0]):
            scale = (np.abs(J[i]) * np.abs(x0)).max()
            err = (np.abs(J[i] - FD[i]) * np.abs(x0)).max()
            assert err <= 1e-6 * scale, (i, s.unknown_names[int(np.argmax(np.abs(J[i] - FD[i])))])


def _burnt() -> list[float]:
    X = list(cb.species.dry_air())
    names = list(cb.species.names)
    X[names.index("CO2")] += 0.04
    X[names.index("H2O")] += 0.08
    tot = sum(X)
    return [v / tot for v in X]


# --- the kernel ---------------------------------------------------------------


def test_the_adiabatic_face_kernel() -> None:
    """No radiation: T_s = T_aw exactly; gas radiation: between T_aw and Tg."""
    r = _solver_tools.wall_adiabatic_face(812.5, 40.0, 0, 0.0, 0.0)
    assert r.T_s == 812.5 and r.q_rad == 0.0 and list(r.dT_s) == [1.0, 0.0, 0.0, 0.0]
    k = 0.5 * cb.STEFAN_BOLTZMANN * 1.7 * 0.3
    g = _solver_tools.wall_adiabatic_face(1500.0, 300.0, 1, 2000.0, k)
    assert 1500.0 < g.T_s < 2000.0
    assert 300.0 * (g.T_s - 1500.0) == pytest.approx(g.q_rad, rel=1e-12)
