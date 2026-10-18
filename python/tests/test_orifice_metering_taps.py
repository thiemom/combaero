"""ISO 5167 metering Cd at a chosen tap, and the orifice drive with a velocity
of approach.

Three defects, each pinned here:

1. Reader-Harris/Gallagher was flange-tap only and lacked ISO 5167-2:2003's
   small-pipe term, +0.011 (0.75 - beta)(2.8 - D/25.4) below D = 71.12 mm.
2. The network's 'Stolz' was a four-term Reader-Harris/Gallagher typed out in
   Python, and its 'Miller' was Stolz; C++'s Miller had no published source.
   Metering Cd now comes from C++'s metering_cd_and_dRe only; Miller is gone.
3. The orifice kernels drove from the upstream STAGNATION pressure and also
   applied E = 1/sqrt(1 - beta^4), counting the approach velocity head twice
   behind a node that carries it: +6.2% flow at beta 0.7 behind a
   MomentumChamberNode. With E applied the drive is now the upstream static.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

M = cb.MeteringCdCorrelation
T = cb.OrificeTap


def _geom(beta: float, D: float, tap: object) -> cb.OrificeGeometry:
    g = cb.OrificeGeometry()
    g.D = D
    g.d = beta * D
    g.tap = tap
    return g


def _iso_rhg(beta: float, Re: float, D: float, L1: float, L2: float) -> float:
    """ISO 5167-2:2003 Eq. (4) with the small-pipe term, typed independently
    of the C++ (the constants are the standard's, not include/orifice.h's)."""
    A = (19000.0 * beta / Re) ** 0.8
    M2 = 2.0 * L2 / (1.0 - beta)
    C = (
        0.5961
        + 0.0261 * beta**2
        - 0.216 * beta**8
        + 0.000521 * (1e6 * beta / Re) ** 0.7
        + (0.0188 + 0.0063 * A) * beta**3.5 * (1e6 / Re) ** 0.3
        + (0.043 + 0.080 * math.exp(-10.0 * L1) - 0.123 * math.exp(-7.0 * L1))
        * (1.0 - 0.11 * A)
        * beta**4
        / (1.0 - beta**4)
        - 0.031 * (M2 - 0.8 * M2**1.1) * beta**1.3
    )
    if D * 1000.0 < 71.12:
        C += 0.011 * (0.75 - beta) * (2.8 - D * 1000.0 / 25.4)
    return C


def _taps(D: float) -> list[tuple[object, tuple[float, float]]]:
    flange = 25.4 / (D * 1000.0)
    return [(T.Corner, (0.0, 0.0)), (T.D_D2, (1.0, 0.47)), (T.Flange, (flange, flange))]


@pytest.mark.parametrize("D", [0.05, 0.07, 0.1, 0.3])
@pytest.mark.parametrize("beta", [0.2, 0.5, 0.7])
@pytest.mark.parametrize("Re", [5e3, 1e5, 3e6])
def test_rhg_is_iso_eq4_at_every_tap(D: float, beta: float, Re: float) -> None:
    for tap, (L1, L2) in _taps(D):
        cd, _ = cb.metering_cd_and_dRe(M.ReaderHarrisGallagher, _geom(beta, D, tap), Re)
        assert cd == pytest.approx(_iso_rhg(beta, Re, D, L1, L2), abs=1e-12), tap


@pytest.mark.parametrize("corr", [M.ReaderHarrisGallagher, M.Stolz])
@pytest.mark.parametrize("D", [0.05, 0.1])
@pytest.mark.parametrize("Re", [5e3, 1e5, 3e6])
def test_the_re_derivative_is_exact_at_every_tap(corr: object, D: float, Re: float) -> None:
    for tap, _ in _taps(D):
        g = _geom(0.6, D, tap)
        _, dcd = cb.metering_cd_and_dRe(corr, g, Re)
        h = 1e-6 * Re
        fd = (
            cb.metering_cd_and_dRe(corr, g, Re + h)[0] - cb.metering_cd_and_dRe(corr, g, Re - h)[0]
        ) / (2 * h)
        assert dcd == pytest.approx(fd, rel=1e-6), tap
        assert dcd != 0.0


def test_the_small_pipe_term_is_present_below_71_mm_and_continuous_there() -> None:
    def cd(D: float) -> float:
        return cb.metering_cd_and_dRe(M.ReaderHarrisGallagher, _geom(0.4, D, T.Corner), 1e5)[0]

    # Corner taps carry no D, so the whole step is the small-pipe term.
    assert cd(0.05) - cd(0.2) == pytest.approx(0.011 * (0.75 - 0.4) * (2.8 - 50.0 / 25.4), rel=1e-9)
    # 2.8 in is exactly 71.12 mm, so the term vanishes where it switches off.
    assert cd(0.07112 * (1 - 1e-9)) == pytest.approx(cd(0.07112 * (1 + 1e-9)), abs=1e-10)


def test_the_flange_default_keeps_the_cd_the_network_already_used() -> None:
    """The tap is a measurement choice; the default stays flange, so no
    network orifice above 71.12 mm moves."""
    g = cb.OrificeGeometry()
    g.d, g.D = 0.05, 0.1
    assert g.tap == T.Flange
    s = cb.OrificeState()
    s.Re_D = 1e5
    L = 25.4 / 100.0
    assert cb.Cd_sharp_thin_plate(g, s) == pytest.approx(_iso_rhg(0.5, 1e5, 0.1, L, L), abs=1e-12)


def test_corner_taps_cross_check_cano_2021() -> None:
    """Cano et al. (2021), Rev. Bras. Eng. Agric. Ambient. 25(1):10-16, Table 2:
    their RHG Cd for beta 0.5 corner-tap plates in 100 and 150 mm water pipes
    averaged 0.606 and 0.605 over Re_D ~ 5e4-3.7e5. A third party evaluating
    the same formula, so an implementation cross-check, not validation."""
    for D, mean in ((0.1, 0.606), (0.15, 0.605)):
        cds = [
            cb.metering_cd_and_dRe(M.ReaderHarrisGallagher, _geom(0.5, D, T.Corner), Re)[0]
            for Re in (5e4, 1e5, 2e5, 3.7e5)
        ]
        assert sum(cds) / len(cds) == pytest.approx(mean, abs=1.5e-3)


def test_stolz_tap_terms_saturate_at_l1_0_4333() -> None:
    """ISO 5167:1980: 0.09 L1 is replaced by 0.039 from L1 = 0.4333 on, which is
    continuous there. D and D/2 taps (L1 = 1) sit on the saturated branch."""
    D = 0.1
    corner = cb.metering_cd_and_dRe(M.Stolz, _geom(0.5, D, T.Corner), 1e5)[0]
    dd2 = cb.metering_cd_and_dRe(M.Stolz, _geom(0.5, D, T.D_D2), 1e5)[0]
    b4 = 0.5**4
    assert dd2 - corner == pytest.approx(0.039 * b4 / (1 - b4) - 0.0337 * 0.47 * 0.5**3, abs=1e-12)


def test_iso_range_flags() -> None:
    rhg = M.ReaderHarrisGallagher
    assert cb.metering_cd_in_range(rhg, _geom(0.5, 0.1, T.Flange), 1e5)
    assert not cb.metering_cd_in_range(rhg, _geom(0.5, 0.04, T.Flange), 1e5)  # D < 50 mm
    assert not cb.metering_cd_in_range(rhg, _geom(0.8, 0.1, T.Flange), 1e5)  # beta > 0.75
    assert not cb.metering_cd_in_range(rhg, _geom(0.5, 0.1, T.Corner), 4e3)  # Re < 5000
    # corner and D-D/2 above beta 0.56: Re_D >= 16000 beta^2 = 7840 at beta 0.7
    assert not cb.metering_cd_in_range(rhg, _geom(0.7, 0.1, T.Corner), 7800)
    assert cb.metering_cd_in_range(rhg, _geom(0.7, 0.1, T.Corner), 7900)
    # flange: Re_D >= 170 beta^2 D[mm] = 41650 at beta 0.7, D 500 mm
    assert not cb.metering_cd_in_range(rhg, _geom(0.7, 0.5, T.Flange), 4.1e4)
    assert cb.metering_cd_in_range(rhg, _geom(0.7, 0.5, T.Flange), 4.2e4)


def test_a_cd_without_a_formula_is_refused() -> None:
    for corr in (M.Constant, M.UserFunction):
        with pytest.raises(ValueError, match="only ReaderHarrisGallagher and Stolz"):
            cb.metering_cd_and_dRe(corr, _geom(0.5, 0.1, T.Flange), 1e5)


def test_miller_is_gone() -> None:
    assert not hasattr(M, "Miller")
    with pytest.raises(ValueError, match="'Miller' correlation was removed"):
        OrificeElement("o", "a", "b", diameter=0.02, correlation="Miller")


def _pipe_orifice_sink(node_cls: type, beta: float, regime: str, correlation: str) -> FlowNetwork:
    """pipe -> a -> orifice -> sink. No pipe beyond the orifice, so no
    pressure recovery (#645): the row is ISO's flow equation on the tapping
    differential, the case this file's drive rule is about. Between pipes
    see test_orifice_permanent_loss.py."""
    D = 0.1
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=300.0))
    g.add_node(node_cls("a"))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_element(ChannelElement("p1", "hi", "a", length=0.5, diameter=D, regime=regime))
    g.add_element(
        OrificeElement(
            "o", "a", "ho", Cd=0.6, diameter=beta * D, regime=regime, correlation=correlation
        )
    )
    return g


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve()
    assert r["__success__"], r["__message__"]
    return r


@pytest.mark.parametrize("node_cls", [PlenumNode, MomentumChamberNode])
@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_a_pipe_fed_orifice_flows_iso_from_its_static_taps(node_cls: type, beta: float) -> None:
    """m = Cd E A sqrt(2 rho (p_a - p_sink)) behind either node. Behind a
    MomentumChamberNode (Pt > P) the flow was high (1.2% at beta 0.5 and 6.2%
    at beta 0.7 between pipes): the stagnation drive counted the approach
    head again."""
    r = _solve(_pipe_orifice_sink(node_cls, beta, "incompressible", "fixed"))
    d = beta * 0.1
    # p1 and rho1 are the approach pipe's static state (#651).
    p1 = r["o.P_approach"]
    rho = cb.density(300.0, p1, cb.species.dry_air())
    E = 1.0 / math.sqrt(1.0 - beta**4)
    m_iso = 0.6 * E * math.pi * d * d / 4 * math.sqrt(2.0 * rho * (p1 - 1.0e5))
    # 1e-4: the network's air vs species.dry_air; the defect was 1.2-6%.
    assert r["o.m_dot"] == pytest.approx(m_iso, rel=1e-4)


def test_a_saved_gui_file_naming_miller_says_what_to_pick() -> None:
    from pydantic import ValidationError

    from gui.backend.schemas import OrificeData

    with pytest.raises(ValidationError, match="Pick 'Stolz' or 'ReaderHarrisGallagher'"):
        OrificeData(correlation="Miller")
    assert OrificeData(correlation="Stolz").correlation == "Stolz"


def _gui_orifice(data: dict) -> OrificeElement:
    from gui.backend.graph_builder import build_network_from_schema
    from gui.backend.schemas import NetworkGraphSchema

    def node(i: str, t: str, d: dict, x: int) -> dict:
        return {"id": i, "type": t, "data": {"id": i, **d}, "position": {"x": x, "y": 0}}

    def edge(i: str, s: str, t: str) -> dict:
        return {
            "id": i,
            "source": s,
            "sourceHandle": "flow-source",
            "target": t,
            "targetHandle": "flow-target",
            "data": {},
        }

    graph = {
        "nodes": [
            node("pb_in", "pressure_boundary", {"Pt": 2e5, "Tt": 300}, 0),
            node("ori", "orifice", {"diameter": 0.01, **data}, 100),
            node("pb_out", "pressure_boundary", {"Pt": 1e5, "Tt": 300}, 200),
        ],
        "edges": [edge("e1", "pb_in", "ori"), edge("e2", "ori", "pb_out")],
    }
    net = build_network_from_schema(NetworkGraphSchema(**graph))
    return net.elements["ori"]


def test_the_gui_passes_the_tap_and_the_bevel_depth() -> None:
    """The GUI's tap picks the reported differential, and its bevel depth
    reaches IdelchikBeveled (it was dropped by the builder before)."""
    o = _gui_orifice({"correlation": "ReaderHarrisGallagher", "tap": "Corner"})
    assert o.tap == "Corner"
    assert _gui_orifice({}).tap == "Flange"
    o = _gui_orifice({"correlation": "IdelchikBeveled", "bevel_depth": 0.002})
    assert o.bevel_depth == 0.002


def test_diagnostics_report_the_approach_along_the_flow() -> None:
    """beta stays the declared approach; beta_flow and E follow the flow, so a
    reversed orifice reports the to-node pipe's (#554)."""
    D_up, D_dn, d = 0.1, 0.08, 0.05
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.0e5, Tt=300.0))
    g.add_node(PlenumNode("a"))
    g.add_node(PlenumNode("b"))
    g.add_node(PressureBoundary("ho", Pt=1.2e5, Tt=300.0))
    g.add_element(ChannelElement("p1", "hi", "a", length=0.5, diameter=D_up))
    g.add_element(OrificeElement("o", "a", "b", Cd=0.6, diameter=d, correlation="fixed"))
    g.add_element(ChannelElement("p2", "b", "ho", length=0.5, diameter=D_dn))
    r = _solve(g)
    assert r["o.m_dot"] < 0.0
    assert r["o.beta"] == pytest.approx(d / D_up, rel=1e-14)
    assert r["o.beta_flow"] == pytest.approx(d / D_dn, rel=1e-14)
    assert r["o.E"] == pytest.approx(cb.velocity_of_approach_factor(d / D_dn), rel=1e-14)
