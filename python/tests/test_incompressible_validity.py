"""The incompressible validity report (#588): one Mach limit, report only."""

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)
from combaero.network.compressibility import INCOMPRESSIBLE_MACH_LIMIT

P_OUT = 101325.0
T_IN = 300.0
D = 0.05
AREA = math.pi / 4 * D * D


def _solve(kind: str, mach: float, regime: str = "incompressible", limit: float | None = None):
    rho = float(cb.density(T_IN, P_OUT, cb.species.dry_air()))
    a = float(cb.speed_of_sound(T_IN, cb.species.dry_air()))
    net = FlowNetwork()
    if kind == "channel":
        m = mach * a * rho * AREA
        el = ChannelElement("el", "inlet", "outlet", length=1.0, diameter=D, regime=regime)
    else:
        # mach is the jet's; the bore runs at 0.65x that.
        m = mach * a * rho * 0.65 * AREA
        el = OrificeElement(
            "el", "inlet", "outlet", Cd=0.65, diameter=D, regime=regime, correlation="fixed"
        )
    net.add_node(MassFlowBoundary("inlet", m_dot=m, Tt=T_IN))
    net.add_node(PressureBoundary("outlet", Pt=P_OUT, Tt=T_IN))
    net.add_element(el)
    s = NetworkSolver(net)
    if limit is not None:
        s.incompressible_mach_limit = limit
    r = s.solve()
    assert r["__success__"], r["__message__"]
    return r


def test_the_default_limit_is_mach_0_2() -> None:
    assert INCOMPRESSIBLE_MACH_LIMIT == 0.2


@pytest.mark.parametrize("kind", ["channel", "orifice"])
def test_below_the_limit_is_ok_and_above_it_warns(kind: str) -> None:
    low = _solve(kind, 0.1)
    v = low["__incompressible_validity__"]["el"]
    assert v["verdict"] == "ok" and v["limit"] == 0.2
    assert low["__incompressible_warnings__"] == []

    high = _solve(kind, 0.3)
    v = high["__incompressible_validity__"]["el"]
    assert v["verdict"] == "warn" and v["mach"] > 0.2
    assert high["__incompressible_warnings__"] == ["el"]
    # A report: the solve still succeeds, with the same solution it would
    # have had without the check.
    assert high["__success__"]


def test_the_channel_reports_its_higher_end() -> None:
    r = _solve("channel", 0.3)
    d = r["__element_diag__"]["el"]
    assert r["__incompressible_validity__"]["el"]["mach"] == pytest.approx(
        max(d["mach_in"], d["mach_out"])
    )


def test_the_orifice_is_judged_by_its_jet_not_its_bore() -> None:
    """Bore Mach ~0.16, vena contracta ~0.25: the jet decides."""
    r = _solve("orifice", 0.25)
    d = r["__element_diag__"]["el"]
    assert max(d["mach_in"], d["mach_out"]) < 0.2
    v = r["__incompressible_validity__"]["el"]
    assert v["mach"] == pytest.approx(d["mach_throat"])
    assert v["verdict"] == "warn"


def test_compressible_elements_are_not_reported() -> None:
    r = _solve("channel", 0.3, regime="compressible")
    assert r["__incompressible_validity__"] == {}
    assert r["__incompressible_warnings__"] == []


def test_the_limit_is_a_setting() -> None:
    r = _solve("channel", 0.3, limit=0.5)
    v = r["__incompressible_validity__"]["el"]
    assert v["verdict"] == "ok" and v["limit"] == 0.5
