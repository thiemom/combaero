"""A chamber or combustor whose heat transfer is not evaluated reports NaN.

With no convective surface (since #1022 a node no longer takes its flow area
as one), MomentumChamberNode and CombustorNode diagnostics reported Nu = St =
htc = 0 and T_aw = Tt. The GUI showed that as "htc 0.0 W/m2/K" (zero heat
transfer) and "T_aw = Tt" (a recovery temperature). The keys stay (#995); the
values are NaN, which the GUI shows as no value. A walled node still reports
its finite htc.
"""

from __future__ import annotations

import math

import pytest

from combaero.network import (
    ChannelElement,
    CombustorNode,
    ConvectiveSurface,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.components import SmoothModel

KEYS = ("Nu", "St", "htc", "T_aw")


def _solve(kind: str, surface: ConvectiveSurface | None) -> dict:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hi", Pt=1.6e5, Tt=600.0))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    if kind == "combustor":
        node = CombustorNode("hn", method="complete", area=3.1416e-4, Dh=0.02, surface=surface)
    else:
        node = MomentumChamberNode("hn", area=3.1416e-4, Dh=0.02, surface=surface)
    g.add_node(node)
    for eid, a, b in (("hin", "hi", "hn"), ("hout", "hn", "ho")):
        g.add_element(ChannelElement(eid, a, b, length=0.3, diameter=0.02, regime="compressible"))
    res = NetworkSolver(g).solve()
    assert res["__success__"], res.get("__message__")
    return res


@pytest.mark.parametrize("kind", ["combustor", "chamber"])
def test_an_unwalled_node_reports_no_heat_transfer_value(kind: str) -> None:
    res = _solve(kind, None)
    for key in KEYS:
        name = f"hn.{key}"
        assert name in res, f"{name} missing: keys follow the network (#995)"
        assert math.isnan(res[name]), f"{name} = {res[name]}: a 0 reads as zero heat transfer"


@pytest.mark.parametrize("kind", ["combustor", "chamber"])
def test_a_walled_node_still_reports_its_htc(kind: str) -> None:
    res = _solve(kind, ConvectiveSurface(model=SmoothModel(), area=0.05))
    assert res["hn.htc"] > 0.0
    assert res["hn.Nu"] > 0.0
    assert math.isfinite(res["hn.T_aw"])
