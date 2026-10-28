"""A bypass at another temperature: each impingement row's reference
temperature through CR-3936's eta_r (#467).

    T_aw = T_j + eta_r(Gc/Gj, z/d) phi_0 (T_0 - T_j)

T_j the supply, T_0 the bypass node, phi_0 the bypass share of the crossflow
approaching the row, eta_r ``cb.impingement_eta_r`` -- a closure FITTED to
Florschuetz and Su (1985) Appendix A (see test_impingement_eta_r_data.py
for how it fits). Only the bypass's share: the wall's own heating of the
spent air is inside Florschuetz's (1981) h already.
"""

from __future__ import annotations

import numpy as np
import pytest
from test_impingement_bypass import (
    SPAN,
    Y_AIR,
    D,
    _array,
    _jacobian_error,
    _pb,
    _set_regime,
    _solve,
    _with_bypass,
)

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    PlenumNode,
    WallLayer,
)


def _rows(diag: dict, arr, key: str) -> list[float]:
    return [float(diag[arr.plate_id(i)][key]) for i in range(1, int(arr.n_rows) + 1)]


def test_a_bypass_at_the_supply_temperature_changes_nothing() -> None:
    """T_aw is the supply's exactly, eta_r reported but without effect."""
    net, arr = _with_bypass(0.004)
    diag = _solve(net)["__element_diag__"]
    assert _rows(diag, arr, "T_aw") == pytest.approx([300.0] * 10, abs=1e-9)
    assert all(e > 0.0 for e in _rows(diag, arr, "eta_r"))


def test_a_hot_bypass_moves_each_row_by_its_share() -> None:
    """Row by row the composition the docstring states: row 1's crossflow
    is all bypass (phi_0 = 1); downstream the spent jets dilute it, so the
    reference falls toward the supply."""
    net, arr = _with_bypass(0.004, Tt=400.0)
    res = _solve(net)
    diag = res["__element_diag__"]
    T0 = float(res["byp.Tt"])
    m_b = float(diag[arr.segment_id(0)]["m_dot"])
    T_aw, eta, phi = (_rows(diag, arr, k) for k in ("T_aw", "eta_r", "bypass_fraction"))
    m_c = [m_b] + [float(diag[arr.segment_id(i)]["m_dot"]) for i in range(1, 10)]
    for i in range(10):
        assert phi[i] == pytest.approx(m_b / m_c[i], rel=1e-9)
        gc = float(diag[arr.plate_id(i + 1)]["Gc_Gj"])
        assert eta[i] == pytest.approx(cb.impingement_eta_r(gc, 5.0, 4.0, 2.0).eta, rel=1e-12)
        assert T_aw[i] == pytest.approx(300.0 + eta[i] * phi[i] * (T0 - 300.0), abs=1e-9)
    assert phi[0] == pytest.approx(1.0)
    assert all(a > b for a, b in zip(T_aw, T_aw[1:], strict=False))
    assert 300.0 < T_aw[-1] < T_aw[0] < T0


def test_outside_the_inline_data_is_flagged() -> None:
    """A staggered array, or a gap taller than the data's z/d 3, with a hot
    bypass: eta_r is extrapolated there and the summary says so."""
    for kw in ({"pattern": "staggered"}, {"z_d": 4.0}):
        net, arr = _with_bypass(0.004, Tt=400.0, **kw) if "pattern" in kw else _z4()
        s = arr.summarize(_solve(net)["__element_diag__"])
        assert s["heat_transfer_outside_source"] == 1.0, kw
    net, arr = _with_bypass(0.004, Tt=400.0)
    assert arr.summarize(_solve(net)["__element_diag__"])["heat_transfer_outside_source"] == 0.0


def _z4():
    from combaero.network import ImpingementArray, MassFlowBoundary, MomentumChamberNode

    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(MassFlowBoundary("feed_in", m_dot=0.004, Tt=400.0, Y=Y_AIR))
    net.add_node(MomentumChamberNode("byp"))
    net.add_element(ChannelElement("feed", "feed_in", "byp", length=0.2, diameter=0.02))
    arr = ImpingementArray(
        "ia", n_rows=10, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=4.0, span=SPAN, plate_thickness=D
    )
    arr.add_to(net, "supply", "exit", bypass="byp")
    return net, arr


def _hot_fed(wall: bool) -> FlowNetwork:
    """A pressure-fed hot bypass through a duct into a plenum, so T_0 is a
    derived node temperature, and its flow an unknown; with a hot-gas wall
    the reference temperature enters the energy balance."""
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    net.add_node(_pb("bsrc", 1.02e5, 420.0))
    net.add_node(PlenumNode("byp"))
    net.add_element(ChannelElement("feed", "bsrc", "byp", length=0.1, diameter=0.02))
    arr = _array(n_rows=5)
    arr.add_to(net, "supply", "exit", bypass="byp")
    if wall:
        net.add_node(_pb("hot_in", 1.02e5, 1200.0))
        net.add_node(_pb("hot_out", 1.0e5, 1200.0))
        L = 5 * 5.0 * D
        net.add_element(
            ChannelElement(
                "hot",
                "hot_in",
                "hot_out",
                length=L,
                diameter=0.05,
                surface=ConvectiveSurface(area=L * SPAN),
            )
        )
        arr.add_wall(net, "w", "hot", [WallLayer(thickness=0.001, conductivity=20.0)])
    return net


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_reference_temperature_is_in_the_jacobian(regime: str) -> None:
    """The wall's heat now moves with the bypass share (neighbour flows) and
    the bypass node's temperature: both relayed, analytic. Measured
    2026-10-23: incompressible 2.7e-9, compressible 4.7e-5 (the Fanno
    iteration's repeatability, as in test_impingement_bypass); with the
    relay removed it fails at 0.26 incompressible, 0.86 compressible."""
    net = _set_regime(_hot_fed(wall=True), regime)
    err, missing, s, x, *_ = _jacobian_error(net)
    tol = 1e-6 if regime == "incompressible" else 1e-4
    assert missing == 0
    assert err < tol, err
    # The case exercises what it claims: a hot bypass reaching the rows.
    res = s.solve(init_strategy="outletref_warmstart")
    d = res["__element_diag__"]
    assert float(res["byp.Tt"]) > 400.0
    assert float(d["ia__p1"]["T_aw"]) > 330.0
    assert np.isfinite(float(d["ia__p5"]["T_aw"]))
