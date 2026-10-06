"""Default convective areas for impingement channels, and the GUI's working example (#462).

Since #460 an impingement element's mass flow IS the jet flow and the convective
area counts the holes, so the default area decides the hole count. A wetted
wall (pi D L) on a 1 m channel spread a default flow over about 820 holes and
put Re_j far below Florschuetz's range. The defaults now follow the model.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    ImpingementModel,
    MassFlowBoundary,
    PressureBoundary,
    SingleJetImpingementModel,
    SmoothModel,
)

Y_AIR = cb.species.from_mapping({"N2": 0.767, "O2": 0.233})


def _channel(model, D: float = 0.05, L: float = 1.0) -> ChannelElement:
    return ChannelElement(
        "ch", "in", "out", length=L, diameter=D, surface=ConvectiveSurface(model=model)
    )


def test_array_row_footprint_matches_the_channel_cross_section() -> None:
    el = _channel(ImpingementModel(d_jet=0.002, xn_d=8.0, yn_d=6.0, z_d=2.0))
    assert el.default_convective_area() == pytest.approx(el.area * 8.0 / 2.0, rel=1e-12)


def test_single_jet_area_is_goldsteins_averaging_disc() -> None:
    el = _channel(SingleJetImpingementModel(d_jet=0.003, R_D=5.0))
    assert el.default_convective_area() == pytest.approx(math.pi * (5.0 * 0.003) ** 2)


def test_other_surfaces_keep_the_wetted_wall() -> None:
    el = _channel(SmoothModel(), D=0.05, L=1.0)
    assert el.default_convective_area() == pytest.approx(math.pi * 0.05 * 1.0)


def test_an_explicit_area_is_never_overridden() -> None:
    el = ChannelElement(
        "ch",
        "in",
        "out",
        length=1.0,
        diameter=0.05,
        surface=ConvectiveSurface(area=0.123, model=ImpingementModel(d_jet=0.002)),
    )
    el._fill_convective_area()
    assert el.surface.area == 0.123


def _solve(model, m_dot: float = 0.01) -> dict:
    """GUI-shaped network: a default 1 m x 0.05 m channel, boundaries either side."""
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("in", m_dot=m_dot, Tt=500.0, Y=Y_AIR))
    net.add_node(PressureBoundary("out", Pt=3.0e5, Tt=500.0, Y=Y_AIR))
    net.add_element(_channel(model))
    res = NetworkSolver(net).solve()
    assert res["__success__"], res.get("__message__")
    return res["__element_diag__"]["ch"]


def test_the_gui_default_array_is_a_working_example() -> None:
    """The GUI's default impingement surface (d 2 mm, xn/d 8, yn/d 6, z/d 2) on
    a default channel lands inside Florschuetz's Re_j range and says so."""
    diag = _solve(ImpingementModel(d_jet=0.002, xn_d=8.0, yn_d=6.0, z_d=2.0, row=1))
    assert 2.5e3 < diag["Re_surface"] < 7.0e4
    assert diag["surface_extrapolated"] == 0.0


def test_the_single_jet_reports_its_reynolds_number() -> None:
    """One jet carries the whole flow; the GUI shows its Re so a user can see
    when the jet diameter and the flow do not match the source."""
    diag = _solve(SingleJetImpingementModel(d_jet=0.003, L_D=7.75, R_D=5.0))
    mu = cb.complete_state(500.0, 3.0e5, cb.mass_to_mole(Y_AIR)).transport.mu
    assert diag["Re_surface"] == pytest.approx(4.0 * 0.01 / (math.pi * 0.003 * mu), rel=0.02)
