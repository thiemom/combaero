"""#895: every convective surface reports its Stanton number St = Nu/(Re Pr),
at the correlation's own Re and Pr, so St is the h it reports made
dimensionless: St = h / (G cp) for a smooth channel."""

from __future__ import annotations

import warnings

import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
)
from combaero.network.schema import HEAT_TRANSFER_FIELDS

Y = list(cb.mole_to_mass(cb.species.dry_air()))


def _channel(surface: ConvectiveSurface | None) -> dict:
    net = FlowNetwork()
    net.add_node(PressureBoundary("a", Pt=2.0e5, Tt=600.0, Y=Y))
    net.add_node(PressureBoundary("b", Pt=1.9e5, Tt=600.0, Y=Y))
    net.add_element(
        ChannelElement("d", "a", "b", length=0.5, diameter=0.02, roughness=0.0, surface=surface)
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return NetworkSolver(net).solve()


def test_st_is_a_heat_transfer_field() -> None:
    assert "St" in HEAT_TRANSFER_FIELDS


def test_smooth_channel_st_is_nu_over_re_pr() -> None:
    sol = _channel(ConvectiveSurface(model=SmoothModel()))
    st = sol["d.St"]
    assert st > 0.0
    # Gnielinski at the surface's own Re and Pr; the reported Pr_ref is at the
    # element's reference state, a hair away from the correlation's (5e-5).
    assert st == pytest.approx(sol["d.Nu"] / (sol["d.Re_surface"] * sol["d.Pr_ref"]), rel=1e-3)
    # The same thing as h / (G cp), G = m / A, to the same 0.1%: cp_ref is at
    # the reference state, the correlation's cp at its own.
    area = 3.141592653589793 * 0.01**2
    assert st == pytest.approx(sol["d.htc"] / (sol["d.m_dot"] / area * sol["d.cp_ref"]), rel=2e-3)
    # A turbulent smooth tube sits near St ~ 0.0023 Re^-0.2 Pr^-0.6 (order of magnitude).
    assert 1e-4 < st < 1e-2


def test_a_bare_channel_reports_the_st_of_the_h_it_reports() -> None:
    """Without a declared surface a channel still reports its smooth-tube Nu
    and h (as diagnostics); its St is theirs."""
    sol = _channel(None)
    assert sol["d.Nu"] > 0.0
    assert sol["d.St"] == pytest.approx(
        sol["d.Nu"] / (sol["d.Re_surface"] * sol["d.Pr_ref"]), rel=1e-3
    )
