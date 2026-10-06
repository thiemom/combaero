"""Jet impingement in ConvectiveSurface: mass-flow bookkeeping and wall coupling.

The correlation library (impingement_correlation.h) only knows Re_j, Gc/Gj and
geometry -- it has no concept of "total mass flow through this element's own
area" or "how many holes are in this row". That bookkeeping lives here, in
ImpingementModel/SingleJetImpingementModel's mass-flow conventions, and is
exactly the kind of composition the model-provenance discipline says needs its
own check: the C++ derivatives (dNu_dRe_j, dNu_dRe) are already
finite-difference verified in tests/test_impingement_correlation.cpp, but the
CHAIN from mdot -> Re_j -> h -> dh_dmdot is new code, not covered there.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network.components import (
    ConvectiveSurface,
    ImpingementModel,
    SingleJetImpingementModel,
)

X_AIR = cb.species.dry_air()
STATE = {"T": 500.0, "P": 3.0e5, "X": X_AIR}
FLOW = {"velocity": 5.0, "diameter": 0.02, "length": 0.05, "T_hot": 900.0}


def _array_surface(**kw) -> ConvectiveSurface:
    model_kw = {"d_jet": 0.002, "xn_d": 8.0, "yn_d": 6.0, "z_d": 2.0, "row": 3}
    model_kw.update(kw)
    # 0.0015 m^2 is about 8 holes at this pitch: the element's flow spread
    # over them gives Re_j near 1e4, inside Florschuetz's range (#460).
    return ConvectiveSurface(area=0.0015, model=ImpingementModel(**model_kw))


def _single_jet_surface(**kw) -> ConvectiveSurface:
    model_kw = {"d_jet": 0.003, "L_D": 7.75, "R_D": 5.0}
    model_kw.update(kw)
    return ConvectiveSurface(area=0.001, model=SingleJetImpingementModel(**model_kw))


def _channel_rho_a() -> float:
    """rho * A_cross: the solver's mass flow per unit velocity."""
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    return rho * math.pi / 4.0 * FLOW["diameter"] ** 2


def _run(surface: ConvectiveSurface, **flow_overrides):
    flow = {**FLOW, **flow_overrides}
    return surface.htc_and_T(**STATE, **flow)


def test_array_default_correlation_set_is_inline() -> None:
    r = _run(_array_surface())
    assert r.h > 0.0
    assert r.Nu > 0.0
    assert r.extrapolated is False


def test_array_first_row_sees_zero_crossflow_and_is_the_least_degraded() -> None:
    """Nu1's own definition: Gc/Gj = 0 at row 1. Every later row must be worse."""
    nu_by_row = {row: _run(_array_surface(row=row)).Nu for row in (1, 3, 6, 10)}
    assert nu_by_row[1] > nu_by_row[3] > nu_by_row[6] > nu_by_row[10]


def test_array_staggered_set_is_selectable() -> None:
    staggered = cb.florschuetz_1981_staggered()
    r = _run(_array_surface(correlation_set=staggered))
    assert r.h > 0.0


def test_array_flags_extrapolation_outside_validity() -> None:
    """A geometry far outside Florschuetz's stated ranges (item 20) should be
    flagged, the same way the underlying correlation flags it."""
    r = _run(_array_surface(xn_d=1.0, yn_d=1.0, z_d=0.1))
    assert r.extrapolated is True


def test_array_h_scales_with_total_mass_flow() -> None:
    """More coolant through the same row area means a faster jet, means a
    higher local Re_j and a higher h -- the whole point of recovering mdot
    from velocity*area rather than treating velocity as the jet's own."""
    low = _run(_array_surface(), velocity=2.0)
    high = _run(_array_surface(), velocity=8.0)
    assert high.h > low.h
    assert high.Re > low.Re


def test_array_wall_coupling_derivative_matches_finite_difference() -> None:
    """The mdot -> Re_j -> h chain is new composition, not covered by the
    C++ derivative test -- falsify it independently here."""
    surface = _array_surface()
    v0 = 5.0
    dv = 1e-4

    r0 = _run(surface, velocity=v0)
    r_plus = _run(surface, velocity=v0 + dv)
    r_minus = _run(surface, velocity=v0 - dv)

    # dh/dmdot is with respect to the SOLVER's mass flow, through the channel
    # cross-section rho v pi D^2/4 -- not the convective area. The first
    # version of this test converted through the convective area too, so it
    # agreed with the bug it should have caught (#456).
    dh_dv_fd = (r_plus.h - r_minus.h) / (2.0 * dv)
    assert r0.dh_dmdot == pytest.approx(dh_dv_fd / _channel_rho_a(), rel=2e-3)


def test_array_reports_extra_fields_a_channel_result_does_not_have() -> None:
    r = _run(_array_surface())
    for attr in ("dh_dmdot", "dh_dT", "dT_aw_dmdot", "dT_aw_dT", "extrapolated"):
        assert hasattr(r, attr)


def test_single_jet_default_correlation_set() -> None:
    r = _run(_single_jet_surface())
    assert r.h > 0.0
    assert r.Nu > 0.0


def test_single_jet_boundary_condition_changes_the_result() -> None:
    heat_flux = _run(_single_jet_surface(bc=cb.ImpingementThermalBC.ConstantHeatFlux))
    wall_temp = _run(_single_jet_surface(bc=cb.ImpingementThermalBC.ConstantWallTemperature))
    assert heat_flux.Nu != wall_temp.Nu


def test_single_jet_wall_coupling_derivative_matches_finite_difference() -> None:
    surface = _single_jet_surface()
    v0 = 5.0
    dv = 1e-4

    r0 = _run(surface, velocity=v0)
    r_plus = _run(surface, velocity=v0 + dv)
    r_minus = _run(surface, velocity=v0 - dv)

    dh_dv_fd = (r_plus.h - r_minus.h) / (2.0 * dv)
    assert r0.dh_dmdot == pytest.approx(dh_dv_fd / _channel_rho_a(), rel=2e-3)


def test_the_element_flow_is_the_jet_flow_and_the_area_counts_holes() -> None:
    """#460. The element's mass flow is the jet flow; the convective area only
    counts holes. So at a fixed element flow, a larger target area means more
    holes, less flow each, and Re_j ~ 1/area. Before #460 the jet flow was
    rho v * area, which made Re_j independent of the area instead.
    """
    a = _array_surface()
    b = _array_surface()
    b.area = 2.0 * a.area
    assert _run(b).Re == pytest.approx(0.5 * _run(a).Re, rel=1e-12)
    assert _run(b).h < _run(a).h


def test_single_jet_does_not_depend_on_the_target_patch_area() -> None:
    """#460. One jet carries the element's whole flow; the target patch the
    reported h applies to does not set it. Before #460 the jet's flow was
    rho v * patch area."""
    a = _single_jet_surface()
    b = _single_jet_surface()
    b.area = 9.0 * a.area
    assert _run(a).Re == pytest.approx(_run(b).Re, rel=1e-12)
    assert _run(a).h == pytest.approx(_run(b).h, rel=1e-12)


def test_jet_reynolds_number_is_the_paper_definition() -> None:
    """Re_j = 4 (m_dot/n) / (pi d mu) with m_dot the element's flow through
    the channel cross-section (Florschuetz's G_j on the hole area)."""
    s = _array_surface()
    r = _run(s)
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    mu = cb.complete_state(STATE["T"], STATE["P"], STATE["X"]).transport.mu
    mdot = rho * FLOW["velocity"] * math.pi / 4.0 * FLOW["diameter"] ** 2
    n = s.area / (8.0 * 6.0 * 0.002**2)
    assert r.Re == pytest.approx(4.0 * (mdot / n) / (math.pi * 0.002 * mu), rel=1e-9)


def test_disabled_surface_returns_none() -> None:
    surface = ConvectiveSurface(area=0.0, model=ImpingementModel(d_jet=0.002))
    assert surface.htc_and_T(**STATE, **FLOW) is None
