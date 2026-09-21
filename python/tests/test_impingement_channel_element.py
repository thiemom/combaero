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
    return ConvectiveSurface(area=0.01, model=ImpingementModel(**model_kw))


def _single_jet_surface(**kw) -> ConvectiveSurface:
    model_kw = {"d_jet": 0.003, "L_D": 7.75, "R_D": 5.0}
    model_kw.update(kw)
    return ConvectiveSurface(area=0.001, model=SingleJetImpingementModel(**model_kw))


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

    # dh/dmdot = dh/dvelocity / (rho * area); recover dh/dvelocity by FD, then
    # convert both sides through the same rho*area factor used inside the
    # element so the comparison is apples to apples.
    dh_dv_fd = (r_plus.h - r_minus.h) / (2.0 * dv)
    rho_val = cb.density(STATE["T"], STATE["P"], STATE["X"])
    area = 0.01
    dh_dmdot_fd = dh_dv_fd / (rho_val * area)

    assert r0.dh_dmdot == pytest.approx(dh_dmdot_fd, rel=2e-3)


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
    rho_val = cb.density(STATE["T"], STATE["P"], STATE["X"])
    area = 0.001
    dh_dmdot_fd = dh_dv_fd / (rho_val * area)

    assert r0.dh_dmdot == pytest.approx(dh_dmdot_fd, rel=2e-3)


def test_disabled_surface_returns_none() -> None:
    surface = ConvectiveSurface(area=0.0, model=ImpingementModel(d_jet=0.002))
    assert surface.htc_and_T(**STATE, **FLOW) is None
