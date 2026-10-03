"""Ratio-form rib sets (#444) through the Python API.

The C++ tests (tests/test_rib_ratio_correlation.cpp) carry the numerics:
exact in-range reproduction, the C1 handover, symmetry and the analytic
derivatives. These check that a user can build, validate and evaluate a set
from Python, and that the guarantees survive the binding.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb


def _set() -> cb.RibRatioSet:
    s = cb.RibRatioSet()
    s.name = "user_rig"
    s.source = "measured, our rig"
    s.C_Nu = 2.5
    s.Nu_Re = cb.RibTerm(-0.2, 1.0e4)
    s.C_f = 8.0
    s.f_Re = cb.RibTerm(0.1, 1.0e4)
    s.Nu0_source = cb.RatioBaseline(0.023, 0.8, 0.4)
    s.f0_source = cb.RatioBaseline(0.046, -0.2)
    s.Re_floor = 1.0e4
    return s


GEOM = cb.RibGeometry(e_D=0.0625, p_e=10.0, W_H=1.0, alpha_deg=90.0)


def test_in_range_is_ratio_times_source_baseline() -> None:
    s = _set()
    cb.validate_rib_ratio_set(s)
    r = cb.evaluate_rib_ratio(s, GEOM, 3.0e4, 0.7)
    expect = 2.5 * (3.0e4 / 1e4) ** -0.2 * 0.023 * 3.0e4**0.8 * 0.7**0.4
    assert r.Nu == pytest.approx(expect, rel=1e-6)
    assert not r.below_floor


def test_below_floor_is_laminar_limited_not_zero() -> None:
    s = _set()
    zero = cb.evaluate_rib_ratio(s, GEOM, 0.0, 0.7)
    assert zero.below_floor and zero.extrapolated
    assert zero.Nu > 3.0, "Gnielinski handover keeps a laminar floor"

    held = cb.RibRatioOptions()
    held.below_floor = cb.RatioBelowFloor.SourceBaseline
    assert cb.evaluate_rib_ratio(s, GEOM, 0.0, 0.7, held).Nu < 0.01 * zero.Nu


@pytest.mark.parametrize("re", [-5.0e4, -3.0, 0.0, 3.0, 5.0e4])
def test_every_state_is_finite_and_even(re: float) -> None:
    s = _set()
    a = cb.evaluate_rib_ratio(s, GEOM, re, 0.7)
    b = cb.evaluate_rib_ratio(s, GEOM, -re, 0.7)
    for v in (a.Nu, a.f, a.dNu_dRe, a.df_dRe):
        assert math.isfinite(v)
    assert a.Nu >= 0.0 and a.f >= 0.0
    assert a.Nu == b.Nu and a.dNu_dRe == -b.dNu_dRe


def test_mistakes_are_rejected() -> None:
    s = _set()
    s.Nu0_source = cb.RatioBaseline()
    with pytest.raises(ValueError, match="Nu0_source"):
        cb.validate_rib_ratio_set(s)
    o = cb.RibRatioOptions()
    o.below_floor = cb.RatioBelowFloor.User
    with pytest.raises(ValueError, match="user_Nu0"):
        cb.validate_rib_ratio_options(o)
