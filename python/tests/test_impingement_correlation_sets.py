"""Jet impingement correlations, from the Python side.

The C++ tests in tests/test_impingement_correlation.cpp cover the formulas,
the guards and the physical invariants. These cover what only matters once
the type crosses pybind11: that a caller can read a set's provenance, that
malformed sets are rejected at the boundary rather than producing numbers,
and that the shipped sets still reproduce the confirmed extraction through
the Python API.
"""

from __future__ import annotations

import pytest

import combaero as cb

_core = cb._core


def test_single_jet_reproduces_the_confirmed_check_point() -> None:
    """Values from validation/cooling/extractions/han_impingement.md item 4,
    not from here. Han's book states these rounded to 60/56."""
    s = _core.goldstein_1986_single_jet()

    nu_q = _core.single_jet_impingement_nu(
        s, _core.ImpingementThermalBC.ConstantHeatFlux, 25000.0, 7.75, 5.0
    )
    assert nu_q == pytest.approx(60.0, abs=0.5)

    nu_t = _core.single_jet_impingement_nu(
        s, _core.ImpingementThermalBC.ConstantWallTemperature, 25000.0, 7.75, 5.0
    )
    assert nu_t == pytest.approx(56.0, abs=0.5)


def test_single_jet_set_carries_its_provenance() -> None:
    s = _core.goldstein_1986_single_jet()
    assert "Goldstein" in s.source
    assert pytest.approx(24.0) == s.A
    assert pytest.approx(533.0) == s.B
    assert pytest.approx(44.0) == s.C
    assert s.Re_exponent == pytest.approx(0.76)


def test_single_jet_validate_rejects_a_malformed_set() -> None:
    s = _core.goldstein_1986_single_jet()
    s.B = -1.0
    with pytest.raises(Exception):  # noqa: B017 - pybind translates std::invalid_argument
        _core.validate_single_jet_set(s)


def test_jet_array_table_coefficients_match_the_confirmed_extraction() -> None:
    """Table 4.1's coefficients, confirmed digit-for-digit against the primary
    paper -- see han_impingement.md items 17/18."""
    inline = _core.florschuetz_1981_inline()
    assert pytest.approx(1.18) == inline.A_fit.C
    assert inline.A_fit.nx == pytest.approx(-0.944)
    assert inline.n_fit.nz == pytest.approx(1.04)

    staggered = _core.florschuetz_1981_staggered()
    assert pytest.approx(1.87) == staggered.A_fit.C
    assert staggered.B_fit.nz == pytest.approx(0.059)


def test_inline_and_staggered_have_different_xn_d_validity() -> None:
    """A stated fact, not a smoothed-over detail: staggered's xn/d bound is
    genuinely tighter (5-10) than inline's (5-15)."""
    inline = _core.florschuetz_1981_inline()
    staggered = _core.florschuetz_1981_staggered()
    assert inline.valid_xn_d.hi == pytest.approx(15.0)
    assert staggered.valid_xn_d.hi == pytest.approx(10.0)


def test_jet_array_evaluates_at_a_real_tested_geometry() -> None:
    inline = _core.florschuetz_1981_inline()
    result = _core.jet_array_impingement_nu(
        inline, Re_j=10_000.0, Gc_Gj=0.3, Pr=0.7, xn_d=10.0, yn_d=6.0, z_d=2.0
    )
    assert result.Nu > 0.0
    assert result.extrapolated is False


def test_jet_array_validate_rejects_a_malformed_set() -> None:
    inline = _core.florschuetz_1981_inline()
    inline.m_fit.nx = float("nan")
    with pytest.raises(Exception):  # noqa: B017
        _core.validate_jet_array_set(inline)


def test_a_user_set_is_structurally_distinguishable() -> None:
    """Same design as the rib sets: a tuned coefficient is a different object,
    not a mutation that hides where it came from."""
    mine = _core.florschuetz_1981_inline()
    mine.name = "my_rig"
    mine.source = "measured on rig 3, 2026-09"
    mine.A_fit.C = 1.5

    stock = _core.florschuetz_1981_inline()
    assert stock.name == "florschuetz_1981_inline"
    assert pytest.approx(1.18) == stock.A_fit.C
    assert pytest.approx(1.5) == mine.A_fit.C


def test_crossflow_ratio_is_zero_at_first_row() -> None:
    """Nu1's own definition (Florschuetz nomenclature): Gc/Gj = 0 at row 1."""
    assert _core.crossflow_to_jet_ratio_at_row(
        yn_d=8.0, z_d=2.0, C_D=cb.FLORSCHUETZ_1981_DEFAULT_CD, row=1
    ) == pytest.approx(0.0)


def test_crossflow_ratio_depends_only_on_yn_d_times_z_d() -> None:
    """Stated explicitly by the source: the flow distribution depends only on
    (yn/d)(z/d), independent of streamwise hole spacing xn/d -- which this
    function does not even take as an argument."""
    a = _core.crossflow_to_jet_ratio_at_row(yn_d=8.0, z_d=2.0, C_D=0.79, row=6)
    b = _core.crossflow_to_jet_ratio_at_row(yn_d=4.0, z_d=4.0, C_D=0.79, row=6)
    assert a == pytest.approx(b, abs=1e-12)


def test_default_discharge_coefficient_matches_the_paper() -> None:
    assert pytest.approx(0.79) == cb.FLORSCHUETZ_1981_DEFAULT_CD
