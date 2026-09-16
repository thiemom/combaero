"""Parametrised rib correlations, from the Python side.

The C++ tests in tests/test_rib_correlation.cpp cover the chain, the guards and
the derivative. These cover what only matters once the type crosses pybind11:
that a caller can build their own set, that provenance stays visible, and that
malformed sets are rejected at the boundary rather than producing numbers.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb

_core = cb._core


def han() -> object:
    return _core.han_1988_orthogonal()


def ref_geometry() -> object:
    return _core.RibGeometry(e_D=0.047, p_e=10.0, W_H=1.0, alpha_deg=90.0)


def test_shipped_set_reproduces_the_confirmed_extraction() -> None:
    """Values from validation/cooling/extractions/han_ribbed.md, not from here."""
    r = _core.evaluate_rib(han(), ref_geometry(), 10000.0)
    assert pytest.approx(3.2000, abs=1e-4) == r.R
    assert r.f == pytest.approx(0.04576, abs=1e-5)
    assert r.e_plus == pytest.approx(71.1, abs=0.1)
    assert pytest.approx(12.21, abs=0.01) == r.G
    assert r.St_r == pytest.approx(0.00968, abs=1e-5)


def test_shipped_set_carries_its_provenance() -> None:
    """A built-in is worth shipping because its origin is visible, not because
    it is accurate: users are expected to replace the coefficients."""
    s = han()
    assert s.provenance == _core.RibProvenance.Extracted
    assert "Han" in s.source and "1988" in s.source
    assert s.accuracy_R == pytest.approx(0.06)
    assert s.accuracy_G == pytest.approx(0.08)
    assert s.valid_Pr == pytest.approx(0.7)


def test_a_user_set_is_structurally_distinguishable() -> None:
    """The point of the provenance class: a tuned coefficient cannot pass for a
    published one, because the field is required rather than conventional."""
    u = han()
    u.name = "my_rig"
    u.source = "measured on rig 3, 2026-09"
    u.provenance = _core.RibProvenance.User
    u.C_G = 4.1

    assert u.provenance == _core.RibProvenance.User
    assert han().provenance == _core.RibProvenance.Extracted
    moved = _core.evaluate_rib(u, ref_geometry(), 10000.0).G
    assert moved > _core.evaluate_rib(han(), ref_geometry(), 10000.0).G


def test_normaliser_is_data_not_convention() -> None:
    """The same correlation under a different reference needs a different
    constant. Mixing them is wrong by a constant factor -- 10^0.35 = 2.24x --
    which never looks like a trend and so survives review."""
    normalised = han()
    raw = han()
    raw.C_R = 3.2 * 10.0**-0.35
    raw.R_pe = _core.RibTerm(exponent=0.35, reference=1.0)

    g = ref_geometry()
    for p_e in (5.0, 10.0, 20.0):
        g.p_e = p_e
        assert (
            pytest.approx(_core.evaluate_rib(raw, g, 3e4).R, rel=1e-12)
            == _core.evaluate_rib(normalised, g, 3e4).R
        )

    mixed = han()
    mixed.R_pe = _core.RibTerm(exponent=0.35, reference=1.0)  # constant left at 3.2
    g.p_e = 10.0
    ratio = _core.evaluate_rib(mixed, g, 3e4).R / _core.evaluate_rib(normalised, g, 3e4).R
    assert ratio == pytest.approx(10.0**0.35, rel=1e-9)


@pytest.mark.parametrize("Re", [-1e6, -1.0, 0.0, 1e-12, 1e8])
def test_guards_survive_states_the_solver_probes(Re: float) -> None:
    r = _core.evaluate_rib(han(), ref_geometry(), Re)
    for field in (r.f, r.G, r.St_r, r.dSt_dRe):
        assert math.isfinite(field)
    assert r.f > 0.0
    assert r.St_r > 0.0


def test_reverse_flow_is_symmetric_in_magnitude() -> None:
    fwd = _core.evaluate_rib(han(), ref_geometry(), 3e4)
    rev = _core.evaluate_rib(han(), ref_geometry(), -3e4)
    assert pytest.approx(rev.G, rel=1e-12) == fwd.G
    assert fwd.St_r == pytest.approx(rev.St_r, rel=1e-12)
    assert fwd.e_plus == pytest.approx(-rev.e_plus, rel=1e-9)


def test_validity_warns_rather_than_refusing() -> None:
    """A band belongs to the source's rig, not the caller's hardware, so it
    reports rather than blocks. Users supplying their own set know their own."""
    g = ref_geometry()
    assert not _core.evaluate_rib(han(), g, 3e4).extrapolated

    g.e_D = 0.20  # far outside Han's 0.047-0.078
    r = _core.evaluate_rib(han(), g, 3e4)
    assert r.extrapolated
    assert math.isfinite(r.St_r) and r.St_r > 0.0


def test_malformed_sets_are_rejected_at_the_boundary() -> None:
    """Bad PARAMETERS are a mistake and fail loudly -- the opposite treatment
    from a bad operating point, which is guarded silently."""
    _core.validate_rib_set(han())

    no_source = han()
    no_source.source = ""
    with pytest.raises(Exception, match="source"):
        _core.validate_rib_set(no_source)

    zero_ref = han()
    zero_ref.R_pe = _core.RibTerm(exponent=0.35, reference=0.0)
    with pytest.raises(Exception, match="reference"):
        _core.validate_rib_set(zero_ref)

    nan_const = han()
    nan_const.C_G = float("nan")
    with pytest.raises(Exception, match="C_G"):
        _core.validate_rib_set(nan_const)


def test_friction_factor_does_not_depend_on_reynolds_number() -> None:
    """Structural to this family: R carries no e+ term. It is what makes
    df/d(mdot) identically zero for the element that will use this."""
    f_ref = _core.evaluate_rib(han(), ref_geometry(), 1e4).f
    for Re in (1.0, 1e3, 1e5, 1e7):
        assert _core.evaluate_rib(han(), ref_geometry(), Re).f == pytest.approx(f_ref, rel=1e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
