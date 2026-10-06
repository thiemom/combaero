"""Pin-fin sets (#335) through the Python API.

The C++ tests (tests/test_pin_fin_correlation.cpp) carry the numerics:
printed-equation reproduction, the friction blend, basis conversion and
analytic derivatives. These check that the guarantees survive the binding and
that the user form behaves as documented.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb

GEOM = cb.PinFinGeometry(S_D=2.5, X_D=2.5, H_D=1.0, N_rows=10)


def test_the_default_sets_evaluate_in_their_box() -> None:
    nu = cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), GEOM, 1.0e4, 0.7)
    f = cb.evaluate_pin_fin_friction(cb.metzger_1982_staggered_friction(), GEOM, 1.0e4)
    x = math.sqrt(1.0e8 + 1.0)
    assert nu.Nu == pytest.approx(0.135 * x**0.69 * 2.5**-0.34, rel=1e-12)
    assert not nu.extrapolated
    assert 0.09 < f.f < 0.1  # 0.0941 at the branch meeting point
    assert not f.extrapolated


def test_a_user_copy_records_that_it_is_the_users() -> None:
    mine = cb.metzger_1986_staggered_nu()
    mine.name = "rig_7"
    mine.source = "measured 2026-10"
    mine.provenance = cb.RibProvenance.User
    mine.C = 0.150
    cb.validate_pin_fin_nu_set(mine)
    ref = cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), GEOM, 2.0e4, 0.7).Nu
    assert cb.evaluate_pin_fin_nu(mine, GEOM, 2.0e4, 0.7).Nu == pytest.approx(
        ref * 0.150 / 0.135, rel=1e-12
    )
    # The shipped factory is unaffected: sets are values, not shared state.
    assert cb.metzger_1986_staggered_nu().provenance == cb.RibProvenance.Extracted


def test_mistakes_raise_at_validation_not_mid_solve() -> None:
    bad = cb.metzger_1982_staggered_friction()
    bad.C1 = 0.0
    with pytest.raises(ValueError):
        cb.validate_pin_fin_friction_set(bad)
    with pytest.raises(ValueError):
        cb.validate_pin_fin_geometry(cb.PinFinGeometry(S_D=1.0))
    # The evaluator itself never throws, even on a bad set.
    assert math.isfinite(cb.evaluate_pin_fin_friction(bad, GEOM, 1.0e4).f)


def test_converters_match_the_unit_cell() -> None:
    g = cb.PinFinGeometry(S_D=4.0, X_D=2.0 * math.sqrt(3.0), H_D=2.0, N_rows=4)
    xs = g.X_D * g.S_D
    assert cb.pin_fin_dprime_over_D(g) == pytest.approx(
        2.0 * (4.0 * xs - math.pi) / (2.0 * xs + 1.5 * math.pi), rel=1e-14
    )
    assert cb.pin_fin_amin_over_afrontal(g) == pytest.approx(3.0 / 4.0, rel=1e-14)
    frac = cb.pin_fin_area_fractions(g)
    assert frac.pin_over_total == pytest.approx(
        frac.pin / (frac.pin + frac.endwall_exposed), rel=1e-14
    )


def test_the_chyu_ratio_carries_no_friction_ratio() -> None:
    mod = cb.chyu_1998_inline_over_staggered()
    assert mod.has_f is False
    inline = cb.PinFinGeometry(arrangement=cb.PinArrangement.Inline, N_rows=7)
    r = cb.evaluate_pin_fin_modifier(mod, inline, 1.0e4)
    assert r.ratio_f == 1.0
    assert 0.8 < r.ratio_Nu < 0.9  # 0.846 at Re 1e4


@pytest.mark.parametrize("arrangement", ["Inline", "Staggered"])
@pytest.mark.parametrize("fillet", [False, True])
def test_chyu_1990_friction_reproduces_its_digitised_points(arrangement: str, fillet: bool) -> None:
    """The Fitted sets must sit on the Fig. 6 points they were fitted to.

    Fit quality, not fidelity: this guards the transcription of the fit into
    the factory, and that the data on disk is what the fit came from.
    """
    from validation.cooling.schema import DATA_ROOT

    arr = getattr(cb.PinArrangement, arrangement)
    tag = f"{arrangement.lower()}_{'fillet' if fillet else 'straight'}"
    path = next((DATA_ROOT / "chyu1990").glob(f"fig6_{tag}_*.csv"))
    rows = [tuple(map(float, ln.split(","))) for ln in path.read_text().splitlines()[1:]]
    s = cb.chyu_1990_friction(arr, fillet)
    g = cb.PinFinGeometry(arrangement=arr, N_rows=7)
    for x, f_chyu in rows:
        f = cb.evaluate_pin_fin_friction(s, g, x * 1.0e4).f
        assert f * 4.0 == pytest.approx(f_chyu, rel=0.04), (tag, x)


def test_fin_efficiency_through_the_binding() -> None:
    e = cb.pin_fin_array_efficiency(h=2.0e3, k_pin=20.0, D=1e-3, H=2e-3, A_f_over_A_t=0.3)
    assert 0.5 < e.eta_fin < 0.9
    assert e.eta_t == pytest.approx(1.0 - 0.3 * (1.0 - e.eta_fin), rel=1e-14)
    assert e.deta_t_dh < 0.0
