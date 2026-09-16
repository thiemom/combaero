"""Ribbed channels in ConvectiveSurface: wall weighting and the documented gap.

The correlation gives the RIBBED SIDE. Combining it with the smooth walls is
the element's job, because how many walls are ribbed is a design choice rather
than a property of the correlation.

Two asymmetries are pinned here because both follow from the source rather than
from convenience:

- friction needs no wall weighting -- the `f` in the roughness function's
  definition is already the four-sided channel value
- heat transfer does
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network.components import ConvectiveSurface, RibbedModel, SmoothModel

X_AIR = cb.species.dry_air()
STATE = {"T": 600.0, "P": 2.0e5, "X": X_AIR}
FLOW = {"velocity": 60.0, "diameter": 0.025, "length": 0.6, "T_hot": 900.0}


def _surface(**kw) -> ConvectiveSurface:
    return ConvectiveSurface(area=0.1, model=RibbedModel(e_D=0.06, p_e=10.0, **kw))


def _run(surface: ConvectiveSurface):
    return surface.htc_and_T(**STATE, **FLOW)


def test_ribbed_side_and_smooth_side_are_both_reported() -> None:
    """The channel average hides a 20% modelling choice; both halves are visible."""
    r = _run(_surface())
    assert r.h_ribbed > r.h_smooth > 0.0
    assert r.h_smooth < r.h < r.h_ribbed


@pytest.mark.parametrize(("n", "expect"), [(1, "least"), (2, "middle"), (4, "all ribbed")])
def test_wall_count_moves_the_average_monotonically(n: int, expect: str) -> None:
    values = {k: _run(_surface(n_ribbed_walls=k)).h for k in (1, 2, 4)}
    assert values[1] < values[2] < values[4]


def test_four_ribbed_walls_is_exactly_the_ribbed_side() -> None:
    """With no smooth wall left, the weighting must be an identity."""
    r = _run(_surface(n_ribbed_walls=4))
    assert r.h == pytest.approx(r.h_ribbed, rel=1e-12)


def test_an_unrepresentable_wall_count_is_rejected() -> None:
    """Three ribbed walls has no unambiguous geometry, so it is not guessed."""
    with pytest.raises(ValueError, match="1, 2 or 4"):
        _run(_surface(n_ribbed_walls=3))


def test_the_documented_gap_and_the_knob_that_closes_it() -> None:
    """Plain smooth walls under-predict Han's own channel average by ~20%.

    Han reports the average directly, implying h_s/h_r ~ 0.70, because ribs
    enhance the adjacent smooth wall by 10-50% as well -- which no correlation
    here covers. Using the plain smooth correlation gives ~0.42.

    This is the regression test for that gap: it pins the size of the
    under-prediction, and that `smooth_wall_Nu_multiplier` is the supported way
    to close it rather than a constant invented inside the library.
    """
    plain = _run(_surface())
    assert plain.h_smooth / plain.h_ribbed == pytest.approx(0.42, abs=0.02)

    tuned = _run(_surface(smooth_wall_Nu_multiplier=1.67))
    assert tuned.h_smooth / tuned.h_ribbed == pytest.approx(0.70, abs=0.02)
    assert tuned.h / plain.h == pytest.approx(1.20, abs=0.02)


def test_friction_comes_from_the_correlation_not_from_a_multiplier() -> None:
    """Correlations own their f. A ribbed channel's f is the four-sided value
    the roughness function is defined against, not pipe friction scaled by
    something -- the round-trip #331 removed moved a drop by -24.7%."""
    r = _run(_surface())
    direct = cb.evaluate_rib(
        cb.han_1988_orthogonal(),
        cb.RibGeometry(e_D=0.06, p_e=10.0, W_H=1.0, alpha_deg=90.0),
        r.Re,
    )
    assert r.f == pytest.approx(direct.f, rel=1e-12)


def test_friction_does_not_depend_on_the_wall_count() -> None:
    """Only heat transfer is weighted: f is already a channel-level quantity."""
    f_values = {n: _run(_surface(n_ribbed_walls=n)).f for n in (1, 2, 4)}
    assert f_values[1] == pytest.approx(f_values[2], rel=1e-12)
    assert f_values[2] == pytest.approx(f_values[4], rel=1e-12)


def test_the_user_knobs_still_apply_on_top() -> None:
    base = _surface()
    scaled = _surface()
    scaled.Nu_multiplier = 2.0
    scaled.f_multiplier = 3.0
    r0, r1 = _run(base), _run(scaled)
    assert r1.h == pytest.approx(2.0 * r0.h, rel=1e-12)
    assert r1.f == pytest.approx(3.0 * r0.f, rel=1e-12)


def test_a_user_supplied_set_is_used_and_G_is_inverse() -> None:
    """A user set replaces the shipped one, and G runs the other way.

    G is the heat-transfer ROUGHNESS function: a larger G means a lower Stanton
    number and so LESS heat transfer. Pinned because the sign is easy to get
    backwards -- this test was written the wrong way round first, and the
    implementation was right.
    """
    worse = cb.han_1988_orthogonal()
    worse.name, worse.source = "rig_3", "measured 2026-09"
    worse.provenance = cb.RibProvenance.User
    worse.C_G = 4.1  # above Han's 3.7

    better = cb.han_1988_orthogonal()
    better.C_G = 3.3  # below it

    default = _run(_surface()).h_ribbed
    assert _run(_surface(correlation_set=worse)).h_ribbed < default
    assert _run(_surface(correlation_set=better)).h_ribbed > default


def test_out_of_range_reports_but_still_returns() -> None:
    """Validity is advisory: the band belongs to the source's rig."""
    r = _run(_surface())
    assert isinstance(r.extrapolated, bool)
    assert r.h > 0.0


def test_smooth_model_is_unaffected() -> None:
    s = ConvectiveSurface(area=0.1, model=SmoothModel())
    assert s.htc_and_T(**STATE, **FLOW).h > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# ChannelElement residual path
# ---------------------------------------------------------------------------

from combaero.network.components import ChannelElement, NetworkMixtureState  # noqa: E402

Y_AIR = cb.mole_to_mass(X_AIR)


def _element(model) -> ChannelElement:
    return ChannelElement(
        id="ch",
        from_node="A",
        to_node="B",
        length=0.6,
        diameter=0.025,
        roughness=1e-5,
        surface=ConvectiveSurface(area=0.1, model=model),
    )


def _element_states(m_dot: float):
    a = NetworkMixtureState(T=600.0, P=2e5, Pt=2.1e5, Tt=610.0, m_dot=m_dot, Y=Y_AIR)
    b = NetworkMixtureState(T=595.0, P=1.95e5, Pt=1.95e5, Tt=605.0, m_dot=m_dot, Y=Y_AIR)
    return a, b


def _rib_element() -> ChannelElement:
    return _element(RibbedModel(e_D=0.06, p_e=10.0, W_H=1.0, n_ribbed_walls=2))


@pytest.mark.parametrize("m_dot", [0.05, 0.35, 1.0])
def test_ribbed_jacobian_matches_central_differences(m_dot: float) -> None:
    el = _rib_element()
    _, jac = el.residuals(*_element_states(m_dot))
    h = 1e-7
    fd = (
        el.residuals(*_element_states(m_dot + h))[0][0]
        - el.residuals(*_element_states(m_dot - h))[0][0]
    ) / (2.0 * h)
    assert jac[0]["ch.m_dot"] == pytest.approx(fd, rel=1e-6)


def test_ribbed_drop_follows_the_flow_but_its_derivative_does_not() -> None:
    """The odd-quantity case: magnitude from the correlation, sign from the
    flow. dP is odd in m_dot, so d(dP)/d(m_dot) is EVEN and must not flip at
    zero. Getting that backwards puts a sign error in the Jacobian exactly at
    the crossing."""
    el = _rib_element()
    out = {}
    for m in (0.35, -0.35):
        a, b = _element_states(m)
        res, jac = el.residuals(a, b)
        out[m] = (a.Pt - b.Pt - res[0], jac[0]["ch.m_dot"])

    assert out[0.35][0] == pytest.approx(-out[-0.35][0], rel=1e-12)
    assert out[0.35][1] == pytest.approx(out[-0.35][1], rel=1e-12)


def test_ribbed_drop_exceeds_smooth_at_the_same_flow() -> None:
    a, b = _element_states(0.35)
    rib = a.Pt - b.Pt - _rib_element().residuals(a, b)[0][0]
    smooth = a.Pt - b.Pt - _element(SmoothModel()).residuals(a, b)[0][0]
    assert rib > 2.0 * smooth


def test_friction_multiplier_scales_the_ribbed_drop() -> None:
    plain = _rib_element()
    scaled = _rib_element()
    scaled.surface.f_multiplier = 2.5
    a, b = _element_states(0.35)
    dp0 = a.Pt - b.Pt - plain.residuals(a, b)[0][0]
    dp1 = a.Pt - b.Pt - scaled.residuals(a, b)[0][0]
    assert dp1 == pytest.approx(2.5 * dp0, rel=1e-12)
