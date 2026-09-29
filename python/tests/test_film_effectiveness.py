"""Baldauf et al. (2002) film-cooling effectiveness, Python surface.

The physics is pinned in tests/test_film_effectiveness.cpp against the
paper's Fig. 14, Fig. 8 and Table 4. These tests cover the binding: that the
Python API exposes the same numbers, the same refusals, and the derivative
tuple in the documented order.
"""

import pytest

import combaero as cb

# Table 3, the paper's worked example and Fig. 14's case.
T3 = {"M": 2.0, "P": 1.2, "alpha_deg": 30.0, "s_over_D": 3.0, "Tu": 0.015}


def test_reproduces_figure_14_plateau() -> None:
    """Fig. 14 plots this exact case; the correlation plateaus at 0.11-0.12."""
    for x in (40.0, 60.0, 80.0):
        eta = cb.film_effectiveness_baldauf_2002(x_over_D=x, **T3)
        assert 0.105 < eta < 0.130, f"x/D = {x}"
    # Rises from the ejection point rather than peaking there, which is the
    # failure mode of the older exponential-decay correlations.
    assert cb.film_effectiveness_baldauf_2002(x_over_D=1.0, **T3) < 0.01


def test_derivative_tuple_order_and_agreement() -> None:
    eta, d_dM, d_dP = cb.film_effectiveness_baldauf_2002_and_derivatives(
        20.0, 2.0, 1.2, 30.0, 3.0, 0.015
    )
    assert eta == pytest.approx(cb.film_effectiveness_baldauf_2002(x_over_D=20.0, **T3), rel=1e-12)
    h = 1e-6
    fd_M = (
        cb.film_effectiveness_baldauf_2002(20.0, 2.0 + h, 1.2, 30.0, 3.0, 0.015)
        - cb.film_effectiveness_baldauf_2002(20.0, 2.0 - h, 1.2, 30.0, 3.0, 0.015)
    ) / (2 * h)
    fd_P = (
        cb.film_effectiveness_baldauf_2002(20.0, 2.0, 1.2 + h, 30.0, 3.0, 0.015)
        - cb.film_effectiveness_baldauf_2002(20.0, 2.0, 1.2 - h, 30.0, 3.0, 0.015)
    ) / (2 * h)
    assert d_dM == pytest.approx(fd_M, rel=1e-5)
    assert d_dP == pytest.approx(fd_P, rel=1e-5)


def test_derivative_changes_sign_through_lift_off() -> None:
    """Below the lift-off peak more blowing helps; above it, it hurts."""
    _, lo, _ = cb.film_effectiveness_baldauf_2002_and_derivatives(20.0, 0.4, 1.2, 30.0, 3.0, 0.015)
    _, hi, _ = cb.film_effectiveness_baldauf_2002_and_derivatives(20.0, 2.0, 1.2, 30.0, 3.0, 0.015)
    assert lo > 0.0
    assert hi < 0.0


def test_composes_with_the_adiabatic_wall_temperature_helper() -> None:
    """Baldauf's Eq. (1) is combaero's own eta convention, so no conversion.

    eta = (T_G - T_AW)/(T_G - T_C). Round-tripping through
    adiabatic_wall_temperature must recover it.
    """
    eta = cb.film_effectiveness_baldauf_2002(x_over_D=20.0, **T3)
    T_hot, T_cool = 1600.0, 800.0
    T_aw = cb.adiabatic_wall_temperature(T_hot, T_cool, eta)
    assert (T_hot - T_aw) / (T_hot - T_cool) == pytest.approx(eta, rel=1e-12)
    assert T_cool < T_aw < T_hot


def test_angle_is_in_degrees_at_the_api_boundary() -> None:
    """The API takes DEGREES; the paper's trigonometry is in radians.

    The unit itself is pinned in C++ against seven of the paper's Table 4
    coefficients. Here the concern is only that the boundary did not change
    unit: alpha is accepted over (0, 90] -- a radian-valued API would reject
    60 and 90 as out of range -- and shallower ejection lays a better film,
    which is the physical ordering the paper reports.
    """
    for alpha in (30.0, 45.0, 60.0, 90.0):
        assert cb.film_effectiveness_baldauf_2002(20.0, 1.0, 1.2, alpha, 3.0, 0.015) > 0.0
    shallow = cb.film_effectiveness_baldauf_2002(20.0, 1.0, 1.2, 30.0, 3.0, 0.015)
    normal = cb.film_effectiveness_baldauf_2002(20.0, 1.0, 1.2, 90.0, 3.0, 0.015)
    assert shallow > normal


@pytest.mark.parametrize(
    "kwargs",
    [
        {"x_over_D": 0.0, "M": 2.0, "P": 1.2, "alpha_deg": 30.0, "s_over_D": 3.0, "Tu": 0.015},
        {"x_over_D": 10.0, "M": 0.0, "P": 1.2, "alpha_deg": 30.0, "s_over_D": 3.0, "Tu": 0.015},
        {"x_over_D": 10.0, "M": 2.0, "P": 1.2, "alpha_deg": 30.0, "s_over_D": 3.0, "Tu": 0.0},
        {"x_over_D": 10.0, "M": 2.0, "P": 1.2, "alpha_deg": 0.0, "s_over_D": 3.0, "Tu": 0.015},
        {"x_over_D": 10.0, "M": 2.0, "P": 1.2, "alpha_deg": 120.0, "s_over_D": 3.0, "Tu": 0.015},
    ],
)
def test_refuses_inputs_the_correlation_cannot_represent(kwargs) -> None:
    with pytest.raises(ValueError):
        cb.film_effectiveness_baldauf_2002(**kwargs)
