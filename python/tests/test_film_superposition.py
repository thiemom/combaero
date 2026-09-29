"""Multi-row film superposition: Sellers, and Gao et al. (2025) Eq. (7).

The algebra is pinned in tests/test_film_superposition.cpp. These cover the
Python surface and the one thing the C++ tests cannot express as neatly: that
the pieces compose into a plate-scale prediction.
"""

import math

import pytest

import combaero as cb


def test_sellers_is_one_minus_the_product() -> None:
    etas = [0.30, 0.22, 0.18, 0.12]
    expected = 1.0 - math.prod(1.0 - e for e in etas)
    assert cb.film_superposition_sellers(etas) == pytest.approx(expected, rel=1e-14)


def test_corrected_reduces_to_sellers_at_alpha_one() -> None:
    etas = [0.30, 0.22, 0.18, 0.12]
    assert cb.film_superposition_corrected(etas, [1.0] * 3) == pytest.approx(
        cb.film_superposition_sellers(etas), rel=1e-14
    )


def test_alpha_below_one_damps_the_accumulation() -> None:
    """Gao's central finding: Sellers overestimates, worse with more rows."""
    etas = [0.30, 0.22, 0.18, 0.12]
    plain = cb.film_superposition_sellers(etas)
    for alpha in (0.95, 0.85, 0.70):
        assert cb.film_superposition_corrected(etas, [alpha] * 3) < plain


def test_the_overestimate_grows_with_row_count() -> None:
    """Why a few-row film model cannot simply be reused for effusion."""
    gaps = []
    for n in (2, 4, 8, 16, 30):
        etas = [0.08] * n
        gaps.append(
            cb.film_superposition_sellers(etas)
            - cb.film_superposition_corrected(etas, [0.9] * (n - 1))
        )
    assert all(b > a for a, b in zip(gaps, gaps[1:], strict=False))


def test_gradients_match_finite_differences() -> None:
    etas = [0.30, 0.22, 0.18, 0.12]
    alphas = [0.95, 0.92, 0.90]
    _, grad = cb.film_superposition_corrected_and_gradient(etas, alphas)
    h = 1e-7
    for i in range(len(etas)):
        up, dn = list(etas), list(etas)
        up[i] += h
        dn[i] -= h
        fd = (
            cb.film_superposition_corrected(up, alphas)
            - cb.film_superposition_corrected(dn, alphas)
        ) / (2 * h)
        assert grad[i] == pytest.approx(fd, abs=1e-6), f"row {i}"


def test_correction_factor_form() -> None:
    """Eq. (5): alpha = a r/(a r + 1) + b. a and b are REQUIRED.

    The paper publishes the form but never its fitted coefficients, so there
    is no default to lean on -- which is the point of making them positional.
    """
    a, b = 8.0, 0.1
    assert cb.mainstream_temperature_correction(0.0, a, b) == pytest.approx(b)
    assert cb.mainstream_temperature_correction(1e9, a, b) == pytest.approx(1.0 + b, rel=1e-6)
    with pytest.raises(TypeError):
        cb.mainstream_temperature_correction(0.05)  # no default for a, b


def test_equivalent_slot_and_blowing_ratio() -> None:
    d, pitch = 1.0e-3, 3.0e-3
    area = math.pi * d * d / 4
    assert cb.equivalent_slot_width(area, pitch) == pytest.approx(area / pitch)
    assert cb.equivalent_blowing_ratio(1.5, 2.0e-6, 4.0e-6) == pytest.approx(0.75)


def test_composes_with_baldauf_into_a_plate_prediction() -> None:
    """The point of all of it: one row's correlation, many rows' answer.

    Andrei's plate geometry -- 18 rows at s_x/d = 9.15, d = 1.5 mm, 30 deg --
    evaluated row by row with Baldauf and superposed. The result must be a
    physical effectiveness that exceeds any single row's contribution.
    """
    s_over_d, rows, alpha_deg, blowing, density = 9.15, 18, 30.0, 2.0, 1.2
    x_over_d_at_trailing_edge = rows * s_over_d

    per_row = []
    for i in range(rows):
        # Distance from row i to the plate's trailing edge, in diameters.
        x = (rows - i) * s_over_d
        per_row.append(
            cb.film_effectiveness_baldauf_2002(x, blowing, density, alpha_deg, 3.0, 0.015)
        )
    total = cb.film_superposition_sellers(per_row)

    assert 0.0 < total < 1.0
    assert total > max(per_row), "superposition must exceed any single row"
    assert x_over_d_at_trailing_edge > 0

    # And the corrected form pulls it back down, as Gao requires.
    damped = cb.film_superposition_corrected(per_row, [0.9] * (rows - 1))
    assert damped < total


@pytest.mark.parametrize(
    "etas,alphas",
    [
        ([], None),
        ([0.3, 1.2], None),
        ([0.3, 0.2], []),
        ([0.3, 0.2], [0.9, 0.9]),
        ([0.3, 0.2], [1.4]),
    ],
)
def test_refuses_malformed_input(etas, alphas) -> None:
    with pytest.raises(ValueError):
        if alphas is None:
            cb.film_superposition_sellers(etas)
        else:
            cb.film_superposition_corrected(etas, alphas)
