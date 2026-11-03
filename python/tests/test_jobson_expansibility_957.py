"""jobson_expansibility converges at high beta Ci (#957).

The approach velocity raises the stagnation state the jet expands from, so
the expansibility solves g = F(g) for the flux per bore area g. F falls with
g (the root is unique), but at high beta Ci its slope at the root is below
-1: the fixed point the function used cycled to its cap and threw on valid
input (beta 0.5-0.9, Ci 0.6-1, bands of r). It is now a bracketed Newton on
F(g) - g. References are the same equation solved by bisection in 40-digit
arithmetic (mpmath), independent of the C++ code.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb

REFERENCE = [
    # r, beta, Ci, n, eps (40-digit bisection)
    (0.6, 0.8, 1.0, 1.4, 0.58374822830105710),
    (0.7, 0.75, 0.9, 1.4, 0.74605852644370772),
    (0.9, 0.8, 0.75, 1.4, 0.94770239954065802),
    (0.5, 0.9, 0.8, 1.667, 0.72983987911184448),
    (0.3, 0.8, 1.0, 1.4, 0.46813846063356227),  # converged before as well
]


@pytest.mark.parametrize(("r", "beta", "Ci", "n", "eps"), REFERENCE)
def test_expansibility_matches_the_40_digit_root(
    r: float, beta: float, Ci: float, n: float, eps: float
) -> None:
    assert cb.jobson_expansibility(r, beta, Ci, n) == pytest.approx(eps, rel=1e-13)


def test_the_formerly_failing_band_converges_and_is_continuous() -> None:
    """beta 0.8, Ci 1, n 1.4: on main r 0.52-0.72 threw. The branch returns
    a continuous curve through it, rising with r. (Its slope drops at r ~
    0.72, where the momentum quadratic's discriminant reaches zero and the
    closure holds it there: the flux kernel's `degenerate` branch, unchanged.)
    """
    rs = np.linspace(0.05, 0.999, 400)
    eps = np.array([cb.jobson_expansibility(float(r), 0.8, 1.0, 1.4) for r in rs])
    assert np.all(np.isfinite(eps))
    assert np.all(np.diff(eps) > 0.0)
    # No jump: no step exceeds three times the median step.
    assert np.max(np.diff(eps)) < 3.0 * np.median(np.diff(eps))
