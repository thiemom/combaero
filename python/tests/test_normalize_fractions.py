"""normalize_fractions keeps absent species absent.

It used to floor every fraction at 1e-15, which put a trace of all 15
species into every mixed stream; the mixture sums then evaluated all of
them (dry air: 15 instead of 4), the Fanno channel kernel 3x slower for it.
"""

from __future__ import annotations

import pytest

import combaero as cb


def test_absent_species_stay_zero_and_the_rest_sums_to_one() -> None:
    X = list(cb.species.dry_air())
    out = cb.normalize_fractions([2.0 * x for x in X])
    assert [x == 0.0 for x in out] == [x == 0.0 for x in X]
    assert sum(out) == pytest.approx(1.0, rel=1e-15)


def test_negative_round_off_clamps_to_zero() -> None:
    X = list(cb.species.dry_air())
    X[5] = -1e-14
    out = cb.normalize_fractions(X)
    assert out[5] == 0.0
    assert min(out) >= 0.0


def test_a_mixed_stream_carries_no_invented_species() -> None:
    """The mixer path (normalize_fractions inside) used to return every
    species at >= 1e-15."""
    X = list(cb.species.dry_air())
    present = {i for i, x in enumerate(X) if x > 0.0}
    out = cb.normalize_fractions(cb.mole_to_mass(X))
    assert {i for i, y in enumerate(out) if y > 0.0} == present
