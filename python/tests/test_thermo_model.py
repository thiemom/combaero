"""The selectable species thermo (#549): NASA-9, or its C2 re-fit.

The numerical proof (continuity at every join, equivalence bounds, h' = cp and
T s' = cp in the bridges) is the C++ test tests/test_thermo_c2_refit.cpp; this
checks the Python switch and that it reaches the mixture functions.
"""

from __future__ import annotations

import pytest

import combaero as cb

AIR = cb.species.dry_air()


@pytest.fixture
def restore_model():
    """Start on the published NASA-9 and leave the model as it was found, so
    no test changes the thermo of the tests after it on the same worker."""
    before = cb.thermo_model()
    cb.set_thermo_model("nasa9")
    yield
    cb.set_thermo_model(before)


def test_the_default_is_the_c2_refit() -> None:
    """The solver's thermo since #549's step 3: C2 at every join, identical to
    NASA-9 outside the bridge windows, bounded tails past the data (#606)."""
    assert cb.thermo_model() == "nasa9_c2"


def test_the_switch_round_trips(restore_model) -> None:
    cb.set_thermo_model("nasa9_c2")
    assert cb.thermo_model() == "nasa9_c2"
    cb.set_thermo_model("nasa9")
    assert cb.thermo_model() == "nasa9"


def test_an_unknown_model_is_refused(restore_model) -> None:
    with pytest.raises(ValueError, match="nasa9_c2"):
        cb.set_thermo_model("nasa7")
    assert cb.thermo_model() == "nasa9"  # unchanged by the refused call


@pytest.mark.parametrize("T", [300.0, 900.0, 1100.0, 2500.0])
def test_identical_outside_the_join_windows(T: float, restore_model) -> None:
    ref = (cb.cp(T, AIR), cb.h(T, AIR), cb.s(T, AIR, 2e5))
    cb.set_thermo_model("nasa9_c2")
    assert (cb.cp(T, AIR), cb.h(T, AIR), cb.s(T, AIR, 2e5)) == ref


def test_inside_a_window_dcp_dt_no_longer_steps(restore_model) -> None:
    """Across 1000 K the mixture's dcp/dT steps on NASA-9 and not on the re-fit;
    cp itself moves by far less than the data's accuracy."""
    e = 1e-6
    step_nasa = abs(cb.dcp_dT(1000.0 + e, AIR) - cb.dcp_dT(1000.0 - e, AIR))
    cp_nasa = cb.cp(1000.0, AIR)
    cb.set_thermo_model("nasa9_c2")
    step_c2 = abs(cb.dcp_dT(1000.0 + e, AIR) - cb.dcp_dT(1000.0 - e, AIR))
    assert step_nasa > 1e-4 * abs(cb.dcp_dT(1000.0, AIR))
    assert step_c2 < 1e-8 * abs(cb.dcp_dT(1000.0, AIR))
    assert cb.cp(1000.0, AIR) == pytest.approx(cp_nasa, rel=1e-4)
