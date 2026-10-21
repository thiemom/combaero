"""Cross-source accuracy on Eames et al. 2013's R245fa chiller
(validation/ejector/eames2013_check.py).

A REPORT, not a target: these tests pin what it says, so a closure change
that moves it is seen and the README's statement is updated with it.
"""

from __future__ import annotations

import pytest

from combaero.network._ejector_huang1999 import (
    EjectorGeometry,
    critical_back_pressure,
    entrainment_ratio,
    huang_critical_back_pressure,
)
from validation.ejector.data import eames2013_fig8 as eames
from validation.ejector.eames2013_check import design_point, report, score


@pytest.fixture(scope="module")
def rows():
    return score()


def _mean(v: list[float]) -> float:
    return sum(v) / len(v)


def test_every_fig8_point_is_scored(rows) -> None:
    assert len(rows) == 10
    assert {r.ejector for r in rows} == {"A", "B"}
    assert report(rows).startswith("CROSS-SOURCE ACCURACY")


def test_table1_is_one_nozzle_two_ducts() -> None:
    """Both area ratios imply the same 10.2 mm primary throat, which is what
    makes Table 1's ratio A_3/A_t."""
    throats = [g["d_mix_mm"] / g["area_ratio_mix"] ** 0.5 for g in eames.GEOMETRY.values()]
    assert throats[0] == pytest.approx(throats[1], rel=0.005)
    assert throats[0] == pytest.approx(10.2, abs=0.05)


def test_the_cop_conversion_reproduces_the_authors_design_pair() -> None:
    """Their code's COP 0.31 maps to omega 0.45 through the baked enthalpies,
    against the 0.44 they print: the conversion, not the kernel, is checked."""
    assert eames.omega_from_cop("design", eames.DESIGN["cop"]) == pytest.approx(
        eames.DESIGN["omega"], rel=0.03
    )


def test_the_missing_nozzle_exit_area_does_not_matter() -> None:
    """Table 1 has no nozzle exit diameter; the kernel's omega and both P_c*
    closures are independent of A_p1/A_t."""
    pid, _ej, t_g, t_e, _tc, _cop = eames.FIG8[0]
    p_g, p_e, _pc, _a, _b, gamma = eames.BAKED[pid]
    out = []
    for nar in (2.0, 4.0, 8.0):
        geom = EjectorGeometry(area_ratio_nozzle=nar, area_ratio_mix=8.6)
        args = (p_g, t_g + 273.15, p_e, t_e + 273.15, geom, gamma)
        out.append(
            (
                entrainment_ratio(*args).omega,
                huang_critical_back_pressure(*args, eames.R_SPEC).p_c_pa,
                critical_back_pressure(*args, eames.R_SPEC).p_c_pa,
            )
        )
    for o in out[1:]:
        assert o == pytest.approx(out[0], rel=1e-12)


def test_the_two_closures_bracket_every_measured_critical_pressure(rows) -> None:
    """Huang's printed chain (the default) is 5-25% low on every point and
    Kracik & Dvorak 3-35% high: the measurement lies between them."""
    for r in rows:
        assert -0.25 < r.pc_huang_err < -0.05, r.point
        assert 0.03 < r.pc_kd_err < 0.35, r.point
    assert abs(_mean([r.pc_huang_err for r in rows])) < abs(_mean([r.pc_kd_err for r in rows]))


def test_the_larger_duct_moves_both_closures_down(rows) -> None:
    """Same nozzle, larger duct (A -> B): both closures' P_c* error drops by
    10 points or more, and the omega overprediction grows. They share only
    Huang's Eqs. 1-8 front end, so the geometry effect sits there."""

    def by(ej: str, attr: str) -> float:
        return _mean([getattr(r, attr) for r in rows if r.ejector == ej])

    assert by("A", "pc_huang_err") - by("B", "pc_huang_err") > 0.10
    assert by("A", "pc_kd_err") - by("B", "pc_kd_err") > 0.10
    assert by("B", "omega_ratio") > by("A", "omega_ratio") + 0.4


def test_omega_is_above_the_cop_derived_value_on_every_point(rows) -> None:
    """1.4-3.3x the omega the plotted COP implies. The authors doubt that COP
    and their own code is about 2x their measurement too, so this is the
    uncertainty they state, not a band."""
    for r in rows:
        assert 1.3 < r.omega_ratio < 3.4, r.point


def test_the_authors_code_is_not_the_huang_front_end() -> None:
    """Model vs model at their design point: the kernel gives omega 0.11 on
    ejector A, the duct their code designed, against the code's 0.44."""
    dp = design_point()
    assert dp["A"][0] == pytest.approx(0.111, abs=0.005)
    assert dp["A"][0] < eames.DESIGN["omega"] / 3.0
