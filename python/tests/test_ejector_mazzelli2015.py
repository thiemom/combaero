"""Cross-source accuracy on the Florence R245fa ejector's third design
(validation/ejector/mazzelli2015_check.py, #697).

A REPORT, not a target: these tests pin what it says, so a closure change
that moves it is seen and the README's statement is updated with it.
"""

from __future__ import annotations

import pytest

from validation.ejector.data import mazzelli2015_tables as mm
from validation.ejector.mazzelli2015_check import (
    PREDICTED_PC_BIAS,
    report,
    score_2014,
    score_critical,
    score_plateau,
)


@pytest.fixture(scope="module")
def plateau():
    return score_plateau()


@pytest.fixture(scope="module")
def critical():
    return score_critical()


def test_the_third_duct_sits_on_the_same_nozzle_between_a_and_b() -> None:
    """31.8 mm on the 10.2 mm throat: A_3/A_t 9.72, between Eames 2013's
    ducts A (8.65) and B (13.2)."""
    assert mm.area_ratio_mix() == pytest.approx(9.72, abs=0.01)
    assert 8.65 < mm.area_ratio_mix() < 13.2


def test_the_plateau_and_the_knee_come_from_the_data_alone() -> None:
    """Within the measurement band of each sweep's largest omega: three
    points at T_e 5 C, two at 10 C; T_c* bracketed by the next point."""
    sw = mm.sweeps()
    assert [r[3] for r in mm.plateau(sw["E5"])] == [27.3, 27.8, 28.3]
    assert [r[3] for r in mm.plateau(sw["E10"])] == [28.7, 29.4]
    assert mm.critical_bracket(sw["E5"]) == (28.3, 28.8)
    assert mm.critical_bracket(sw["E10"]) == (29.4, 30.2)


def test_the_measured_omega_band_is_the_flow_bands_in_quadrature() -> None:
    w, u = mm.omega(mm.TABLES[0])
    assert w == pytest.approx(0.19 / 0.35)
    assert u == pytest.approx(w * ((0.03 / 0.35) ** 2 + (0.01 / 0.19) ** 2) ** 0.5)
    assert 0.09 < u / w < 0.11


def test_the_kernel_omega_is_low_on_every_plateau_point(plateau) -> None:
    """25-30 % below the measured m_s/m_p, outside its 95 % band on all
    five: the first omega on this rig scored against a measurement band.
    #687's COP-derived omega had the kernel about 2x HIGH on ducts A and B;
    against a trusted measurement it reads low."""
    assert len(plateau) == 5
    for r in plateau:
        assert -0.31 < r.omega_err < -0.24, r.point
        assert not r.omega_in_band, r.point


def test_the_primary_flow_is_a_sonic_throat(plateau) -> None:
    """The real-gas isentropic choke of saturated vapour at T_g lands inside
    the measured 0.35 +- 0.03 kg/s; the kernel's ideal-gas Eq. 1 at the
    omega-convention gamma is 4.5-9 % below the measurement, at the band's edge."""
    for r in plateau:
        assert abs(r.m_p_real - r.m_p_meas) <= r.m_p_band, r.point
        assert -0.09 < r.m_p_kernel / r.m_p_meas - 1.0 < -0.045, r.point


def test_the_two_closures_bracket_the_measured_critical_pressure(critical) -> None:
    """Huang's chain (the default) below the bracket, Kracik & Dvorak above
    it, both outside -- the pattern of Eames 2013 and the air rigs."""
    for c in critical:
        assert -0.20 < c.pc_huang_err < -0.12, c.sweep
        assert 0.08 < c.pc_kd_err < 0.15, c.sweep


def test_the_pre_registered_duct_trend_is_missed(critical) -> None:
    """#687's A -> B shift, interpolated to 9.72, predicted Huang -11.5 % and
    K&D +25.7 %. Both closures land 6-13 points lower: the third duct reads
    closer to B than its area ratio. Reported as found, not refitted."""
    h = sum(c.pc_huang / c.pc_mid - 1.0 for c in critical) / len(critical)
    k = sum(c.pc_kd / c.pc_mid - 1.0 for c in critical) / len(critical)
    assert h < PREDICTED_PC_BIAS["huang"] - 0.04
    assert k < PREDICTED_PC_BIAS["kd"] - 0.10
    assert h == pytest.approx(-0.176, abs=0.005)
    assert k == pytest.approx(0.130, abs=0.005)


def test_the_original_nozzle_point_reads_low_on_omega_too() -> None:
    """2014, exit 28.8 mm, T_g 100 C: omega 25 % below its COP's, like the
    replacement nozzle's measured flows; the kernel cannot tell the nozzles
    apart."""
    s = score_2014()
    assert s["omega_model"] / s["omega_from_cop"] == pytest.approx(0.75, abs=0.03)
    assert s["pc_huang_err"] < -0.25
    assert 0.0 < s["pc_kd_err"] < 0.05


def test_the_report_runs(plateau) -> None:
    assert report().startswith("CROSS-SOURCE ACCURACY")
