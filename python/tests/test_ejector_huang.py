"""Validate the ejector reference model against the golden fixture.

The single source of truth is validation/ejector/data/huang1999_reference.json,
generated from the Python reference model with gamma baked in from the CoolProp
EOS. The SAME fixture (and its C++ header sibling huang1999_reference_data.h) is
what the future C++ implementation, pybind layer, and GUI validate against, so
none of them ever needs a CoolProp environment.

Checks per case, for both the critical-mode entrainment ratio omega (Huang's
Eqs. 1-8) and the critical back pressure P_c* (Kracik & Dvorak's mixing
closure, their Eqs. 7-13 -- see huang1999.py's module docstring for why P_c*
no longer uses Huang's original Eqs. 9-18):
  * the live model still reproduces the stored reference values (regression)
  * the reference omega and P_c* match their respective published theory
    columns (physics; P_c* to a looser tolerance -- see huang1999_reference.json's
    "p_c_note" for the full comparison against alternatives)
  * the stored d(omega)/d(A_3/A_t) matches the closed-form analytic derivative
    (demonstrates the analytic-vs-golden-Jacobian pattern the C++ side will use)
  * P_c* decreases monotonically as A_3/A_t increases at fixed inlets (Fig.
    9/10 trend)

References: docs/ejector/Huang_1d_analysis_ejector.pdf,
docs/ejector/Development_of_an_Analytical_Method_for.pdf (Kracik & Dvorak)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import validation.ejector.data as _data_pkg
from validation.ejector.models.huang1999 import (
    ETA_P,
    ETA_S,
    CriticalPressureResult,
    EjectorGeometry,
    EjectorResult,
    critical_back_pressure,
    entrainment_ratio,
)

_FIXTURE_PATH = Path(_data_pkg.__file__).parent / "huang1999_reference.json"
_FIXTURE = json.loads(_FIXTURE_PATH.read_text())
_CASES = _FIXTURE["cases"]
_OMEGA_REL_TOL = _FIXTURE["omega_rel_tol"]
_PC_REL_TOL = _FIXTURE["pc_rel_tol"]

# Tight tolerance for reproducing the stored reference (same math, so this only
# catches an actual change in the model or a stale fixture).
_REGRESSION_RTOL = 1.0e-9


def _case_id(case: dict) -> str:
    return case["id"]


def _model_result(case: dict) -> EjectorResult:
    i = case["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    return entrainment_ratio(i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"])


def _model_omega(case: dict) -> float:
    return _model_result(case).omega


def _model_pc(case: dict) -> CriticalPressureResult:
    i = case["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    return critical_back_pressure(
        i["p_g_pa"],
        i["t_g_k"],
        i["p_e_pa"],
        i["t_e_k"],
        geom,
        i["gamma"],
        i["r_gas"],
        i["recovery_efficiency"],
    )


def test_fixture_is_populated() -> None:
    assert len(_CASES) == 31
    assert _FIXTURE["fluid"]["name"] == "R141b"


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_model_reproduces_reference(case: dict) -> None:
    """Live model still matches the checked-in golden omega."""
    assert _model_omega(case) == pytest.approx(case["reference"]["omega"], rel=_REGRESSION_RTOL)


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_reference_matches_paper_theory(case: dict) -> None:
    """Golden omega reproduces Huang's published theory column."""
    assert case["reference"]["omega"] == pytest.approx(
        case["paper_theory"]["omega"], rel=_OMEGA_REL_TOL
    )


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_stored_jacobian_matches_analytic(case: dict) -> None:
    """The stored FD d(omega)/d(A_3/A_t) matches the closed-form derivative.

    omega = (P_e/P_g)(A_3/A_t - A_py/A_t) sqrt(T_g/T_e) sqrt(eta_s/eta_p), and
    A_py/A_t is independent of A_3/A_t, so
        d(omega)/d(A_3/A_t) = (P_e/P_g) sqrt(T_g/T_e) sqrt(eta_s/eta_p).
    This is exactly the kind of analytic-vs-golden check the C++ analytic
    Jacobian will run against huang1999_reference_data.h.
    """
    i = case["inputs"]
    analytic = (
        (i["p_e_pa"] / i["p_g_pa"]) * math.sqrt(i["t_g_k"] / i["t_e_k"]) * math.sqrt(ETA_S / ETA_P)
    )
    stored = case["jacobian"]["domega_d_area_ratio_mix"]
    assert stored == pytest.approx(analytic, rel=1.0e-5)


def test_entrainment_decreases_with_back_pressure_matched_geometry() -> None:
    """Sanity: at fixed inlets a smaller A_3/A_t (higher critical back
    pressure) entrains less flow -- the monotonic trend of Fig. 9/10."""
    i = _CASES[0]["inputs"]
    geom_large = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    geom_small = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"] - 2.0)
    w_large = entrainment_ratio(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom_large, i["gamma"]
    ).omega
    w_small = entrainment_ratio(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom_small, i["gamma"]
    ).omega
    assert w_small < w_large


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_model_reproduces_reference_pc(case: dict) -> None:
    """Live model still matches the checked-in golden P_c* chain outputs."""
    pc = _model_pc(case)
    r = case["reference"]
    assert pc.p_c_pa == pytest.approx(r["p_c_pa"], rel=_REGRESSION_RTOL)
    assert pc.p_mixed_stagnation_pa == pytest.approx(
        r["p_mixed_stagnation_pa"], rel=_REGRESSION_RTOL
    )
    assert pc.temp_mixed_stagnation_k == pytest.approx(
        r["temp_mixed_stagnation_k"], rel=_REGRESSION_RTOL
    )
    assert pc.lambda_mixed == pytest.approx(r["lambda_mixed"], rel=_REGRESSION_RTOL)
    assert pc.mach_mixed == pytest.approx(r["mach_mixed"], rel=_REGRESSION_RTOL)


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_reference_matches_paper_theory_pc(case: dict) -> None:
    """Golden P_c* (recovery_efficiency=1.0) reproduces Huang's T_c*-derived
    theory column to PC_REL_TOL -- see huang1999_reference.json's "p_c_note"
    for the full comparison (Huang's original chain: mean 25%, max 35%;
    Kracik & Dvorak's closure used here: mean 6.2%, max 13.2%)."""
    assert case["reference"]["p_c_pa"] == pytest.approx(
        case["paper_theory"]["p_c_pa"], rel=_PC_REL_TOL
    )


def test_recovery_efficiency_scales_p_c_linearly() -> None:
    """recovery_efficiency=1.0 (the default) must reproduce the lossless
    mixed stagnation pressure p03 exactly, and P_c* must scale linearly with
    it -- it is a pure multiplier on p03, not folded into the mixing solve
    itself (lambda3, p03, T03 are all independent of recovery_efficiency)."""
    i = _CASES[0]["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    lossless = critical_back_pressure(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"], 1.0
    )
    assert lossless.p_c_pa == pytest.approx(lossless.p_mixed_stagnation_pa, rel=_REGRESSION_RTOL)

    derated = critical_back_pressure(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"], 0.9
    )
    assert derated.p_c_pa == pytest.approx(
        0.9 * lossless.p_mixed_stagnation_pa, rel=_REGRESSION_RTOL
    )
    assert derated.p_mixed_stagnation_pa == pytest.approx(
        lossless.p_mixed_stagnation_pa, rel=_REGRESSION_RTOL
    )
    assert derated.lambda_mixed == pytest.approx(lossless.lambda_mixed, rel=_REGRESSION_RTOL)


@pytest.mark.parametrize("case", _CASES, ids=_case_id)
def test_pc_is_physically_bounded(case: dict) -> None:
    """P_c* must sit strictly between the suction and primary pressures: the
    diffuser recovers some of the kinetic energy lost to mixing, but not
    enough to reach the primary supply pressure back."""
    i = case["inputs"]
    p_c = case["reference"]["p_c_pa"]
    assert i["p_e_pa"] < p_c < i["p_g_pa"]


def test_critical_back_pressure_decreases_with_area_ratio_mix() -> None:
    """Sanity: at fixed inlets, a larger A_3/A_t (Fig. 9/10's independent
    variable) trades higher omega for a LOWER achievable critical back
    pressure -- the two halves of the same critical-mode plateau. Holds for
    Kracik & Dvorak's closure too, even though A_3/A_t only enters indirectly
    (through omega feeding the mixed-flow lambda3/z3 solve, not directly)."""
    i = _CASES[0]["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    geom_larger = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"] + 2.0)
    pc = critical_back_pressure(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"]
    ).p_c_pa
    pc_larger = critical_back_pressure(
        i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom_larger, i["gamma"], i["r_gas"]
    ).p_c_pa
    assert pc_larger < pc


def test_from_count_one_matches_direct_construction() -> None:
    """A single nozzle is the trivial case: from_count(count=1, ...) must
    reduce exactly to the direct EjectorGeometry constructor."""
    geom = EjectorGeometry.from_count(
        count=1, area_throat_single=2.905e-6, area_exit_single=8.44e-6, area_mix=2.6e-5
    )
    direct = EjectorGeometry(area_ratio_nozzle=8.44e-6 / 2.905e-6, area_ratio_mix=2.6e-5 / 2.905e-6)
    assert geom.area_ratio_nozzle == pytest.approx(direct.area_ratio_nozzle)
    assert geom.area_ratio_mix == pytest.approx(direct.area_ratio_mix)


def test_from_count_nozzle_ratio_invariant_to_count() -> None:
    """A_p1/A_t is a per-nozzle property: every nozzle contributes
    proportionally to both areas, so it must not depend on how many
    identical nozzles are aggregated."""
    area_throat_single = 2.905e-6
    area_exit_single = 8.44e-6
    geoms = [
        EjectorGeometry.from_count(
            count=n,
            area_throat_single=area_throat_single,
            area_exit_single=area_exit_single,
            area_mix=2.6e-5,
        )
        for n in (1, 2, 4, 12)
    ]
    ratios = [g.area_ratio_nozzle for g in geoms]
    assert ratios == pytest.approx([ratios[0]] * len(ratios))


def test_from_count_mix_ratio_shrinks_with_count() -> None:
    """A_3/A_t must shrink as nozzle count grows: the aggregate throat area
    scales with count, but the shared mixing chamber does not (Kracik &
    Dvorak's peripheral multi-nozzle topology)."""
    area_mix_ratios = [
        EjectorGeometry.from_count(
            count=n, area_throat_single=2.905e-6, area_exit_single=8.44e-6, area_mix=2.6e-5
        ).area_ratio_mix
        for n in (1, 2, 4, 12)
    ]
    assert area_mix_ratios == sorted(area_mix_ratios, reverse=True)
    assert area_mix_ratios[0] == pytest.approx(12 * area_mix_ratios[-1])
