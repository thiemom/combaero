"""Generate the language-neutral Huang-1999 ejector validation fixture.

Emits two artifacts from the single Python reference model, so the future C++
implementation (analytic Jacobian, ctests), the pybind layer, and the GUI can
all validate against identical golden values WITHOUT recreating a CoolProp
environment -- gamma is already baked in per suction condition.

  huang1999_reference.json       canonical, human-readable, language-neutral
  huang1999_reference_data.h     constexpr array for zero-dependency ctests

Each case carries:
  * inputs           -- everything the model needs, incl. the baked gamma
  * reference        -- this model's omega + intermediates (regression golden;
                        the C++ forward solve must match these tightly)
  * jacobian         -- central-difference d(omega)/d(input) (the golden target
                        the C++ ANALYTIC Jacobian must match to FD tolerance)
  * paper_theory     -- Huang's published theory column (the physics check)

Regenerate (no CoolProp needed; gamma comes from the checked-in table):

    uv run python validation/ejector/data/generate_reference.py

gamma itself is regenerated separately and rarely, via generate_gamma.py, which
is the only step that needs CoolProp.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable
from pathlib import Path

# Allow running as a plain script (python validation/.../generate_reference.py):
# put the repo root on sys.path so the validation package imports resolve.
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from validation.ejector.data.huang1999_tables import ALL_ROWS, HuangRow  # noqa: E402
from validation.ejector.models.huang1999 import (  # noqa: E402
    ETA_P,
    ETA_S,
    PHI_P,
    EjectorGeometry,
    blended_entrainment_ratio,
    critical_back_pressure,
    dead_head_back_pressure,
    entrainment_ratio,
    subcritical_entrainment_ratio,
)

# R141b constants from the CoolProp EOS (baked so this generator needs no
# CoolProp): molar mass 116.9496 g/mol -> R_specific = 8.314462618 / M.
R141B_MOLAR_MASS_KG_PER_MOL = 0.1169496
R141B_R_SPECIFIC_J_PER_KGK = 71.094

# Relative tolerance on omega vs the paper theory column (physics check).
OMEGA_REL_TOL = 0.03

# Relative tolerance on P_c* (recovery_efficiency=1.0) vs the T_c*-derived
# theory column. Measured mean 6.2%, max 13.2% across all 31 rows (Kracik &
# Dvorak's mixing closure) -- comfortable margin above the observed max.
PC_REL_TOL = 0.15

_HERE = Path(__file__).resolve().parent
_JSON_PATH = _HERE / "huang1999_reference.json"
_HEADER_PATH = _HERE / "huang1999_reference_data.h"

# Continuous inputs the finite-difference Jacobian is taken over. The physical
# solver unknowns are the pressures/temperatures; geometry, gamma, r_gas and
# recovery_efficiency are parameters but are included so the C++ side can
# validate sensitivities too. (Some entries are naturally zero for omega,
# which doesn't depend on r_gas or recovery_efficiency -- kept for a uniform
# Jacobian shape between the two outputs.)
_JAC_INPUTS = (
    "p_g_pa",
    "t_g_k",
    "p_e_pa",
    "t_e_k",
    "area_ratio_nozzle",
    "area_ratio_mix",
    "gamma",
    "r_gas",
    "recovery_efficiency",
)


def _omega_from_vector(v: dict[str, float]) -> float:
    geom = EjectorGeometry(v["area_ratio_nozzle"], v["area_ratio_mix"])
    return entrainment_ratio(
        v["p_g_pa"], v["t_g_k"], v["p_e_pa"], v["t_e_k"], geom, v["gamma"]
    ).omega


def _pc_from_vector(v: dict[str, float]) -> float:
    geom = EjectorGeometry(v["area_ratio_nozzle"], v["area_ratio_mix"])
    return critical_back_pressure(
        v["p_g_pa"],
        v["t_g_k"],
        v["p_e_pa"],
        v["t_e_k"],
        geom,
        v["gamma"],
        v["r_gas"],
        v["recovery_efficiency"],
    ).p_c_pa


def _pb0_from_vector(v: dict[str, float]) -> float:
    geom = EjectorGeometry(v["area_ratio_nozzle"], v["area_ratio_mix"])
    return dead_head_back_pressure(
        v["p_g_pa"],
        v["t_g_k"],
        v["p_e_pa"],
        v["t_e_k"],
        geom,
        v["gamma"],
        v["r_gas"],
        v["recovery_efficiency"],
    ).p_b0_pa


def _central_diff(f: Callable[[dict[str, float]], float], v: dict[str, float], key: str) -> float:
    x = v[key]
    h = 1.0e-6 * abs(x) if x != 0.0 else 1.0e-6
    hi = dict(v)
    lo = dict(v)
    hi[key] = x + h
    lo[key] = x - h
    return (f(hi) - f(lo)) / (2.0 * h)


def _case(row: HuangRow) -> dict:
    inputs = {
        "p_g_pa": row.p_g_pa,
        "t_g_k": row.t_g_k,
        "p_e_pa": row.p_e_pa,
        "t_e_k": row.t_e_k,
        "area_ratio_nozzle": row.area_ratio_nozzle,
        "area_ratio_mix": row.a3_at_theory,
        "gamma": row.gamma,
        "r_gas": R141B_R_SPECIFIC_J_PER_KGK,
        "recovery_efficiency": 1.0,
    }
    geom = EjectorGeometry(inputs["area_ratio_nozzle"], inputs["area_ratio_mix"])
    res = entrainment_ratio(
        inputs["p_g_pa"], inputs["t_g_k"], inputs["p_e_pa"], inputs["t_e_k"], geom, inputs["gamma"]
    )
    pc_res = critical_back_pressure(
        inputs["p_g_pa"],
        inputs["t_g_k"],
        inputs["p_e_pa"],
        inputs["t_e_k"],
        geom,
        inputs["gamma"],
        inputs["r_gas"],
        inputs["recovery_efficiency"],
    )
    dh_res = dead_head_back_pressure(
        inputs["p_g_pa"],
        inputs["t_g_k"],
        inputs["p_e_pa"],
        inputs["t_e_k"],
        geom,
        inputs["gamma"],
        inputs["r_gas"],
        inputs["recovery_efficiency"],
    )
    omega_jac = {k: _central_diff(_omega_from_vector, inputs, k) for k in _JAC_INPUTS}
    pc_jac = {k: _central_diff(_pc_from_vector, inputs, k) for k in _JAC_INPUTS}
    pb0_jac = {k: _central_diff(_pb0_from_vector, inputs, k) for k in _JAC_INPUTS}
    return {
        "id": f"{row.table}-{row.ejector}-Pg{row.p_g_mpa:g}",
        "table": row.table,
        "ejector": row.ejector,
        "inputs": inputs,
        "reference": {
            "omega": res.omega,
            "mach_nozzle_exit": res.mach_nozzle_exit,
            "mach_hypothetical_throat": res.mach_hypothetical_throat,
            "area_ratio_primary_core": res.area_ratio_primary_core,
            "area_ratio_entrained": res.area_ratio_entrained,
            "p_c_pa": pc_res.p_c_pa,
            "p_mixed_stagnation_pa": pc_res.p_mixed_stagnation_pa,
            "temp_mixed_stagnation_k": pc_res.temp_mixed_stagnation_k,
            "lambda_mixed": pc_res.lambda_mixed,
            "mach_mixed": pc_res.mach_mixed,
            "p_b0_pa": dh_res.p_b0_pa,
            "lambda_dead_head": dh_res.lambda_mixed,
        },
        "jacobian": {
            **{f"domega_d_{k}": omega_jac[k] for k in _JAC_INPUTS},
            **{f"dp_c_d_{k}": pc_jac[k] for k in _JAC_INPUTS},
            **{f"dp_b0_d_{k}": pb0_jac[k] for k in _JAC_INPUTS},
        },
        "paper_theory": {
            "omega": row.omega_theory,
            "area_ratio_mix": row.a3_at_theory,
            "tc_star_c": row.tc_star_c,
            "p_c_pa": row.pc_star_theory_pa,
        },
    }


# Fractions of the [P_c*, P_b0] window swept for the subcritical droop case.
# A plateau point (f < 0, clipped to omega_crit), the corner (f = 0 = P_c*),
# and the droop toward dead-head (f -> 1). eps_frac is the smooth-min rounding.
_SUBCRITICAL_FRACS = (-0.10, 0.0, 0.05, 0.20, 0.40, 0.60, 0.80, 0.95)
_SUBCRITICAL_EPS_FRAC = 1.0e-3


def _subcritical_sweep(row: HuangRow) -> dict:
    """Back-pressure sweep of the blended entrainment omega_eff for one
    geometry: the golden target for the Phase 2/3 R1 residual and its new
    d(omega_eff)/d(P_b) coupling (OPERATING_REGIMES_DESIGN.md sec 3.2)."""
    inputs = {
        "p_g_pa": row.p_g_pa,
        "t_g_k": row.t_g_k,
        "p_e_pa": row.p_e_pa,
        "t_e_k": row.t_e_k,
        "area_ratio_nozzle": row.area_ratio_nozzle,
        "area_ratio_mix": row.a3_at_theory,
        "gamma": row.gamma,
        "r_gas": R141B_R_SPECIFIC_J_PER_KGK,
        "recovery_efficiency": 1.0,
    }
    geom = EjectorGeometry(inputs["area_ratio_nozzle"], inputs["area_ratio_mix"])
    omega_crit = entrainment_ratio(
        inputs["p_g_pa"], inputs["t_g_k"], inputs["p_e_pa"], inputs["t_e_k"], geom, inputs["gamma"]
    ).omega
    p_c = critical_back_pressure(
        inputs["p_g_pa"], inputs["t_g_k"], inputs["p_e_pa"], inputs["t_e_k"],
        geom, inputs["gamma"], inputs["r_gas"],
    ).p_c_pa
    p_b0 = dead_head_back_pressure(
        inputs["p_g_pa"], inputs["t_g_k"], inputs["p_e_pa"], inputs["t_e_k"],
        geom, inputs["gamma"], inputs["r_gas"],
    ).p_b0_pa

    sweep = []
    for frac in _SUBCRITICAL_FRACS:
        p_b = p_c + frac * (p_b0 - p_c)
        omega_eff = blended_entrainment_ratio(omega_crit, p_c, p_b0, p_b, _SUBCRITICAL_EPS_FRAC)
        h = 1.0e-6 * abs(p_b)
        hi = blended_entrainment_ratio(omega_crit, p_c, p_b0, p_b + h, _SUBCRITICAL_EPS_FRAC)
        lo = blended_entrainment_ratio(omega_crit, p_c, p_b0, p_b - h, _SUBCRITICAL_EPS_FRAC)
        sweep.append(
            {
                "frac": frac,
                "p_b_pa": p_b,
                "omega_eff": omega_eff,
                "omega_sub": subcritical_entrainment_ratio(omega_crit, p_c, p_b0, p_b),
                "domega_eff_d_p_b": (hi - lo) / (2.0 * h),
            }
        )
    return {
        "id": f"{row.table}-{row.ejector}-Pg{row.p_g_mpa:g}-subcritical",
        "inputs": inputs,
        "anchors": {
            "omega_crit": omega_crit,
            "p_c_pa": p_c,
            "p_b0_pa": p_b0,
            "eps_frac": _SUBCRITICAL_EPS_FRAC,
        },
        "sweep": sweep,
    }


def build_fixture() -> dict:
    return {
        "schema": "combaero.ejector.huang1999.v1",
        "reference": (
            "B.J. Huang, J.M. Chang, C.P. Wang, V.A. Petrenko, "
            "A 1-D analysis of ejector performance, "
            "Int. J. Refrigeration 22 (1999) 354-364"
        ),
        "model": (
            "critical-mode 1-D ejector: Huang's entrainment ratio (Eqs. 1-8) + "
            "Kracik & Dvorak's mixing closure for P_c* (their Eqs. 7-13)"
        ),
        "generated_by": "validation/ejector/data/generate_reference.py",
        "fluid": {
            "name": "R141b",
            "molar_mass_kg_per_mol": R141B_MOLAR_MASS_KG_PER_MOL,
            "R_specific_J_per_kgK": R141B_R_SPECIFIC_J_PER_KGK,
            "gamma_note": (
                "ideal-gas gamma at the entrained-flow choking plane from the "
                "CoolProp EOS; per suction condition, baked in generate_gamma.py"
            ),
        },
        "coefficients": {"eta_p": ETA_P, "eta_s": ETA_S, "phi_p": PHI_P},
        "omega_rel_tol": OMEGA_REL_TOL,
        "pc_rel_tol": PC_REL_TOL,
        "p_c_note": (
            "P_c* (recovery_efficiency=1.0) uses Kracik & Dvorak (2016)'s mixing "
            "closure -- mass, momentum and energy solved together directly from "
            "the y-y state to the fully-mixed subsonic state, no phi_m loss "
            "coefficient and no separate shock -- in place of Huang's original "
            "Eqs. 9-18 (phi_m-weighted momentum + Rankine-Hugoniot shock + "
            "isentropic diffuser), which was systematically low (mean 25%, max "
            "35%) against the T_c*-derived target across all 31 rows for reasons "
            "confirmed NOT to be gamma (monotonic in gamma over its whole valid "
            "range, gamma -> 1 limit still short) or an equation error (hand- "
            "verified; reproduces classical normal-shock tables exactly). The "
            "new closure cuts this to mean 6.2%, max 13.2% -- PC_REL_TOL=0.15 "
            "reflects that with margin. Independently, S. Akbarnejad and "
            "M. Ziabasharhagh's CFD (J. Braz. Soc. Mech. Sci. Eng. 47:253, 2025, "
            "DOI 10.1007/s40430-025-05536-7) measured a naive momentum-balance "
            "error of up to 30% at critical back pressure -- close to Huang's "
            "gap here -- but their own proposed fix (drop momentum, impose an "
            "exact-sonic condition) is a geometry-free DESIGN tool and does not "
            "transplant onto this given-geometry forward analysis (tested: mean "
            "32-60% depending on the ambiguous mixing-pressure convention used, "
            "worse than Huang). recovery_efficiency defaults to 1.0 (no "
            "artificial loss) rather than a value fitted to this dataset -- see "
            "README.md's Accuracy section for the full comparison."
        ),
        "jacobian_method": "central difference, h = 1e-6 * |x|",
        "p_b0_note": (
            "P_b0 (dead-head back pressure) is the zero-entrainment upper anchor "
            "of the subcritical droop -- Kracik & Dvorak's SAME mixing closure at "
            "omega = 0 (gam = 0), no new constant. Present on every case; the "
            "subcritical_cases section sweeps the back pressure between P_c* and "
            "P_b0 and records the blended entrainment omega_eff plus its "
            "d(omega_eff)/d(P_b) coupling -- the Phase 2/3 R1 golden target. See "
            "validation/ejector/OPERATING_REGIMES_DESIGN.md sec 2.2, 3.2."
        ),
        "cases": [_case(r) for r in ALL_ROWS],
        "subcritical_cases": [_subcritical_sweep(ALL_ROWS[0])],
    }


def _fmt(x: float) -> str:
    return repr(float(x))


def _emit_header(fixture: dict) -> str:
    cases = fixture["cases"]
    lines: list[str] = []
    lines.append("// AUTO-GENERATED by validation/ejector/data/generate_reference.py")
    lines.append("// Do not edit by hand. Golden validation data for the Huang-1999")
    lines.append("// critical-mode 1-D ejector model (see huang1999_reference.json).")
    lines.append("#pragma once")
    lines.append("")
    lines.append("#include <array>")
    lines.append("")
    lines.append("namespace combaero::validation::ejector {")
    lines.append("")
    lines.append("struct HuangCase {")
    lines.append("  const char *id;")
    lines.append("  // Inputs (gamma baked from CoolProp EOS; no runtime dependency).")
    lines.append("  double p_g_pa, t_g_k, p_e_pa, t_e_k;")
    lines.append("  double area_ratio_nozzle, area_ratio_mix, gamma, r_gas;")
    lines.append("  double recovery_efficiency;")
    lines.append("  // Reference-model forward outputs (C++ forward solve must match).")
    lines.append("  double omega, mach_nozzle_exit, mach_hypothetical_throat;")
    lines.append("  double area_ratio_primary_core, area_ratio_entrained;")
    lines.append("  double p_c_pa, p_mixed_stagnation_pa, temp_mixed_stagnation_k;")
    lines.append("  double lambda_mixed, mach_mixed;")
    lines.append("  // Dead-head (omega -> 0) upper anchor P_b0 and its mixed lambda3.")
    lines.append("  double p_b0_pa, lambda_dead_head;")
    lines.append("  // Central-difference d(omega)/d(input): analytic Jacobian target.")
    lines.append("  double domega_dp_g, domega_dt_g, domega_dp_e, domega_dt_e;")
    lines.append("  double domega_dar_nozzle, domega_dar_mix, domega_dgamma, domega_dr_gas;")
    lines.append("  double domega_drecovery_eff;")
    lines.append("  // Central-difference d(P_c*)/d(input): analytic Jacobian target.")
    lines.append("  double dpc_dp_g, dpc_dt_g, dpc_dp_e, dpc_dt_e;")
    lines.append("  double dpc_dar_nozzle, dpc_dar_mix, dpc_dgamma, dpc_dr_gas;")
    lines.append("  double dpc_drecovery_eff;")
    lines.append("  // Central-difference d(P_b0)/d(input): analytic Jacobian target.")
    lines.append("  double dpb0_dp_g, dpb0_dt_g, dpb0_dp_e, dpb0_dt_e;")
    lines.append("  double dpb0_dar_nozzle, dpb0_dar_mix, dpb0_dgamma, dpb0_dr_gas;")
    lines.append("  double dpb0_drecovery_eff;")
    lines.append("  // Huang published theory column (physics check).")
    lines.append("  double omega_theory;")
    lines.append("  // P_c* target from T_c* (theory), via CoolProp saturation pressure.")
    lines.append("  // See kHuangPcRelTol and README.md's Accuracy section.")
    lines.append("  double tc_star_c, p_c_theory_pa;")
    lines.append("};")
    lines.append("")
    lines.append(
        "inline constexpr double kHuangOmegaRelTol = " + _fmt(fixture["omega_rel_tol"]) + ";"
    )
    lines.append("inline constexpr double kHuangPcRelTol = " + _fmt(fixture["pc_rel_tol"]) + ";")
    lines.append("")
    lines.append(f"inline constexpr std::array<HuangCase, {len(cases)}> kHuangCases = {{{{")
    for c in cases:
        i = c["inputs"]
        r = c["reference"]
        j = c["jacobian"]
        vals = [
            f'"{c["id"]}"',
            _fmt(i["p_g_pa"]),
            _fmt(i["t_g_k"]),
            _fmt(i["p_e_pa"]),
            _fmt(i["t_e_k"]),
            _fmt(i["area_ratio_nozzle"]),
            _fmt(i["area_ratio_mix"]),
            _fmt(i["gamma"]),
            _fmt(i["r_gas"]),
            _fmt(i["recovery_efficiency"]),
            _fmt(r["omega"]),
            _fmt(r["mach_nozzle_exit"]),
            _fmt(r["mach_hypothetical_throat"]),
            _fmt(r["area_ratio_primary_core"]),
            _fmt(r["area_ratio_entrained"]),
            _fmt(r["p_c_pa"]),
            _fmt(r["p_mixed_stagnation_pa"]),
            _fmt(r["temp_mixed_stagnation_k"]),
            _fmt(r["lambda_mixed"]),
            _fmt(r["mach_mixed"]),
            _fmt(r["p_b0_pa"]),
            _fmt(r["lambda_dead_head"]),
            _fmt(j["domega_d_p_g_pa"]),
            _fmt(j["domega_d_t_g_k"]),
            _fmt(j["domega_d_p_e_pa"]),
            _fmt(j["domega_d_t_e_k"]),
            _fmt(j["domega_d_area_ratio_nozzle"]),
            _fmt(j["domega_d_area_ratio_mix"]),
            _fmt(j["domega_d_gamma"]),
            _fmt(j["domega_d_r_gas"]),
            _fmt(j["domega_d_recovery_efficiency"]),
            _fmt(j["dp_c_d_p_g_pa"]),
            _fmt(j["dp_c_d_t_g_k"]),
            _fmt(j["dp_c_d_p_e_pa"]),
            _fmt(j["dp_c_d_t_e_k"]),
            _fmt(j["dp_c_d_area_ratio_nozzle"]),
            _fmt(j["dp_c_d_area_ratio_mix"]),
            _fmt(j["dp_c_d_gamma"]),
            _fmt(j["dp_c_d_r_gas"]),
            _fmt(j["dp_c_d_recovery_efficiency"]),
            _fmt(j["dp_b0_d_p_g_pa"]),
            _fmt(j["dp_b0_d_t_g_k"]),
            _fmt(j["dp_b0_d_p_e_pa"]),
            _fmt(j["dp_b0_d_t_e_k"]),
            _fmt(j["dp_b0_d_area_ratio_nozzle"]),
            _fmt(j["dp_b0_d_area_ratio_mix"]),
            _fmt(j["dp_b0_d_gamma"]),
            _fmt(j["dp_b0_d_r_gas"]),
            _fmt(j["dp_b0_d_recovery_efficiency"]),
            _fmt(c["paper_theory"]["omega"]),
            _fmt(c["paper_theory"]["tc_star_c"]),
            _fmt(c["paper_theory"]["p_c_pa"]),
        ]
        lines.append("    {" + ", ".join(vals) + "},")
    lines.append("}};")
    lines.append("")
    lines.extend(_emit_subcritical(fixture))
    lines.append("}  // namespace combaero::validation::ejector")
    lines.append("")
    return "\n".join(lines)


def _emit_subcritical(fixture: dict) -> list[str]:
    """Subcritical back-pressure sweep: blended entrainment omega_eff(P_b) and
    its d(omega_eff)/d(P_b) coupling -- the Phase 2/3 R1 golden target."""
    sub_cases = fixture["subcritical_cases"]
    n_points = len(sub_cases[0]["sweep"])
    lines: list[str] = []
    lines.append("// Subcritical droop sweep (Phase 1A): omega_eff = smooth-min of the")
    lines.append("// critical plateau and the Tier-1 linear droop between P_c* and P_b0.")
    lines.append("struct SubcriticalPoint {")
    lines.append("  double frac, p_b_pa, omega_eff, omega_sub, domega_eff_dp_b;")
    lines.append("};")
    lines.append("struct SubcriticalCase {")
    lines.append("  const char *id;")
    lines.append("  double omega_crit, p_c_pa, p_b0_pa, eps_frac;")
    lines.append(f"  std::array<SubcriticalPoint, {n_points}> sweep;")
    lines.append("};")
    lines.append("")
    lines.append(
        f"inline constexpr std::array<SubcriticalCase, {len(sub_cases)}> "
        "kSubcriticalCases = {{"
    )
    for c in sub_cases:
        a = c["anchors"]
        head = [
            f'"{c["id"]}"',
            _fmt(a["omega_crit"]),
            _fmt(a["p_c_pa"]),
            _fmt(a["p_b0_pa"]),
            _fmt(a["eps_frac"]),
        ]
        pts = []
        for p in c["sweep"]:
            pts.append(
                "{"
                + ", ".join(
                    _fmt(p[k]) for k in ("frac", "p_b_pa", "omega_eff", "omega_sub", "domega_eff_d_p_b")
                )
                + "}"
            )
        lines.append("    {" + ", ".join(head) + ", {{" + ", ".join(pts) + "}}},")
    lines.append("}};")
    lines.append("")
    return lines


def main() -> None:
    fixture = build_fixture()
    _JSON_PATH.write_text(json.dumps(fixture, indent=2) + "\n")
    _HEADER_PATH.write_text(_emit_header(fixture))
    print(f"wrote {_JSON_PATH.relative_to(_HERE.parents[2])} ({len(fixture['cases'])} cases)")
    print(f"wrote {_HEADER_PATH.relative_to(_HERE.parents[2])}")


if __name__ == "__main__":
    main()
