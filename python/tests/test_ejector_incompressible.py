"""The incompressible ejector: Sanger's jet pump, the ejector's low-Mach limit (#588).

Provenance chain (agreed on #588):
1. Sanger's theory against Sanger's data: NASA TN D-4759 (1968) Table II, his
   12 tabulated head-rise points at s/d_t = 0, scored with his own
   separately calibrated loss coefficients (Appendix B) -- not fitted.
2. The closure reproduces Sanger's Eq. 11 (TN D-4445 Appendix B) to machine
   precision: fidelity.
3. The compressible ejector reduces to it as Ma^2 at matched losses, equal
   and unequal stream temperatures: the incompressible regime is the
   compressible one's limit.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

from combaero import _solver_tools as st
from combaero.network._ejector_huang1999 import ETA_P, ETA_S
from combaero.network.ejector_element import EjectorElement

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# Sanger, NASA TN D-4759 (1968), Table II: water jet pumps, throat 1.35 in,
# 5.66 diameters long, fully inserted nozzle; secondary inlet 15 psia.
# (R, Q1 gal/min, Q2 gal/min, head rise H5 - H2 ft)
TABLE_II = [
    (0.066, 28, 33, 16.1),
    (0.066, 28, 99, 10.4),
    (0.066, 28, 144, 5.3),
    (0.066, 35, 33, 24.8),
    (0.066, 35, 89, 19.2),
    (0.066, 35, 140, 13.1),
    (0.197, 63, 31, 24.0),
    (0.197, 63, 88, 14.2),
    (0.197, 63, 135, 4.8),
    (0.197, 83, 42, 41.7),
    (0.197, 83, 91, 30.6),
    (0.197, 83, 133, 19.3),
]
# D-4759 Appendix B: component calibrations, independent of the overall data.
K_S = {0.066: 0.09, 0.197: 0.14}
K_T = 0.076
K_D = 0.079

D_T = 1.35 * 0.0254
A_T = math.pi / 4 * D_T**2
GPM = 6.30902e-5
G = 9.80665
FT = 0.3048
RHO = 998.0


def _eq11_head_ft(r: float, q1: float, q2: float, ks: float) -> float:
    """Sanger's Eq. 11, written out independently of the closure."""
    m = q2 / q1
    vn = q1 * GPM / (r * A_T)
    bracket = (
        2 * r
        + 2 * r**2 * m**2 / (1 - r)
        - (1 + K_T + K_D) * r**2 * (1 + m) ** 2
        - (1 + ks) * r**2 * m**2 / (1 - r) ** 2
    )
    return vn**2 / (2 * G) * bracket / FT


def _closure_head_ft(r: float, q1: float, q2: float, ks: float) -> float:
    v = 1.0 / RHO
    mp, ms = q1 * GPM * RHO, q2 * GPM * RHO
    ae, as_, a3 = r * A_T, (1 - r) * A_T, A_T
    pt_s = 1.03e5
    p3 = pt_s - (1 + ks) * 0.5 * v * ms * ms / (as_ * as_)  # R1 = 0
    res = st.ejector_jetpump_incompressible_residuals_and_jacobian(
        mp,
        ms,
        mp + ms,
        0.0,
        v,
        pt_s,
        v,
        0.0,
        p3,
        ae,
        as_,
        a3,
        k_nozzle=0.0,
        k_suction=ks,
        k_throat=K_T,
        k_diffuser=K_D,
        recovery=1.0,
    )
    pt_o = -res.residuals[3]  # R3 at Pt_o = 0 is -(the outlet Pt the closure sets)
    return (pt_o - pt_s) / (RHO * G) / FT


def test_the_closure_reproduces_sangers_eq_11() -> None:
    for r, q1, q2, _ in TABLE_II:
        assert _closure_head_ft(r, q1, q2, K_S[r]) == pytest.approx(
            _eq11_head_ft(r, q1, q2, K_S[r]), rel=1e-12
        )


def test_sangers_theory_scores_his_own_data_at_9_percent() -> None:
    """Accuracy of the THEORY on the author's own data, with his calibrated
    K's: a recorded number, not a test bound tuned to pass. Under 5% at low
    and mid flow ratio; it over-predicts toward the high end, where the head
    rise is small (+33% on 4.8 ft)."""
    errs = [(_closure_head_ft(r, q1, q2, K_S[r]) - h) / h for r, q1, q2, h in TABLE_II]
    mae = sum(abs(e) for e in errs) / len(errs)
    assert mae == pytest.approx(0.093, abs=0.001)


def test_lossless_is_worse_and_the_score_moves() -> None:
    """Falsifies the score: the same points without friction losses."""
    errs = []
    for r, q1, q2, h in TABLE_II:
        v = 1.0 / RHO
        mp, ms = q1 * GPM * RHO, q2 * GPM * RHO
        ae, as_, a3 = r * A_T, (1 - r) * A_T, A_T
        p3 = 1.03e5 - 0.5 * v * ms * ms / (as_ * as_)
        res = st.ejector_jetpump_incompressible_residuals_and_jacobian(
            mp,
            ms,
            mp + ms,
            0.0,
            v,
            1.03e5,
            v,
            0.0,
            p3,
            ae,
            as_,
            a3,
            k_nozzle=0.0,
            k_suction=0.0,
            k_throat=0.0,
            k_diffuser=0.0,
        )
        errs.append(((-res.residuals[3] - 1.03e5) / (RHO * G) / FT - h) / h)
    assert sum(abs(e) for e in errs) / len(errs) > 0.4


def test_the_default_losses_are_the_compressible_nozzle_efficiencies() -> None:
    ej = EjectorElement("e", "p", "s", "o", 3.14e-5, 1e-4, 8e-4, regime="incompressible")
    assert ej.k_nozzle == pytest.approx(1.0 / ETA_P - 1.0)
    assert ej.k_suction == pytest.approx(1.0 / ETA_S - 1.0)
    assert ej.k_throat == 0.0 and ej.k_diffuser == 0.0


def test_the_closure_jacobian_matches_central_differences() -> None:
    args = {
        "mp": 0.01,
        "ms": 0.04,
        "mo": 0.05,
        "pt_p": 1.03e5,
        "v_p": 0.8,
        "pt_s": 1.01e5,
        "v_s": 0.9,
        "pt_o": 1.015e5,
        "p3": 1.0e5,
        "area_nozzle_exit": 1e-4,
        "area_secondary": 7e-4,
        "area_mixing": 8e-4,
        "k_nozzle": 0.05,
        "k_suction": 0.17,
        "k_throat": 0.08,
        "k_diffuser": 0.1,
        "recovery": 0.97,
    }
    order = ["mp", "ms", "mo", "pt_p", "v_p", "pt_s", "v_s", "pt_o", "p3"]
    for sign in (1.0, -1.0):  # forward and reversed secondary
        a = {**args, "ms": sign * args["ms"]}
        base = st.ejector_jetpump_incompressible_residuals_and_jacobian(**a)
        for k, name in enumerate(order):
            h = 1e-6 * max(abs(a[name]), 1.0)
            up = st.ejector_jetpump_incompressible_residuals_and_jacobian(
                **{**a, name: a[name] + h}
            )
            dn = st.ejector_jetpump_incompressible_residuals_and_jacobian(
                **{**a, name: a[name] - h}
            )
            for row in range(4):
                fd = (up.residuals[row] - dn.residuals[row]) / (2 * h)
                assert base.jacobian[row][k] == pytest.approx(fd, rel=1e-6, abs=1e-6), (row, name)


# ---- R2: the compressible ejector's low-Mach limit -------------------------


def _jetpump_run(regime: str, m: float, t_p: float) -> tuple[float, float]:
    import test_gui_ejector as t

    from gui.backend.runner import NetworkRunner

    g = t._jetpump_ejector_schema()
    for n in g["nodes"]:
        if n["id"] == "mb_p":
            n["data"]["m_dot"] = m
            n["data"]["Tt"] = t_p
    g["solver_settings"] = {"global_regime": regime}
    res = NetworkRunner.from_dict(g).solve(timeout=120.0)
    raw = res._raw
    assert raw["__success__"], (regime, m, raw["__message__"])
    df = res.to_dataframe()
    omega = float(df[df["id"] == "ej1"].iloc[0]["omega"])
    return omega, raw["mc_p.Pt"] - 101325.0


@pytest.mark.parametrize("t_primary", [300.0, 450.0])
def test_is_the_low_mach_limit_of_the_compressible_ejector(t_primary: float) -> None:
    """Same network, both regimes, default losses: the gap in the primary's
    drive quarters as the jet's Mach halves -- equal and unequal stream
    temperatures (the volume-additive mixed density is the right limit)."""
    gaps = []
    for m in (0.0051, 0.00255, 0.001275):
        _, dp_i = _jetpump_run("incompressible", m, t_primary)
        _, dp_c = _jetpump_run("compressible", m, t_primary)
        gaps.append((dp_i - dp_c) / dp_c)
    assert abs(gaps[-1]) < 2e-3, gaps
    for hi, lo in zip(gaps, gaps[1:], strict=False):
        assert 3.7 < hi / lo < 4.3, gaps


def test_the_validity_report_judges_the_jet() -> None:
    """An incompressible ejector is judged at its primary jet's Mach (#588)."""
    import test_gui_ejector as t

    from gui.backend.runner import NetworkRunner

    g = t._jetpump_ejector_schema()
    for n in g["nodes"]:
        if n["id"] == "mb_p":
            n["data"]["m_dot"] = 0.01
    g["solver_settings"] = {"global_regime": "incompressible"}
    raw = NetworkRunner.from_dict(g).solve(timeout=120.0)._raw
    assert raw["__success__"], raw["__message__"]
    v = raw["__incompressible_validity__"]["ej1"]
    assert v["mach"] > 0.2 and v["verdict"] == "warn"
