"""The "huang1999" anchor closures of the ejector (P_c* by Huang's printed
chain, P_b0 by the sudden expansion): C++ against the Python reference, their
analytic Jacobians against central differences, their smoothness where a branch
switches, and the whole element's Jacobian with the closure selected.

Accuracy and fidelity of the closures themselves are scored in
validation/ejector (anchor_scorecard.py, printed_chain_check.py) and pinned in
test_ejector_anchor_scorecard.py / test_ejector_chen2013.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import combaero as cb
from combaero.network._ejector_huang1999 import (
    ETA_P,
    ETA_S,
    EjectorGeometry,
    huang_critical_back_pressure,
    huang_mixing_coefficient,
    sudden_expansion_dead_head,
)

st = cb._solver_tools
DATA = Path(__file__).resolve().parents[2] / "validation" / "ejector" / "data"
HUANG_CASES = json.loads((DATA / "huang1999_reference.json").read_text())["cases"]

# Air cases spanning the rigs scored cross-source (A_3/A_t 4.5-33).
AIR = [
    (pp, 293.15, 101325.0, 293.15, 1.376, a3, 1.4, 287.05)
    for pp in (3e5, 5e5, 2e6)
    for a3 in (4.46, 5.30, 33.3)
]
R141B = [
    (
        c["inputs"]["p_g_pa"],
        c["inputs"]["t_g_k"],
        c["inputs"]["p_e_pa"],
        c["inputs"]["t_e_k"],
        c["inputs"]["area_ratio_nozzle"],
        c["inputs"]["area_ratio_mix"],
        c["inputs"]["gamma"],
        c["inputs"]["r_gas"],
    )
    for c in HUANG_CASES
]
CASES = AIR + R141B


def _py(fn, case):
    p_g, t_g, p_e, t_e, arn, arm, g, r = case
    return fn(p_g, t_g, p_e, t_e, EjectorGeometry(arn, arm), g, r)


@pytest.mark.parametrize("area_ratio_mix", [6.0, 6.9, 7.5, 8.3, 9.0])
def test_mixing_coefficient_steps_at_huangs_area_ratios(area_ratio_mix: float) -> None:
    assert st.ejector_huang_mixing_coefficient(area_ratio_mix) == huang_mixing_coefficient(
        area_ratio_mix
    )


def test_cpp_matches_the_python_reference() -> None:
    for case in CASES:
        pc = st.ejector_huang_critical_back_pressure_and_jacobian(*case, 1.0).value
        ref = _py(huang_critical_back_pressure, case)
        assert pc.p_c_pa == pytest.approx(ref.p_c_pa, rel=1e-12)
        assert pc.mach_mixed == pytest.approx(ref.mach_mixed, rel=1e-12)
        pb = st.ejector_sudden_expansion_dead_head_and_jacobian(*case, 1.0).value
        ref_b = _py(sudden_expansion_dead_head, case)
        assert pb.p_c_pa == pytest.approx(ref_b.p_b0_pa, rel=1e-12)
        assert pb.mach_mixed == pytest.approx(ref_b.mach_mixed, rel=1e-12)


def _cd(f, x0: list[float], i: int, eps_rel: float = 1e-6) -> float:
    eps = eps_rel * abs(x0[i])
    xp, xm = list(x0), list(x0)
    xp[i] += eps
    xm[i] -= eps
    return (f(xp) - f(xm)) / (2.0 * eps)


@pytest.mark.parametrize(
    "binding",
    [
        "ejector_huang_critical_back_pressure_and_jacobian",
        "ejector_sudden_expansion_dead_head_and_jacobian",
    ],
)
def test_closure_jacobians_match_central_differences(binding: str) -> None:
    fn = getattr(st, binding)
    for case in CASES:
        x0 = list(case[:4])
        rest = case[4:]

        def value(x, fn=fn, rest=rest):
            return fn(*x, *rest, 1.0).value.p_c_pa

        res = fn(*x0, *rest, 1.0)
        analytic = [res.dpc_dp_g, res.dpc_dt_g, res.dpc_dp_e, res.dpc_dt_e]
        scale = res.value.p_c_pa
        for i in range(4):
            # J_ij |x_j| against the value: the residual-row-scaling metric
            fd = _cd(value, x0, i)
            assert abs(analytic[i] - fd) * abs(x0[i]) < 1e-6 * scale, (binding, case, i)


def test_huangs_chain_is_smooth_through_a_sonic_mixed_flow() -> None:
    """The shock branch switches on at M_m = 1. A normal shock's stagnation
    loss is O((M - 1)^3), so P_c* and its slope are continuous there. The
    post-shock Mach M_3 peaks at exactly 1 at the switch (M_3 = M_m below it,
    the shock's subsonic conjugate above), which locates it."""
    arn, arm, g, r = 1.376, 4.46, 1.4, 287.05

    def f(p_g: float):
        return st.ejector_huang_critical_back_pressure_and_jacobian(
            p_g, 293.15, 101325.0, 293.15, arn, arm, g, r, 1.0
        )

    lo, hi = 1.3e5, 1.2e6
    for _ in range(200):  # golden-section search for the M_3 peak
        a, b = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        lo, hi = (a, hi) if f(a).value.mach_mixed < f(b).value.mach_mixed else (lo, b)
    p_switch = 0.5 * (lo + hi)
    assert f(p_switch).value.mach_mixed == pytest.approx(1.0, abs=1e-6)
    assert 1.3e5 * 1.01 < p_switch < 1.2e6 * 0.99  # both branches are in range
    below, above = f(p_switch * (1 - 1e-6)), f(p_switch * (1 + 1e-6))
    assert below.value.p_c_pa == pytest.approx(above.value.p_c_pa, rel=1e-5)
    assert below.dpc_dp_g == pytest.approx(above.dpc_dp_g, rel=1e-3)


def test_the_dead_head_reduces_to_a_normal_shock_when_the_tube_is_the_nozzle_exit() -> None:
    """With A_3 = A_p1 there is no annulus: the sudden expansion is a normal
    shock at the nozzle-exit Mach M_p1, and P_b0 is the stagnation pressure
    behind it (normal-shock tables: p02/p01 at M_p1)."""
    g, arn = 1.4, 1.86
    res = st.ejector_sudden_expansion_dead_head_and_jacobian(
        5e5, 293.15, 101325.0, 293.15, arn, arn, g, 287.05, 1.0
    )
    m1 = st.ejector_entrainment_ratio_and_jacobian(
        5e5, 293.15, 1e5, 293.15, arn, 5.3, g
    ).value.mach_nozzle_exit
    m1s = m1 * m1
    p02_p01 = ((g + 1) * m1s / ((g - 1) * m1s + 2)) ** (g / (g - 1)) * (
        (g + 1) / (2 * g * m1s - (g - 1))
    ) ** (1 / (g - 1))
    assert res.value.p_c_pa == pytest.approx(5e5 * p02_p01, rel=1e-9)
    m2s = (1 + 0.5 * (g - 1) * m1s) / (g * m1s - 0.5 * (g - 1))
    assert res.value.mach_mixed == pytest.approx(m2s**0.5, rel=1e-9)


# Henry-like air ejector (d_t 3.3, d_1 4.5, d_3 7.6 mm) for the whole element.
A_T = 3.14159265 * 0.25 * 3.3e-3**2
A_E = 3.14159265 * 0.25 * 4.5e-3**2
A_3 = 3.14159265 * 0.25 * 7.6e-3**2


@pytest.mark.parametrize(
    ("p_out", "regime"),
    [(1.05e5, "plateau"), (1.45e5, "droop"), (2.4e5, "pinned above P_b0")],
)
def test_element_jacobian_with_huang_anchors(p_out: float, regime: str) -> None:
    g, r = 1.4, 287.05
    p_g, t = 5.0e5, 293.15
    cap = st.ejector_choked_mass_flow_and_jacobian(p_g, t, A_T, g, r, ETA_P)[0]
    x0 = [cap * 0.995, 0.6 * cap, 1.6 * cap, p_g, t, 1.01e5, t, p_out, 0.5e5]
    geo = (
        A_E / A_T,
        A_3 / A_T,
        A_T,
        A_E,
        A_3 - A_E,
        g,
        r,
        ETA_P,
        ETA_S,
        1.0,
        1e-3,
        0.95,
        0.999,
        0.98,
        1.02,
    )

    def rows(x):
        return st.ejector_element_residuals_and_jacobian(*x, *geo, anchor_closure="huang1999")

    res = rows(x0)
    for i in range(4):
        scale = max(abs(res.residuals[i]), max(abs(res.jacobian[i][j] * x0[j]) for j in range(9)))
        for j in range(9):
            fd = _cd(lambda x, i=i: rows(x).residuals[i], x0, j)
            assert abs(res.jacobian[i][j] - fd) * abs(x0[j]) < 1e-5 * scale, (regime, i, j)


def test_element_rejects_an_unknown_anchor_closure() -> None:
    geo = (
        1.86,
        5.3,
        A_T,
        A_E,
        A_3 - A_E,
        1.4,
        287.05,
        ETA_P,
        ETA_S,
        1.0,
        1e-3,
        0.95,
        0.999,
        0.98,
        1.02,
    )
    with pytest.raises(ValueError, match="anchor_closure"):
        st.ejector_element_residuals_and_jacobian(
            1e-3, 5e-4, 1.5e-3, 5e5, 293.0, 1e5, 293.0, 1.2e5, 5e4, *geo, anchor_closure="lossless"
        )


def test_the_production_chain_is_the_one_the_fidelity_gap_was_measured_on() -> None:
    """printed_chain_check.printed_chain (where the fidelity gap and the
    hypotheses live) and the production closure are the same chain."""
    from validation.ejector.printed_chain_check import printed_chain

    for case in CASES:
        p_g, t_g, p_e, t_e, arn, arm, g, r = case
        geom = EjectorGeometry(arn, arm)
        assert huang_critical_back_pressure(p_g, t_g, p_e, t_e, geom, g, r).p_c_pa == pytest.approx(
            printed_chain(p_g, t_g, p_e, t_e, geom, g, r), rel=1e-13
        )
