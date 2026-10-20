"""Huang 1999's printed P_c* chain against its authors' OWN computed outputs.

Two authors evaluated the same chain (Huang Eqs. 9-18 with phi_m, Eq. 19):

* Huang 1999, Tables 3-4: the theory A_3/A_t that his code solved to hit each
  row's measured T_c*. Run forward at that area, the chain should return the
  T_c*-derived P_c* (``paper_theory.p_c_pa`` in the golden fixture).
* Chen 2013, Figs. 5a and 7a: his model curves (``chen2013_figs``, kind
  "model"). Their plateau end brackets his computed P_c*.

These are model-against-model comparisons (fidelity), so no rig bias enters:
the authors' codes took the same inputs we do. Printed as one table: the
printed chain, then single-change hypotheses for what their codes did instead.
None closes both authors at once (README "Why the printed chain misses its
authors' own numbers"). Also printed: Chen's subcritical branch as printed,
against his own droop.

Run: ``uv run python validation/ejector/printed_chain_check.py``
"""

from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from combaero.network._ejector_huang1999 import (
    ETA_P,
    ETA_S,
    PHI_P,
    EjectorGeometry,
    _area_over_astar,
    _mach_from_area_ratio_supersonic,
    entrainment_ratio,
)
from validation.ejector.data import chen2013_figs as chen

FIXTURE = Path(__file__).parent / "data" / "huang1999_reference.json"
AIR_GAMMA, AIR_R, CHEN_T_K, CHEN_PS_PA = 1.4, 287.05, 293.15, 1.0e5


def phi_m(area_ratio_mix: float) -> float:
    """Huang Eq. 19 (Chen Eq. 28)."""
    if area_ratio_mix > 8.3:
        return 0.80
    return 0.82 if area_ratio_mix >= 6.9 else 0.84


def printed_chain(
    p_g: float,
    t_g: float,
    p_e: float,
    t_e: float,
    geom: EjectorGeometry,
    g: float,
    r: float,
    hyp: str = "printed",
) -> float:
    """P_c* [Pa] from the critical y-y state; ``hyp`` swaps one step."""
    e = entrainment_ratio(p_g, t_g, p_e, t_e, geom, g)
    om, mpy, ppy = e.omega, e.mach_hypothetical_throat, e.p_mixing_pa
    tpy = t_g / (1 + 0.5 * (g - 1) * mpy**2)
    tsy = t_e * 2 / (g + 1)
    vpy = mpy * math.sqrt(g * r * tpy)
    vsy = math.sqrt(g * r * tsy)
    if hyp == "a_py at the stagnation T_g":
        vpy = mpy * math.sqrt(g * r * t_g)
    if hyp == "a_sy at the stagnation T_e":
        vsy = math.sqrt(g * r * t_e)
    if hyp == "V_py from the nozzle-exit state":
        mp1 = e.mach_nozzle_exit
        vpy = mp1 * math.sqrt(g * r * t_g / (1 + 0.5 * (g - 1) * mp1**2))
    psi = phi_m(geom.area_ratio_mix)
    if hyp == "psi as kinetic-energy efficiency":
        psi = math.sqrt(psi)
    if hyp == "psi on the primary only":
        vm = (psi * vpy + om * vsy) / (1 + om)
    elif hyp == "psi on the slip velocity only":
        vm = (vpy + om * vsy) / (1 + om) - (1 - psi) * (vpy - vsy) * om / (1 + om)
    else:
        vm = psi * (vpy + om * vsy) / (1 + om)
    cp = g * r / (g - 1)
    tm = (cp * (t_g + om * t_e) / (1 + om) - 0.5 * vm * vm) / cp
    if hyp == "mixed T from static enthalpies":
        tm = (tpy + om * tsy) / (1 + om)
    if hyp == "a_m at the mixed stagnation T":
        tm = (t_g + om * t_e) / (1 + om)
    mm = vm / math.sqrt(g * r * tm)
    if hyp == "shock at the primary-core Mach":
        mm = mpy
    if hyp == "no shock (isentropic from the mix)":
        return ppy * (1 + 0.5 * (g - 1) * mm**2) ** (g / (g - 1))
    if mm > 1:
        p3 = ppy * (1 + 2 * g / (g + 1) * (mm**2 - 1))
        m3sq = (1 + 0.5 * (g - 1) * mm**2) / (g * mm**2 - 0.5 * (g - 1))
    else:
        p3, m3sq = ppy, mm**2
    if hyp == "incompressible diffuser":
        return p3 + 0.5 * g * p3 * m3sq
    if hyp == "diffuser exponent inverted":
        return p3 * (1 + 0.5 * (g - 1) * m3sq) ** ((g - 1) / g)
    if hyp == "M3 read as M3^2 (Eq. 17)":
        return p3 * (1 + 0.5 * (g - 1) * m3sq**2) ** (g / (g - 1))
    return p3 * (1 + 0.5 * (g - 1) * m3sq) ** (g / (g - 1))


HYPOTHESES = [
    "printed",
    "psi as kinetic-energy efficiency",
    "psi on the primary only",
    "psi on the slip velocity only",
    "mixed T from static enthalpies",
    "shock at the primary-core Mach",
    "no shock (isentropic from the mix)",
    "incompressible diffuser",
    "diffuser exponent inverted",
    "M3 read as M3^2 (Eq. 17)",
    "a_py at the stagnation T_g",
    "a_sy at the stagnation T_e",
    "V_py from the nozzle-exit state",
    "a_m at the mixed stagnation T",
]


def chen_curves() -> list[tuple[str, float, EjectorGeometry, float, float, list[tuple]]]:
    """(fig, P_p, geometry, P_c* bracket lo/hi [bar], model points) per Chen curve."""
    out = []
    for fig, data in (("5a", chen.FIG5A), ("7a", chen.FIG7A)):
        a1, a2 = chen.AREA_RATIOS[fig]
        geom = EjectorGeometry(area_ratio_nozzle=a1, area_ratio_mix=a2)
        for pp in sorted({r[2] for r in data}):
            m = sorted((r for r in data if r[1] == "model" and r[2] == pp), key=lambda r: r[3])
            n = next(i for i, r in enumerate(m) if abs(r[4] - m[0][4]) > 0.01)
            out.append((fig, pp, geom, m[n - 1][3], m[n][3], m))
    return out


def chen_subcritical_omega(pp: float, geom: EjectorGeometry, pc: float) -> float:
    """Chen's subcritical branch as printed (his Eqs. 12-27): above P_c* the
    entrained flow unchokes, and P_sy rises until the chain's back pressure
    equals P_c. Air, P_p absolute, suction 1 bar, 20 C, as Chen ran it."""
    g, r, t0, ps = AIR_GAMMA, AIR_R, CHEN_T_K, CHEN_PS_PA
    gm1, cp = g - 1, g * r / (g - 1)
    mp1 = _mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, g)
    pp1 = pp / (1 + 0.5 * gm1 * mp1**2) ** (g / gm1)
    mp = pp / math.sqrt(t0) * math.sqrt(g / r * (2 / (g + 1)) ** ((g + 1) / gm1) * ETA_P)

    def at(p_sy: float) -> tuple[float, float]:
        core = (1 + 0.5 * gm1 * mp1**2) * (p_sy / pp1) ** (-gm1 / g)
        mpy = math.sqrt(max((core - 1) / (0.5 * gm1), 0.0))
        a_sy = geom.area_ratio_mix - PHI_P * _area_over_astar(mpy, g)
        tsy = t0 * (p_sy / ps) ** (gm1 / g)
        vsy = math.sqrt(max(2 * cp * (t0 - tsy), 0.0))
        om = p_sy / (r * tsy) * vsy * a_sy * math.sqrt(ETA_S) / mp
        tpy = t0 / (1 + 0.5 * gm1 * mpy**2)
        vpy = mpy * math.sqrt(g * r * tpy)
        vm = phi_m(geom.area_ratio_mix) * (vpy + om * vsy) / (1 + om)
        mm = vm / math.sqrt(g * r * (t0 - 0.5 * vm * vm / cp))
        if mm > 1:
            p2 = p_sy * (1 + 2 * g / (g + 1) * (mm**2 - 1))
            m2sq = (1 + 0.5 * gm1 * mm**2) / (g * mm**2 - 0.5 * gm1)
        else:
            p2, m2sq = p_sy, mm**2
        return om, p2 * (1 + 0.5 * gm1 * m2sq) ** (g / gm1)

    lo, hi = ps / ((g + 1) / 2) ** (g / gm1), ps * (1 - 1e-9)
    om_c, pc_star = at(lo)
    if pc <= pc_star:
        return om_c
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if at(mid)[1] > pc else (mid, hi)
    return at(0.5 * (lo + hi))[0]


def main() -> None:
    huang = json.loads(FIXTURE.read_text())["cases"]
    curves = chen_curves()
    print("P_c* of the printed chain (and single-change hypotheses) against each author's own numbers")
    print("(a real explanation also collapses Huang's row scatter to table rounding, about 0.5%)")
    print(f"{'hypothesis':38s} {'Huang 31 rows: mean (range)':>28s} {'std':>5s}   {'Chen 8 curves: mean (range)':>28s} in bracket")
    for h in HYPOTHESES:
        eh = []
        for c in huang:
            i = c["inputs"]
            geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
            pc = printed_chain(i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"], h)
            eh.append(pc / c["paper_theory"]["p_c_pa"] - 1)
        ec, inside = [], 0
        for _, pp, geom, lo, hi, _ in curves:
            pc = printed_chain(pp * 1e5, CHEN_T_K, CHEN_PS_PA, CHEN_T_K, geom, AIR_GAMMA, AIR_R, h) / 1e5
            ec.append(pc / (0.5 * (lo + hi)) - 1)
            inside += lo <= pc <= hi
        f = lambda v: f"{100 * sum(v) / len(v):+6.1f}% ({100 * min(v):+5.1f}..{100 * max(v):+5.1f})"  # noqa: E731
        print(f"{h:38s} {f(eh):>28s} {100 * statistics.pstdev(eh):4.1f}%   {f(ec):>28s} {inside}/8")
    print("\nUnit slips: does a constant pressure offset (gauge/absolute, +1 bar, 0.1 MPa) explain it?")
    pairs = {
        "Huang": [
            (
                printed_chain(
                    c["inputs"]["p_g_pa"], c["inputs"]["t_g_k"], c["inputs"]["p_e_pa"], c["inputs"]["t_e_k"],
                    EjectorGeometry(c["inputs"]["area_ratio_nozzle"], c["inputs"]["area_ratio_mix"]),
                    c["inputs"]["gamma"], c["inputs"]["r_gas"],
                ),
                c["paper_theory"]["p_c_pa"],
            )
            for c in huang
        ],
        "Chen": [
            (printed_chain(pp * 1e5, CHEN_T_K, CHEN_PS_PA, CHEN_T_K, geom, AIR_GAMMA, AIR_R), 0.5 * (lo + hi) * 1e5)
            for _, pp, geom, lo, hi, _ in curves
        ],
    }
    for name, pr in pairs.items():
        d = [(y - x) / 1e5 for x, y in pr]
        mean = sum(d) / len(d)
        sd = (sum((v - mean) ** 2 for v in d) / len(d)) ** 0.5
        print(f"  {name}: author - printed = {mean:+.3f} +- {sd:.3f} bar (a slip would be a constant 1.0 bar)")
    i = huang[0]["inputs"]
    geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
    args = (i["p_e_pa"],)
    k = printed_chain(i["p_g_pa"], i["t_g_k"], *args, i["t_e_k"], geom, i["gamma"], i["r_gas"])
    rk = printed_chain(i["p_g_pa"], 1.8 * i["t_g_k"], *args, 1.8 * i["t_e_k"], geom, i["gamma"], i["r_gas"])
    print(f"  temperature in Rankine instead of Kelvin: P_c* {rk / k - 1:+.1e} (only Mach numbers and T_e/T_g enter)")

    print("\nChen's subcritical branch as printed, against his own model droop (omega_ours - omega_Chen)")
    for fig, pp, geom, lo, _, m in curves:
        drop = [r for r in m if r[3] > lo]
        d = [chen_subcritical_omega(pp * 1e5, geom, r[3] * 1e5) - r[4] for r in drop]
        print(f"  Fig. {fig} P_p {pp:g} bar: " + ", ".join(f"{v:+.2f}" for v in d))


if __name__ == "__main__":
    main()
