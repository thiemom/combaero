"""Does real R141b explain Huang 1999's P_c* chain landing ~23% low on his own rows?

Huang's original critical-back-pressure chain (Eqs. 9-18): the y-y state from
Eqs. 1-8, constant-pressure mixing at P_py with phi_m (Eq. 19), energy (Eq. 12),
a normal shock (Eqs. 16-17) and an isentropic diffuser (Eq. 18). It is run three
ways against each row's P_c* (the saturation pressure of the row's measured
T_c*, PC_STAR_PA_BY_TC_STAR_C):

  ideal    : ideal gas, constant gamma (the fixture's choke-plane gamma)
  rg_mix   : y-y pressure and omega from the ideal chain (they reproduce
             Huang's theory omega to 0.8%); stream velocities, mixing, shock
             and diffuser on the CoolProp R141b EOS
  rg_full  : additionally the entrained flow's choke pressure and both mass
             fluxes from the EOS

The vapour is held metastable (gas phase imposed): on CoolProp's R141b EOS the
saturated-vapour line is slightly wet at 8-12 C, so an equilibrium isentrope
reaches about 1% liquid by the y-y plane. Each row's equilibrium quality there
is reported.

Not a test dependency. Run with CoolProp available:

    uv run --with CoolProp python validation/ejector/data/huang_pc_real_gas.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import CoolProp.CoolProp as CP

from combaero.network._ejector_huang1999 import (
    ETA_P,
    ETA_S,
    PHI_P,
    EjectorGeometry,
    entrainment_ratio,
)

FIXTURE = Path(__file__).with_name("huang1999_reference.json")

AS = CP.AbstractState("HEOS", "R141b")  # metastable vapour (gas phase imposed)
AS.specify_phase(CP.iphase_gas)
EQ = CP.AbstractState("HEOS", "R141b")  # phase equilibrium


def phi_m(area_ratio_mix: float) -> float:
    """Huang's Eq. 19 mixing-loss coefficient."""
    if area_ratio_mix > 8.3:
        return 0.80
    return 0.82 if area_ratio_mix > 6.9 else 0.84


def ideal_chain(
    p_g: float, t_g: float, p_e: float, t_e: float, geom: EjectorGeometry, g: float, r: float
) -> float:
    """Eqs. 9-18 on an ideal gas with constant gamma."""
    e = entrainment_ratio(p_g, t_g, p_e, t_e, geom, g)
    om, mpy, ppy = e.omega, e.mach_hypothetical_throat, e.p_mixing_pa
    tpy = t_g / (1.0 + 0.5 * (g - 1.0) * mpy * mpy)
    tsy = t_e * 2.0 / (g + 1.0)
    vpy = mpy * math.sqrt(g * r * tpy)
    vsy = math.sqrt(g * r * tsy)
    vm = phi_m(geom.area_ratio_mix) * (vpy + om * vsy) / (1.0 + om)
    cp = g * r / (g - 1.0)
    h0 = (cp * tpy + 0.5 * vpy**2 + om * (cp * tsy + 0.5 * vsy**2)) / (1.0 + om)
    mm = vm / math.sqrt(g * r * (h0 - 0.5 * vm * vm) / cp)
    if mm > 1.0:
        p3 = ppy * (1.0 + 2.0 * g / (g + 1.0) * (mm * mm - 1.0))
        m3 = math.sqrt((1.0 + 0.5 * (g - 1.0) * mm * mm) / (g * mm * mm - 0.5 * (g - 1.0)))
    else:
        p3, m3 = ppy, mm
    return p3 * (1.0 + 0.5 * (g - 1.0) * m3 * m3) ** (g / (g - 1.0))


def _vapour_at(p: float, key: str, target: float) -> None:
    """Set AS to the metastable vapour at p where AS.<key>() == target. P-T
    updates honour the imposed phase; CoolProp's P-s and P-h flashes do not."""
    EQ.update(CP.PQ_INPUTS, p, 1.0)
    lo, hi = EQ.T() - 30.0, 500.0
    for _ in range(70):
        t = 0.5 * (lo + hi)
        AS.update(CP.PT_INPUTS, p, t)
        if getattr(AS, key)() < target:
            lo = t
        else:
            hi = t
    AS.update(CP.PT_INPUTS, p, 0.5 * (lo + hi))


def _ps(p: float, s: float) -> tuple[float, float, float]:
    _vapour_at(p, "smass", s)
    return AS.hmass(), AS.rhomass(), AS.speed_sound()


def _ph(p: float, h: float) -> tuple[float, float, float]:
    _vapour_at(p, "hmass", h)
    return AS.rhomass(), AS.speed_sound(), AS.smass()


def _sonic(h0: float, s0: float, p0: float) -> tuple[float, float]:
    """Static pressure and mass flux where V = a along the isentrope."""
    lo, hi = 0.2 * p0, 0.999 * p0
    for _ in range(80):
        p = 0.5 * (lo + hi)
        h, _, a = _ps(p, s0)
        if math.sqrt(max(2.0 * (h0 - h), 0.0)) > a:
            lo = p
        else:
            hi = p
    _, rho, a = _ps(p, s0)
    return p, rho * a


def _normal_shock(p1: float, h1: float, rho1: float, v1: float) -> tuple[float, float, float]:
    """Rankine-Hugoniot on the EOS: bisect the downstream density."""
    flux, mom, ht = rho1 * v1, p1 + rho1 * v1 * v1, h1 + 0.5 * v1 * v1

    def state(rho2: float) -> tuple[float, float, float, float]:
        v2 = flux / rho2
        p2, h2 = mom - flux * v2, ht - 0.5 * v2 * v2
        return _ph(p2, h2)[0] - rho2, p2, h2, v2

    lo, hi = 1.0001 * rho1, 30.0 * rho1
    f_lo = state(lo)[0]
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f_mid = state(mid)[0]
        if (f_mid > 0.0) == (f_lo > 0.0):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    _, p2, h2, v2 = state(0.5 * (lo + hi))
    return p2, h2, v2


def real_gas_chain(
    p_g: float, p_e: float, t_e: float, geom: EjectorGeometry, g: float, full: bool
) -> tuple[float, float, float, float]:
    """Eqs. 9-18 on the EOS. Returns (P_c*, omega, primary and entrained
    equilibrium quality at y-y)."""
    EQ.update(CP.PQ_INPUTS, p_g, 1.0)
    h_g, s_g, t_g = EQ.hmass(), EQ.smass(), EQ.T()
    EQ.update(CP.PQ_INPUTS, p_e, 1.0)
    h_e, s_e = EQ.hmass(), EQ.smass()
    if full:
        p_py, flux_sy = _sonic(h_e, s_e, p_e)
        _, flux_t = _sonic(h_g, s_g, p_g)
        h_py, rho_py, _ = _ps(p_py, s_g)
        area_py = PHI_P * flux_t / (rho_py * math.sqrt(2.0 * (h_g - h_py)))  # Eq. 5
        om = (geom.area_ratio_mix - area_py) * math.sqrt(ETA_S) * flux_sy / (math.sqrt(ETA_P) * flux_t)
    else:
        e = entrainment_ratio(p_g, t_g, p_e, t_e, geom, g)
        p_py, om = e.p_mixing_pa, e.omega
    v_sy = math.sqrt(2.0 * (h_e - _ps(p_py, s_e)[0]))
    v_py = math.sqrt(2.0 * (h_g - _ps(p_py, s_g)[0]))
    vm = phi_m(geom.area_ratio_mix) * (v_py + om * v_sy) / (1.0 + om)
    hm = (h_g + om * h_e) / (1.0 + om) - 0.5 * vm * vm
    rho_m, a_m, s_m = _ph(p_py, hm)
    if vm > a_m:
        p3, h3, v3 = _normal_shock(p_py, hm, rho_m, vm)
        s3 = _ph(p3, h3)[2]
    else:
        p3, h3, v3, s3 = p_py, hm, vm, s_m
    h03 = h3 + 0.5 * v3 * v3
    lo, hi = p3, 10.0 * p3
    for _ in range(70):
        p = 0.5 * (lo + hi)
        if _ps(p, s3)[0] < h03:
            lo = p
        else:
            hi = p
    EQ.update(CP.PSmass_INPUTS, p_py, s_g)
    q_p = EQ.Q()
    EQ.update(CP.PSmass_INPUTS, p_py, s_e)
    return 0.5 * (lo + hi), om, q_p, EQ.Q()


def main() -> None:
    cases = json.loads(FIXTURE.read_text())["cases"]
    errs: dict[str, list[float]] = {"ideal": [], "rg_mix": [], "rg_full": []}
    om_errs, q_p, q_s = [], [], []
    print(f"{'row':16s} {'ideal':>7s} {'rg_mix':>7s} {'rg_full':>8s} {'Z_g':>6s}")
    for c in cases:
        i = c["inputs"]
        geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
        target = c["paper_theory"]["p_c_pa"]
        args = (i["p_e_pa"], i["t_e_k"], geom, i["gamma"])
        pc_i = ideal_chain(i["p_g_pa"], i["t_g_k"], *args[:3], i["gamma"], i["r_gas"])
        pc_m, _, qp, qs = real_gas_chain(i["p_g_pa"], *args, full=False)
        pc_f, om_f, _, _ = real_gas_chain(i["p_g_pa"], *args, full=True)
        EQ.update(CP.PQ_INPUTS, i["p_g_pa"], 1.0)
        z_g = EQ.p() / (EQ.rhomass() * i["r_gas"] * EQ.T())
        for key, pc in (("ideal", pc_i), ("rg_mix", pc_m), ("rg_full", pc_f)):
            errs[key].append(pc / target - 1.0)
        om_errs.append(om_f / c["paper_theory"]["omega"] - 1.0)
        q_p.append(qp)
        q_s.append(qs)
        print(
            f"{c['id']:16s} {100 * (pc_i / target - 1):+6.1f}% {100 * (pc_m / target - 1):+6.1f}%"
            f" {100 * (pc_f / target - 1):+7.1f}% {z_g:6.3f}"
        )
    for key, v in errs.items():
        print(f"P_c* {key:8s}: mean {100 * sum(v) / len(v):+5.1f}%  range {100 * min(v):+5.1f}..{100 * max(v):+5.1f}%")
    print(
        f"omega rg_full vs Huang theory: mean {100 * sum(om_errs) / len(om_errs):+5.1f}%"
        f"  range {100 * min(om_errs):+5.1f}..{100 * max(om_errs):+5.1f}%"
    )
    print(f"equilibrium quality at y-y: primary {min(q_p):.3f}-{max(q_p):.3f}, entrained {min(q_s):.3f}-{max(q_s):.3f}")
    property_state_check(cases)


def property_state_check(cases: list[dict]) -> None:
    """Neither Huang nor Chen states at which state C_p and gamma were taken
    ("computed by the NIST database", Chen). Take them as real values at each
    candidate state, C_p then inconsistent with gamma R/(gamma - 1) in Eq. 12."""
    print("\nC_p and gamma at a real state (C_p used in the energy equation, gamma in the mixing)")
    for label, p_key in (("generator", "p_g_pa"), ("evaporator", "p_e_pa"), ("condenser", None)):
        errs = []
        for c in cases:
            i = c["inputs"]
            p = c["paper_theory"]["p_c_pa"] if p_key is None else i[p_key]
            EQ.update(CP.PQ_INPUTS, p, 1.0)
            cp, g = EQ.cpmass(), EQ.cpmass() / EQ.cvmass()
            geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
            e = entrainment_ratio(i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"])
            om, mpy, ppy, r = e.omega, e.mach_hypothetical_throat, e.p_mixing_pa, i["r_gas"]
            tpy = i["t_g_k"] / (1 + 0.5 * (g - 1) * mpy**2)
            tsy = i["t_e_k"] * 2 / (g + 1)
            vpy, vsy = mpy * math.sqrt(g * r * tpy), math.sqrt(g * r * tsy)
            vm = phi_m(geom.area_ratio_mix) * (vpy + om * vsy) / (1 + om)
            tm = ((tpy + om * tsy) * cp + 0.5 * (vpy**2 + om * vsy**2) - 0.5 * (1 + om) * vm**2) / ((1 + om) * cp)
            mm = vm / math.sqrt(g * r * tm)
            p3 = ppy * (1 + 2 * g / (g + 1) * (mm**2 - 1)) if mm > 1 else ppy
            m3sq = (1 + 0.5 * (g - 1) * mm**2) / (g * mm**2 - 0.5 * (g - 1)) if mm > 1 else mm**2
            errs.append(p3 * (1 + 0.5 * (g - 1) * m3sq) ** (g / (g - 1)) / c["paper_theory"]["p_c_pa"] - 1)
        print(f"  {label:10s} sat. vapour: P_c* mean {100 * sum(errs) / len(errs):+5.1f}%"
              f"  range {100 * min(errs):+5.1f}..{100 * max(errs):+5.1f}%")


if __name__ == "__main__":
    main()
