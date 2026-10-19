"""CROSS-SOURCE ACCURACY of the ejector's operating anchors (envelope work after #598).

What the network EjectorElement predicts, on three rigs that are not sources
of its closures:

* the critical-mode entrainment ratio omega (Huang 1999, Eqs. 1-8);
* the critical back pressure P_c* (Kracik & Dvorak 2016's mixing closure);
* the dead-head back pressure P_b0 (the same closure at omega = 0), which is
  the far end of the linear subcritical droop.

The rigs:

* **Henry et al. 2007** (HEFAT2007, Fig. 5): an air ejector, MEASURED.
  A_p1/A_t = 1.86, A_3/A_t = 5.30. The knee and breakdown are read off the
  curve (``FIG5_ANCHORS``). P1 is scored both as absolute and as gauge
  (``P1_READINGS``).
* **Chen et al. 2013** (Int. J. Refrig. 36, Fig. 7a): a large rectangular
  air ejector, MEASURED. A_1/A_t = 1.38, A_2/A_t = 4.46, P_p 3.0-5.0 bar read
  as absolute (the efficiency-free ratio test, ``chen2013_figs``). The knee is
  the midpoint between the last plateau point and the first droop point; the
  curves stop at omega 0.13-0.45, so the breakdown is the last two points
  extrapolated linearly to omega = 0.
* **Akbarnejad & Ziabasharhagh 2025** (Thermal Eng. 72(1), Tables 1-3): steam,
  CFD (not measured). E1 has A_3/A_t = 5.44; E2 has 33.3, the range where the
  #598 stress draws were refused. Stream temperatures are assumed saturated.

A CANDIDATE anchor is scored beside them, and it is not in the model:
``candidate_dead_head``. It treats the dead-head as a sudden expansion: the
primary leaves the nozzle (M_p1, P_p1 on A_p1) into the mixing tube, the
dead-headed suction annulus A_3 - A_p1 sits at the secondary pressure P_e,
and mass, momentum and energy over A_3 give a subsonic exit, recovered
losslessly. This is Keenan & Neumann's zero-secondary-flow constant-area
mixer. It has no fitted constant, and it depends on A_3/A_t, which today's
anchor does not. See the README for why it is not adopted yet.

This is ACCURACY, not fidelity: neither rig belongs to Huang or to Kracik &
Dvorak, so a miss is the model's limitation, not an implementation bug
(docs/VALIDATION_POLICY.md). Nothing is fitted. ``recovery_efficiency`` stays
at its unfitted 1.0; matching a rig is the user's job.

Run: ``uv run python validation/ejector/anchor_scorecard.py``
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb
from combaero.network._ejector_huang1999 import (
    EjectorGeometry,
    _mach_from_area_ratio_supersonic,
    critical_back_pressure,
    dead_head_back_pressure,
    entrainment_ratio,
)
from combaero.network.ejector_element import choke_plane_gamma
from validation.ejector.data import akbarnejad2025_tables as akb
from validation.ejector.data import chen2013_figs as chen
from validation.ejector.data import henry2007_fig5 as henry
from validation.ejector.printed_chain_check import printed_chain

ATM_PA = 101325.0
#: Henry's laboratory: both streams at room temperature (PT100-measured,
#: not printed per point).
HENRY_T_K = 293.15
#: Chen's rectangular rig draws its secondary from ambient; Chen states 20 C.
CHEN_T_K = 293.15


@dataclass(frozen=True)
class AnchorRow:
    source: str  # "Henry 2007 (measured, air, P1 absolute)" etc.
    kind: str  # "measured" or "CFD"
    case: str
    area_ratio_mix: float
    omega_model: float
    omega_ref: float
    pc_model: float  # Pa
    pc_ref: float
    pb0_model: float
    pb0_ref: float
    pb0_candidate: float  # candidate_dead_head, not in the model
    pc_candidate: float  # Huang's printed chain (printed_chain_check), not in the model
    #: Whole measured curve omega(P_b), RMS as a fraction of the plateau, for the
    #: model's linear droop and for the candidate pair (None for CFD rows,
    #: which give only the two anchors).
    curve_rms: float | None = None
    curve_rms_candidate: float | None = None

    @property
    def omega_err(self) -> float:
        return self.omega_model / self.omega_ref - 1.0

    @property
    def pc_err(self) -> float:
        return self.pc_model / self.pc_ref - 1.0

    @property
    def pb0_err(self) -> float:
        return self.pb0_model / self.pb0_ref - 1.0

    @property
    def pb0_candidate_err(self) -> float:
        return self.pb0_candidate / self.pb0_ref - 1.0

    @property
    def pc_candidate_err(self) -> float:
        return self.pc_candidate / self.pc_ref - 1.0

    @property
    def window_model(self) -> float:
        """P_b0 / P_c*: how far the droop runs above the knee."""
        return self.pb0_model / self.pc_model

    @property
    def window_ref(self) -> float:
        return self.pb0_ref / self.pc_ref


def _geometry(d_throat: float, d_exit: float, d_mix: float) -> EjectorGeometry:
    return EjectorGeometry(
        area_ratio_nozzle=(d_exit / d_throat) ** 2, area_ratio_mix=(d_mix / d_throat) ** 2
    )


def _gas(x: list[float], t_e: float) -> tuple[float, float]:
    """gamma and R as the element takes them: at the entrained-flow choke plane."""
    return choke_plane_gamma(t_e, x), float(cb.specific_gas_constant(x))


def candidate_dead_head(
    p_g: float, t_g: float, p_e: float, geom: EjectorGeometry, gamma: float, r_gas: float
) -> float:
    """CANDIDATE P_b0 [Pa], not in the model: the zero-entrainment sudden
    expansion over A_3 (see the module docstring). Areas in units of A_t."""
    g = gamma
    m_p1 = _mach_from_area_ratio_supersonic(geom.area_ratio_nozzle, g)
    p_p1 = p_g / (1.0 + 0.5 * (g - 1.0) * m_p1 * m_p1) ** (g / (g - 1.0))
    a1, a3 = geom.area_ratio_nozzle, geom.area_ratio_mix
    # choked primary through the unit throat
    m = p_g * math.sqrt(g / (r_gas * t_g)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    impulse = p_p1 * a1 * (1.0 + g * m_p1 * m_p1) + p_e * (a3 - a1)

    def per_mass(mach: float) -> float:
        t = t_g / (1.0 + 0.5 * (g - 1.0) * mach * mach)
        return (1.0 + g * mach * mach) / (mach * math.sqrt(g / (r_gas * t)))

    target = impulse / m
    if target <= per_mass(1.0):
        mach = 1.0  # no subsonic exit: the mixed flow chokes
    else:
        lo, hi = 1e-9, 1.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if per_mass(mid) > target else (lo, mid)
        mach = 0.5 * (lo + hi)
    t3 = t_g / (1.0 + 0.5 * (g - 1.0) * mach * mach)
    p3 = m * math.sqrt(r_gas * t3 / g) / (a3 * mach)
    return p3 * (1.0 + 0.5 * (g - 1.0) * mach * mach) ** (g / (g - 1.0))


def _model(
    p_g: float, t_g: float, p_e: float, t_e: float, geom: EjectorGeometry, x: list[float]
) -> tuple[float, float, float, float]:
    gamma, r = _gas(x, t_e)
    omega = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega
    pc = critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_c_pa
    pb0 = dead_head_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_b0_pa
    cand = candidate_dead_head(p_g, t_g, p_e, geom, gamma, r)
    pc_cand = printed_chain(p_g, t_g, p_e, t_e, geom, gamma, r)
    return omega, pc, pb0, cand, pc_cand


def _droop(p_b: float, omega: float, p_c: float, p_b0: float) -> float:
    """The element's characteristic: plateau to P_c*, straight chord to P_b0."""
    if p_b <= p_c:
        return omega
    return max(0.0, omega * (p_b0 - p_b) / (p_b0 - p_c))


def _curve_rms(
    points: list[tuple[float, float]], omega: float, p_c: float, p_b0: float, plateau: float
) -> float:
    """RMS of omega over a measured curve [(P_b in Pa, omega)], over its plateau."""
    e = [_droop(pb, omega, p_c, p_b0) - om for pb, om in points]
    return math.sqrt(sum(v * v for v in e) / len(e)) / plateau


def henry_rows() -> list[AnchorRow]:
    g = henry.GEOMETRY_MM
    geom = _geometry(g["d_throat"], g["d_exit"], g["D_mix"])
    x = list(cb.species.dry_air())
    rows = []
    for reading in henry.P1_READINGS:
        offset = ATM_PA if reading == "gauge" else 0.0
        for bar in henry.PRIMARY_PRESSURES_BAR:
            a = henry.FIG5_ANCHORS[bar]
            om, pc, pb0, cand, pc_cand = _model(
                bar * 1e5 + offset, HENRY_T_K, ATM_PA, HENRY_T_K, geom, x
            )
            pts = [(pb * ATM_PA, w) for pb, w in henry.FIG5_OMEGA_VS_PB[bar]]
            plateau = a["omega_plateau"]
            rows.append(
                AnchorRow(
                    source=f"Henry 2007 (air, P1 {reading})",
                    kind="measured",
                    case=f"P1 = {bar} bar",
                    area_ratio_mix=geom.area_ratio_mix,
                    omega_model=om,
                    omega_ref=a["omega_plateau"],
                    pc_model=pc,
                    pc_ref=a["pb_knee"] * ATM_PA,
                    pb0_model=pb0,
                    pb0_ref=a["pb_breakdown"] * ATM_PA,
                    pb0_candidate=cand,
                    pc_candidate=pc_cand,
                    curve_rms=_curve_rms(pts, om, pc, pb0, plateau),
                    curve_rms_candidate=_curve_rms(pts, om, pc_cand, cand, plateau),
                )
            )
    return rows


def chen_rect_anchors(pp_bar: float) -> tuple[float, float, float]:
    """(omega plateau, knee P_c [bar], breakdown P_b0 [bar]) of one measured
    Fig. 7a curve, by the rules in the module docstring."""
    pts = sorted(
        (r for r in chen.FIG7A if r[1] == "measured" and r[2] == pp_bar), key=lambda r: r[3]
    )
    plateau = pts[0][4]
    n = next(i for i, r in enumerate(pts) if abs(r[4] - plateau) > 0.012)
    knee = 0.5 * (pts[n - 1][3] + pts[n][3])
    (x1, y1), (x2, y2) = (pts[-2][3], pts[-2][4]), (pts[-1][3], pts[-1][4])
    breakdown = x2 + y2 * (x2 - x1) / (y1 - y2)
    return plateau, knee, breakdown


def chen_rect_rows() -> list[AnchorRow]:
    a1, a2 = chen.AREA_RATIOS["7a"]
    geom = EjectorGeometry(area_ratio_nozzle=a1, area_ratio_mix=a2)
    x = list(cb.species.dry_air())
    rows = []
    for pp in sorted({r[2] for r in chen.FIG7A}):
        plateau, knee, breakdown = chen_rect_anchors(pp)
        om, pc, pb0, cand, pc_cand = _model(pp * 1e5, CHEN_T_K, ATM_PA, CHEN_T_K, geom, x)
        pts = [(r[3] * 1e5, r[4]) for r in chen.FIG7A if r[1] == "measured" and r[2] == pp]
        rows.append(
            AnchorRow(
                source="Chen 2013 rect (air, P_p absolute)",
                kind="measured",
                case=f"P_p = {pp:g} bar",
                area_ratio_mix=geom.area_ratio_mix,
                omega_model=om,
                omega_ref=plateau,
                pc_model=pc,
                pc_ref=knee * 1e5,
                pb0_model=pb0,
                pb0_ref=breakdown * 1e5,
                pb0_candidate=cand,
                pc_candidate=pc_cand,
                curve_rms=_curve_rms(pts, om, pc, pb0, plateau),
                curve_rms_candidate=_curve_rms(pts, om, pc_cand, cand, plateau),
            )
        )
    return rows


def akbarnejad_rows() -> list[AnchorRow]:
    x = [0.0] * len(cb.species.dry_air())
    x[cb.species.indices["H2O"]] = 1.0
    rows = []
    for ej, pp, ps, er, pc_ref, pb_ref in akb.CFD_ROWS:
        geom = _geometry(*akb.GEOMETRY_MM[ej])
        tp = akb.T_SAT_C[pp] + 273.15
        ts = akb.T_SAT_C[ps] + 273.15
        om, pc, pb0, cand, pc_cand = _model(pp * 1e3, tp, ps * 1e3, ts, geom, x)
        rows.append(
            AnchorRow(
                source=f"Akbarnejad 2025 {ej} (steam)",
                kind="CFD",
                case=f"pp{pp:g}-ps{ps:g}",
                area_ratio_mix=geom.area_ratio_mix,
                omega_model=om,
                omega_ref=er,
                pc_model=pc,
                pc_ref=pc_ref * 1e3,
                pb0_model=pb0,
                pb0_ref=pb_ref * 1e3,
                pb0_candidate=cand,
                pc_candidate=pc_cand,
            )
        )
    return rows


def score() -> list[AnchorRow]:
    return henry_rows() + chen_rect_rows() + akbarnejad_rows()


def report(rows: list[AnchorRow]) -> str:
    def rms(v: float | None) -> str:
        return "  -  " if v is None else f"{100 * v:4.0f}%"

    lines = [
        "CROSS-SOURCE ACCURACY -- ejector anchors (no rig is a source of the closures)",
        "model = today's closures (Kracik & Dvorak P_c* and P_b0); candidate pair = Huang's",
        "printed P_c* chain + the sudden-expansion P_b0 (neither in the model).",
        f"{'source':34s} {'kind':8s} {'case':16s} {'A3/At':>6s} "
        f"{'omega':>7s} {'P_c*':>7s} {'P_b0':>7s} {'Pb0/Pc*':>8s} {'ref':>5s}"
        f" | {'cand P_c*':>9s} {'cand P_b0':>9s} {'Pb0/Pc*':>8s}"
        f" | {'curve RMS':>9s} {'cand':>5s}",
    ]
    for r in rows:
        lines.append(
            f"{r.source:34s} {r.kind:8s} {r.case:16s} {r.area_ratio_mix:6.2f} "
            f"{100 * r.omega_err:+6.0f}% {100 * r.pc_err:+6.0f}% {100 * r.pb0_err:+6.0f}% "
            f"{r.window_model:8.2f} {r.window_ref:5.2f}"
            f" | {100 * r.pc_candidate_err:+8.0f}% {100 * r.pb0_candidate_err:+8.0f}%"
            f" {r.pb0_candidate / r.pc_candidate:8.2f}"
            f" | {rms(r.curve_rms):>9s} {rms(r.curve_rms_candidate):>5s}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(report(score()))
