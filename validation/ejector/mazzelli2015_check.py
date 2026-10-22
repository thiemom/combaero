"""CROSS-SOURCE ACCURACY of the Huang kernel on the Florence R245fa ejector's
third design (Mazzelli & Milazzo 2015, Tables 3-4).

The rig is Eames et al. 2013's (eames2013_check.py) with a third duct,
A_3/A_t = 9.72, on the same 10.2 mm primary throat. Neither ejector belongs to
Huang 1999 or Kracik & Dvorak 2016, so a miss is the model's limitation, not
an implementation bug (docs/VALIDATION_POLICY.md). Nothing is fitted; the
stream states are saturated (data/mazzelli2015_tables.py).

Scored:

* omega on each plateau point, against the measured m_s / m_p with the
  paper's 95 % band -- the first omega on this rig with a measurement band;
* the primary flow, against the measured m_p (0.35 +- 0.03 kg/s): the
  kernel's choked flow (Eq. 1) at the omega-convention gamma, beside the
  real-gas isentropic choke of saturated vapour (a yardstick);
* P_c*, against the bracket the sweep gives (last plateau T_c, first T_c
  below it), for both ``EjectorElement.anchor_closure`` values.

PRE-REGISTERED (before this file was run): interpolating #687's A -> B shift
linearly in A_3/A_t to 9.72 predicts a P_c* bias of about -11.5 % for Huang's
chain and +25.7 % for Kracik & Dvorak. The conditions differ (T_g 89-90 C
here, 100-110 C in Eames 2013), so a miss tests the trend, not the kernel.

The 2014 original-nozzle point (T_g 100 C, T_c* 36.2 C, COP 0.23) is
reported beside it: the kernel cannot tell the two nozzles apart.

Run: ``uv run python validation/ejector/mazzelli2015_check.py``
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from combaero.network._ejector_huang1999 import (
    ETA_P,
    EjectorGeometry,
    choked_mass_flow,
    critical_back_pressure,
    entrainment_ratio,
    huang_critical_back_pressure,
)
from validation.ejector.data import mazzelli2015_tables as mm

#: The kernel's outputs are independent of A_p1/A_t (pinned in
#: test_ejector_eames2013); the third design's real one is (20.2/10.2)^2.
AREA_RATIO_NOZZLE = (mm.GEOMETRY["d_nozzle_exit_mm"] / mm.GEOMETRY["d_throat_mm"]) ** 2

#: Pre-registered P_c* biases at A_3/A_t 9.72 (module docstring).
PREDICTED_PC_BIAS = {"huang": -0.115, "kd": 0.257}


@dataclass(frozen=True)
class PlateauRow:
    point: str
    omega_model: float
    omega_meas: float
    omega_band: float
    m_p_kernel: float  # kg/s, Eq. 1 at the omega-convention gamma
    m_p_real: float  # kg/s, real-gas isentropic choke (yardstick)
    m_p_meas: float
    m_p_band: float

    @property
    def omega_err(self) -> float:
        return self.omega_model / self.omega_meas - 1.0

    @property
    def omega_in_band(self) -> bool:
        return abs(self.omega_model - self.omega_meas) <= self.omega_band


@dataclass(frozen=True)
class CriticalRow:
    sweep: str
    tc_lo: float
    tc_hi: float
    pc_lo: float  # Pa
    pc_hi: float  # Pa
    pc_huang: float  # Pa
    pc_kd: float  # Pa

    @staticmethod
    def _err(pc: float, lo: float, hi: float) -> float:
        """Relative distance outside [lo, hi]; 0 inside."""
        if pc < lo:
            return pc / lo - 1.0
        if pc > hi:
            return pc / hi - 1.0
        return 0.0

    @property
    def pc_mid(self) -> float:
        return 0.5 * (self.pc_lo + self.pc_hi)

    @property
    def pc_huang_err(self) -> float:
        return self._err(self.pc_huang, self.pc_lo, self.pc_hi)

    @property
    def pc_kd_err(self) -> float:
        return self._err(self.pc_kd, self.pc_lo, self.pc_hi)


def _kernel(point_id: str, t_g_c: float, t_e_c: float) -> tuple[float, float, float, float]:
    """(omega, P_c* Huang, P_c* K&D, Eq.-1 primary flow) at a baked state."""
    p_g, p_e, gamma, _flux = mm.BAKED[point_id]
    t_g, t_e = t_g_c + 273.15, t_e_c + 273.15
    geom = EjectorGeometry(area_ratio_nozzle=AREA_RATIO_NOZZLE, area_ratio_mix=mm.area_ratio_mix())
    om = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega
    pc_h = huang_critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, mm.R_SPEC).p_c_pa
    pc_kd = critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, mm.R_SPEC).p_c_pa
    m_p = choked_mass_flow(p_g, t_g, mm.throat_area(), gamma, mm.R_SPEC, ETA_P)
    return om, pc_h, pc_kd, m_p


def score_plateau() -> list[PlateauRow]:
    rows = []
    for sweep in mm.sweeps().values():
        for row in mm.plateau(sweep):
            pid, t_g, t_e, _tc, m_p, u_p, _ms, _us = row
            om, _h, _k, m_kernel = _kernel(pid, t_g, t_e)
            w, u = mm.omega(row)
            rows.append(
                PlateauRow(
                    point=pid,
                    omega_model=om,
                    omega_meas=w,
                    omega_band=u,
                    m_p_kernel=m_kernel,
                    m_p_real=mm.BAKED[pid][3] * mm.throat_area(),
                    m_p_meas=m_p,
                    m_p_band=u_p,
                )
            )
    return rows


def score_critical() -> list[CriticalRow]:
    """Per sweep, both closures at the last plateau point's state against
    the bracket."""
    rows = []
    for name, sweep in mm.sweeps().items():
        lo, hi = mm.critical_bracket(sweep)
        last = next(r for r in sweep if r[3] == lo)
        _om, pc_h, pc_kd, _m = _kernel(last[0], last[1], last[2])
        rows.append(
            CriticalRow(
                sweep=name,
                tc_lo=lo,
                tc_hi=hi,
                pc_lo=mm.P_SAT_COND[lo],
                pc_hi=mm.P_SAT_COND[hi],
                pc_huang=pc_h,
                pc_kd=pc_kd,
            )
        )
    return rows


def score_2014() -> dict[str, float]:
    """The original-nozzle point: kernel vs its printed T_c* and COP."""
    p = mm.POINT_2014
    om, pc_h, pc_kd, _m = _kernel("2014", p["t_g"], p["t_e"])
    pc_ref = mm.P_SAT_COND[p["t_c"]]
    return {
        "omega_model": om,
        "omega_from_cop": mm.omega_from_cop_2014(),
        "pc_huang_err": pc_h / pc_ref - 1.0,
        "pc_kd_err": pc_kd / pc_ref - 1.0,
    }


def report() -> str:
    plat, crit, old = score_plateau(), score_critical(), score_2014()
    lines = [
        "CROSS-SOURCE ACCURACY -- Mazzelli & Milazzo 2015, R245fa third design,",
        f"A_3/A_t {mm.area_ratio_mix():.2f} on the 10.2 mm throat. Measured flows, 95 % bands.",
        f"{'point':9s} {'omega':>6s} {'meas':>6s} {'band':>6s} {'err':>7s}"
        f" {'m_p Eq1':>8s} {'real':>6s} {'meas':>5s}",
    ]
    for r in plat:
        lines.append(
            f"{r.point:9s} {r.omega_model:6.3f} {r.omega_meas:6.3f} {r.omega_band:6.3f}"
            f" {100 * r.omega_err:+6.1f}% {r.m_p_kernel:8.3f} {r.m_p_real:6.3f}"
            f" {r.m_p_meas:5.2f}"
        )
    errs = [r.omega_err for r in plat]
    lines.append(
        f"omega: bias {100 * sum(errs) / len(errs):+.1f}%,"
        f" RMS {100 * math.sqrt(sum(e * e for e in errs) / len(errs)):.1f}%,"
        f" {sum(r.omega_in_band for r in plat)}/{len(plat)} inside the measured band"
    )
    for c in crit:
        lines.append(
            f"P_c* {c.sweep}: T_c* in [{c.tc_lo}, {c.tc_hi}] C = [{c.pc_lo / 1e3:.1f},"
            f" {c.pc_hi / 1e3:.1f}] kPa; Huang {c.pc_huang / 1e3:.1f} kPa"
            f" ({100 * c.pc_huang_err:+.1f}% outside), K&D {c.pc_kd / 1e3:.1f} kPa"
            f" ({100 * c.pc_kd_err:+.1f}% outside)"
        )
    mid_h = [c.pc_huang / c.pc_mid - 1.0 for c in crit]
    mid_k = [c.pc_kd / c.pc_mid - 1.0 for c in crit]
    lines.append(
        f"P_c* vs bracket midpoint: Huang bias {100 * sum(mid_h) / len(mid_h):+.1f}%"
        f" (pre-registered {100 * PREDICTED_PC_BIAS['huang']:+.1f}%),"
        f" K&D {100 * sum(mid_k) / len(mid_k):+.1f}%"
        f" (pre-registered {100 * PREDICTED_PC_BIAS['kd']:+.1f}%)"
    )
    lines.append(
        f"2014 original nozzle (exit {mm.D_NOZZLE_EXIT_2014_MM} mm, T_g 100 C):"
        f" P_c* Huang {100 * old['pc_huang_err']:+.1f}%, K&D {100 * old['pc_kd_err']:+.1f}%;"
        f" omega {old['omega_model']:.3f} vs {old['omega_from_cop']:.3f} from its COP"
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
