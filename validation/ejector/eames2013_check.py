"""CROSS-SOURCE ACCURACY of the Huang kernel on Eames et al. 2013's R245fa chiller.

Neither ejector belongs to Huang 1999 or to Kracik & Dvorak 2016, so a miss
is the model's limitation, not an implementation bug (docs/VALIDATION_POLICY.md).
Nothing is fitted; the stream states are saturated (data/eames2013_fig8.py).

Scored per Fig. 8 point:

* P_c*: the measured T_c*'s saturation pressure, against both
  ``EjectorElement.anchor_closure`` values: Huang's printed chain
  (``huang_critical_back_pressure``, the default) and Kracik & Dvorak's
  lossless mixing (``critical_back_pressure``).
* omega: the critical-mode entrainment ratio against the omega the plotted
  COP implies. The authors doubt their COP (data module docstring), so this
  column carries their uncertainty, not a measurement band.

The authors' design point is their own 1-D code, so it is reported as a
model-to-model comparison, never scored.

Run: ``uv run python validation/ejector/eames2013_check.py``
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from combaero.network._ejector_huang1999 import (
    EjectorGeometry,
    critical_back_pressure,
    entrainment_ratio,
    huang_critical_back_pressure,
)
from validation.ejector.data import eames2013_fig8 as eames

#: Not in Table 1 and not used by the kernel (its outputs are independent of
#: A_p1/A_t); any supersonic value gives the same numbers.
AREA_RATIO_NOZZLE = 4.0


@dataclass(frozen=True)
class EamesRow:
    point: str
    ejector: str
    area_ratio_mix: float
    omega_model: float
    omega_from_cop: float
    pc_ref: float  # Pa, saturation pressure of the measured T_c*
    pc_huang: float  # Pa, huang_critical_back_pressure (the default anchor)
    pc_kd: float  # Pa, critical_back_pressure (Kracik & Dvorak)

    @property
    def pc_huang_err(self) -> float:
        return self.pc_huang / self.pc_ref - 1.0

    @property
    def pc_kd_err(self) -> float:
        return self.pc_kd / self.pc_ref - 1.0

    @property
    def omega_ratio(self) -> float:
        return self.omega_model / self.omega_from_cop


def _evaluate(
    point_id: str, t_g_c: float, t_e_c: float, area_ratio_mix: float
) -> tuple[float, float, float]:
    p_g, p_e, _pc, _dhg, _dhe, gamma = eames.BAKED[point_id]
    t_g, t_e = t_g_c + 273.15, t_e_c + 273.15
    geom = EjectorGeometry(area_ratio_nozzle=AREA_RATIO_NOZZLE, area_ratio_mix=area_ratio_mix)
    omega = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega
    r = eames.R_SPEC
    pc_h = huang_critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_c_pa
    pc_kd = critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_c_pa
    return omega, pc_h, pc_kd


def score() -> list[EamesRow]:
    rows = []
    for pid, ej, t_g, t_e, _tc, cop in eames.FIG8:
        ar = eames.GEOMETRY[ej]["area_ratio_mix"]
        omega, pc_h, pc_kd = _evaluate(pid, t_g, t_e, ar)
        rows.append(
            EamesRow(
                point=pid,
                ejector=ej,
                area_ratio_mix=ar,
                omega_model=omega,
                omega_from_cop=eames.omega_from_cop(pid, cop),
                pc_ref=eames.BAKED[pid][2],
                pc_huang=pc_h,
                pc_kd=pc_kd,
            )
        )
    return rows


def design_point() -> dict[str, tuple[float, float, float]]:
    """Per ejector, the kernel at the authors' design conditions:
    (omega, P_c* Huang, P_c* Kracik & Dvorak) [-, Pa, Pa]."""
    d = eames.DESIGN
    return {
        ej: _evaluate("design", d["t_g"], d["t_e"], g["area_ratio_mix"])
        for ej, g in eames.GEOMETRY.items()
    }


def _stats(v: list[float]) -> tuple[float, float]:
    return sum(v) / len(v), math.sqrt(sum(x * x for x in v) / len(v))


def report(rows: list[EamesRow]) -> str:
    lines = [
        "CROSS-SOURCE ACCURACY -- Eames et al. 2013, R245fa (measured T_c*; COP doubted by",
        "the authors). P_c* error vs the measured T_c*'s saturation pressure.",
        f"{'point':9s} {'A3/At':>6s} {'Pc Huang':>9s} {'Pc K&D':>7s}"
        f" {'omega':>6s} {'w(COP)':>7s} {'ratio':>6s}",
    ]
    for r in rows:
        lines.append(
            f"{r.point:9s} {r.area_ratio_mix:6.1f} {100 * r.pc_huang_err:+8.1f}%"
            f" {100 * r.pc_kd_err:+6.1f}% {r.omega_model:6.3f} {r.omega_from_cop:7.3f}"
            f" {r.omega_ratio:6.2f}"
        )
    for ej in ("A", "B", None):
        sel = [r for r in rows if ej is None or r.ejector == ej]
        bh, rh = _stats([r.pc_huang_err for r in sel])
        bk, rk = _stats([r.pc_kd_err for r in sel])
        wr = sum(r.omega_ratio for r in sel) / len(sel)
        lines.append(
            f"{'all' if ej is None else 'ejector ' + ej:10s} P_c* Huang bias {100 * bh:+.1f}%"
            f" RMS {100 * rh:.1f}% | K&D bias {100 * bk:+.1f}% RMS {100 * rk:.1f}%"
            f" | omega/omega(COP) {wr:.2f}"
        )
    d = eames.DESIGN
    lines.append(
        f"design point (authors' 1-D code, MODEL vs MODEL): T_g {d['t_g']:g} C, T_e {d['t_e']:g} C,"
        f" T_c {d['t_c']:g} C, their omega {d['omega']}"
    )
    for ej, (om, pc_h, pc_kd) in design_point().items():
        lines.append(
            f"  ejector {ej}: kernel omega {om:.3f}; P_c* Huang {pc_h / 1e3:.0f} kPa,"
            f" K&D {pc_kd / 1e3:.0f} kPa (design condenser {eames.BAKED['design'][2] / 1e3:.0f} kPa)"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(report(score()))
