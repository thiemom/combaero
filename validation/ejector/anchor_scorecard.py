"""CROSS-SOURCE ACCURACY of the ejector's operating anchors (envelope work after #598).

What the network EjectorElement predicts, on two rigs that are not sources of
its closures:

* the critical-mode entrainment ratio omega (Huang 1999, Eqs. 1-8);
* the critical back pressure P_c* (Kracik & Dvorak 2016's mixing closure);
* the dead-head back pressure P_b0 (the same closure at omega = 0), which is
  the far end of the linear subcritical droop.

The rigs:

* **Henry et al. 2007** (HEFAT2007, Fig. 5): an air ejector, MEASURED.
  A_p1/A_t = 1.86, A_3/A_t = 5.30. The knee and breakdown are read off the
  curve (``FIG5_ANCHORS``). P1 is scored both as absolute and as gauge
  (``P1_READINGS``).
* **Akbarnejad & Ziabasharhagh 2025** (Thermal Eng. 72(1), Tables 1-3): steam,
  CFD (not measured). E1 has A_3/A_t = 5.44; E2 has 33.3, the range where the
  #598 stress draws were refused. Stream temperatures are assumed saturated.

This is ACCURACY, not fidelity: neither rig belongs to Huang or to Kracik &
Dvorak, so a miss is the model's limitation, not an implementation bug
(docs/VALIDATION_POLICY.md). Nothing is fitted. ``recovery_efficiency`` stays
at its unfitted 1.0; matching a rig is the user's job.

Run: ``uv run python validation/ejector/anchor_scorecard.py``
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import combaero as cb
from combaero.network._ejector_huang1999 import (
    EjectorGeometry,
    critical_back_pressure,
    dead_head_back_pressure,
    entrainment_ratio,
)
from combaero.network.ejector_element import choke_plane_gamma
from validation.ejector.data import akbarnejad2025_tables as akb
from validation.ejector.data import henry2007_fig5 as henry

ATM_PA = 101325.0
#: Henry's laboratory: both streams at room temperature (PT100-measured,
#: not printed per point).
HENRY_T_K = 293.15


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


def _model(
    p_g: float, t_g: float, p_e: float, t_e: float, geom: EjectorGeometry, x: list[float]
) -> tuple[float, float, float]:
    gamma, r = _gas(x, t_e)
    omega = entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega
    pc = critical_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_c_pa
    pb0 = dead_head_back_pressure(p_g, t_g, p_e, t_e, geom, gamma, r).p_b0_pa
    return omega, pc, pb0


def henry_rows() -> list[AnchorRow]:
    g = henry.GEOMETRY_MM
    geom = _geometry(g["d_throat"], g["d_exit"], g["D_mix"])
    x = list(cb.species.dry_air())
    rows = []
    for reading in henry.P1_READINGS:
        offset = ATM_PA if reading == "gauge" else 0.0
        for bar in henry.PRIMARY_PRESSURES_BAR:
            a = henry.FIG5_ANCHORS[bar]
            om, pc, pb0 = _model(bar * 1e5 + offset, HENRY_T_K, ATM_PA, HENRY_T_K, geom, x)
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
        om, pc, pb0 = _model(pp * 1e3, tp, ps * 1e3, ts, geom, x)
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
            )
        )
    return rows


def score() -> list[AnchorRow]:
    return henry_rows() + akbarnejad_rows()


def report(rows: list[AnchorRow]) -> str:
    lines = [
        "CROSS-SOURCE ACCURACY -- ejector anchors (neither rig is a source of the closures)",
        f"{'source':34s} {'kind':8s} {'case':16s} {'A3/At':>6s} "
        f"{'omega':>7s} {'P_c*':>7s} {'P_b0':>7s} {'Pb0/Pc* model':>14s} {'ref':>5s}",
    ]
    for r in rows:
        lines.append(
            f"{r.source:34s} {r.kind:8s} {r.case:16s} {r.area_ratio_mix:6.2f} "
            f"{100 * r.omega_err:+6.0f}% {100 * r.pc_err:+6.0f}% {100 * r.pb0_err:+6.0f}% "
            f"{r.window_model:14.2f} {r.window_ref:5.2f}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(report(score()))
