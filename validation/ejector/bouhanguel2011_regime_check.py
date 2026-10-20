"""Bouhanguel, Desevaux & Gavignet 2011 (Int J Refrig 34:1633): which anchor
pair agrees with the flow regime the authors saw. Accuracy, cross-source.

Their air ejector (Table 1): throat d* 8 mm, nozzle exit d 12 mm (A_p1/A_t
2.25, M_p1 2.3), mixing tube D 24 mm (A_3/A_t 9) or 16 mm (A_3/A_t 4),
discharging to atmosphere. Laser-sheet visualisation, so the data are regime
statements, not curves:

- A = 9 at P1 = 4, 5 and 6 bar: no choking of the secondary flow (Fig. 4),
  i.e. subcritical at the atmospheric back pressure (P_c* < Pa);
- A = 4 at P1 = 5.5 bar: the secondary chokes (Fig. 5), i.e. P_c* > Pa;
- measured entrainment U = 0.47 at A = 9, P1 = 5 bar (Fig. 3) and U = 0.50
  at A = 9, P1 = 4.5 bar (Fig. 10).

The secondary pressure P2 is not given: the induced air is drawn from
ambient through a valve, so P2 <= Pa. A measured U bounds it from below. A
point that is not choked sits on or below the critical plateau at its true
P2, so P2 >= P2_min, where the plateau at P2_min equals U. P_c* rises with P2,
so a closure is consistent with "not choked" only if P_c*(P2) < Pa somewhere
in (P2_min, Pa]: that is P_c*(P2_min) < Pa. P1 is read both as absolute and
as gauge, since the paper does not say. Stream temperatures are taken as
ambient (293.15 K, not printed).

Run: uv run python validation/ejector/bouhanguel2011_regime_check.py
"""

from __future__ import annotations

from dataclasses import dataclass

import combaero as cb

PA = 101325.0
T_AMBIENT = 293.15
GAMMA, R_AIR = 1.4, 287.05
AREA_RATIO_NOZZLE = (12.0 / 8.0) ** 2  # 2.25
CLOSURES = {
    "kracik_dvorak": "ejector_critical_back_pressure_and_jacobian",
    "huang1999": "ejector_huang_critical_back_pressure_and_jacobian",
}
# (A_3/A_t, P1 bar, measured U) of the two A = 9 points the paper quantifies
MEASURED_U = [(9.0, 5.0, 0.47), (9.0, 4.5, 0.50)]
# (A_3/A_t, P1 bar) the paper calls choked, with the outlet at Pa
CHOKED = [(4.0, 5.5)]


@dataclass(frozen=True)
class RegimeRow:
    area_ratio_mix: float
    p1_bar: float
    reading: str  # "absolute" or "gauge"
    closure: str
    p2_min_pa: float  # NaN for a choked statement (P2 = Pa there)
    pc_over_pa_lo: float  # P_c*/Pa at P2_min (choked: at Pa)
    pc_over_pa_hi: float  # P_c*/Pa at P2 = Pa
    stated: str  # "subcritical" or "critical"

    @property
    def consistent(self) -> bool:
        if self.stated == "subcritical":
            return self.pc_over_pa_lo < 1.0
        return self.pc_over_pa_hi > 1.0


def _p1(p1_bar: float, reading: str) -> float:
    return p1_bar * 1e5 + (PA if reading == "gauge" else 0.0)


def _plateau_omega(p1: float, p2: float, area_ratio_mix: float) -> float:
    return cb._solver_tools.ejector_entrainment_ratio_and_jacobian(
        p1, T_AMBIENT, p2, T_AMBIENT, AREA_RATIO_NOZZLE, area_ratio_mix, GAMMA
    ).value.omega


def _pc(closure: str, p1: float, p2: float, area_ratio_mix: float) -> float:
    fn = getattr(cb._solver_tools, CLOSURES[closure])
    return fn(
        p1, T_AMBIENT, p2, T_AMBIENT, AREA_RATIO_NOZZLE, area_ratio_mix, GAMMA, R_AIR, 1.0
    ).value.p_c_pa


def _p2_min(p1: float, area_ratio_mix: float, u: float) -> float:
    """The secondary pressure at which the critical plateau equals U."""
    lo, hi = 1e3, PA
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _plateau_omega(p1, mid, area_ratio_mix) < u else (lo, mid)
    return 0.5 * (lo + hi)


def regime_rows() -> list[RegimeRow]:
    rows = []
    for reading in ("absolute", "gauge"):
        for arm, p1_bar, u in MEASURED_U:
            p1 = _p1(p1_bar, reading)
            p2 = _p2_min(p1, arm, u)
            for closure in CLOSURES:
                rows.append(
                    RegimeRow(
                        arm,
                        p1_bar,
                        reading,
                        closure,
                        p2,
                        _pc(closure, p1, p2, arm) / PA,
                        _pc(closure, p1, PA, arm) / PA,
                        "subcritical",
                    )
                )
        for arm, p1_bar in CHOKED:
            p1 = _p1(p1_bar, reading)
            for closure in CLOSURES:
                pc = _pc(closure, p1, PA, arm) / PA
                rows.append(
                    RegimeRow(arm, p1_bar, reading, closure, float("nan"), pc, pc, "critical")
                )
    return rows


def main() -> None:
    print("A_3/A_t  P1 bar  reading   closure        stated       P2_min/Pa  P_c*/Pa      verdict")
    for r in regime_rows():
        p2 = "-" if r.stated == "critical" else f"{r.p2_min_pa / PA:.3f}"
        span = (
            f"{r.pc_over_pa_hi:.2f}"
            if r.stated == "critical"
            else f"{r.pc_over_pa_lo:.2f}..{r.pc_over_pa_hi:.2f}"
        )
        verdict = "consistent" if r.consistent else "CONTRADICTS"
        print(
            f"{r.area_ratio_mix:7.0f}  {r.p1_bar:6.1f}  {r.reading:8s}  {r.closure:13s}  "
            f"{r.stated:11s}  {p2:>9s}  {span:11s}  {verdict}"
        )


if __name__ == "__main__":
    main()
