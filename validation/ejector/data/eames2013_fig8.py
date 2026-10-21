"""Digitized Fig. 8 and Table 1 of Eames, Milazzo, Paganini & Livi (2013).

I.W. Eames, A. Milazzo, D. Paganini, M. Livi, "The design, manufacture and
testing of a jet-pump chiller for air conditioning and industrial
application", Applied Thermal Engineering 58 (2013) 234-240,
DOI 10.1016/j.applthermaleng.2013.04.028. Not committed (copyrighted).

A 40 kW R245fa jet-pump chiller (Frigel, Florence). Two mixer/diffuser ducts
share ONE movable primary nozzle (Table 1):

* ejector A: area ratio 8.6, diffuser throat 30 mm, carbon fibre;
* ejector B: area ratio 13.2, diffuser throat 37.1 mm, aluminium, throat
  32.5 mm further downstream (302.5 vs 270 mm).

Both imply the same primary throat, 30/sqrt(8.6) = 10.23 mm and
37.1/sqrt(13.2) = 10.21 mm, so Table 1's "area ratio" is A_3/A_t. The nozzle
exit diameter is not given; the Huang kernel does not depend on it (omega and
both P_c* closures are exactly independent of A_p1/A_t). Every point is at the
optimum nozzle position NXP = 0 (Fig. 5).

**Fig. 8** plots the critical condenser temperature T_c* against the
plateau COP for each (T_g, T_e) the rig ran: 5 points per ejector. Generator
and evaporator temperatures are saturation temperatures ("C_sat"). T_g follows
the dotted lines (constant T_g), T_e the dash-dot lines (constant T_e). Panel
a's horizontal dash-dot line at COP 0.135 ends on A-105-13 and marks no other
point; it is not used.

Digitization (2026-10-21): 300 dpi render, each panel calibrated on all of its
gridlines (worst tick residual 0.004 C and 0.0002 COP). Each point is the
centroid of its diamond marker. The user reviewed the points on a scan overlay.

**What is measured, and how well** (the authors' own words, Section 3.1 and
the Conclusions):

* T_c* is the saturation temperature of the condenser pressure at the knee
  of Figs. 6-7. It is a pressure measurement (Keller piezo-resistive
  transmitters), and the cleaner of the two quantities.
* COP is the evaporator duty over the generator's electrical heater power.
  The authors call it uncertain: "because of the rapid variation in the
  chilled water temperature difference across the evaporator, precise
  evaluation of COP is difficult", and the conclusions record "some doubts
  about the accuracy of the experimental data caused by the fluctuations in
  chilled water temperature and relatively small changes in water
  temperature across the evaporator". They also suspect reduced evaporator
  heat exchange from uneven flow distribution. Measured COP is about half
  their 1-D design code's (0.15 average, 0.17 maximum, against 0.31).
* The text quotes those averages for ejector A, but Fig. 8a peaks at 0.135
  and Fig. 8b at 0.25. That inconsistency is the paper's and is left as is.

``omega_from_cop`` converts COP to an entrainment ratio, ignoring the feed
pump (about 2.5 kW against 100+ kW of generator heat):

    omega = COP (h_g - h_l) / (h_e - h_l)

with h_g and h_e saturated vapour at T_g and T_e, and h_l saturated liquid at
T_c*. Heater power that never reaches the refrigerant counts as generator
heat, so this omega is, if anything, low.

**Design point** (Section 2): the authors' 1-D simulation code (Grazzini,
Milazzo & Paganini 2012), generator 106 C, condenser 39 C, evaporator 4 C,
gives omega 0.44 and COP 0.31; ejector A was designed from it. It is a
MODEL, not a measurement: a model-to-model comparison only. Section 3.1
quotes the same COP 0.31 at evaporator 5 C and condenser 40 C.

R245fa values are baked from CoolProp by ``generate_eames2013.py``; gamma
follows huang1999_tables.py's convention (ideal-gas, at the entrained-flow
choking plane).
"""

from __future__ import annotations

#: Table 1: A_3/A_t and the diffuser (mixing) throat diameter [mm].
GEOMETRY: dict[str, dict[str, float]] = {
    "A": {"area_ratio_mix": 8.6, "d_mix_mm": 30.0, "throat_position_mm": 270.0},
    "B": {"area_ratio_mix": 13.2, "d_mix_mm": 37.1, "throat_position_mm": 302.5},
}

#: Fig. 8: (id, ejector, T_g [C], T_e [C], T_c* [C], COP).
FIG8: list[tuple[str, str, float, float, float, float]] = [
    ("A-105-13", "A", 105.0, 13.0, 36.26, 0.1350),
    ("A-105-10", "A", 105.0, 10.0, 35.18, 0.0980),
    ("A-105-7", "A", 105.0, 7.0, 34.48, 0.0850),
    ("A-100-10", "A", 100.0, 10.0, 32.03, 0.1070),
    ("A-100-7", "A", 100.0, 7.0, 31.88, 0.0810),
    ("B-105-13", "B", 105.0, 13.0, 32.62, 0.2480),
    ("B-105-10", "B", 105.0, 10.0, 32.47, 0.1980),
    ("B-105-7", "B", 105.0, 7.0, 31.88, 0.1483),
    ("B-110-13", "B", 110.0, 13.0, 37.31, 0.1701),
    ("B-110-10", "B", 110.0, 10.0, 36.92, 0.0951),
]

#: The authors' design point (their 1-D code, ejector A): T_g, T_c, T_e [C],
#: omega, COP.
DESIGN = {"t_g": 106.0, "t_c": 39.0, "t_e": 4.0, "omega": 0.44, "cop": 0.31}

#: Specific gas constant of R245fa [J/(kg K)].
R_SPEC = 62.0260

#: Per point (and "design"), from CoolProp: P_g, P_e, P_c* [Pa] (saturation at
#: T_g, T_e, T_c*), h_g - h_l, h_e - h_l [J/kg], choke-plane gamma.
BAKED: dict[str, tuple[float, float, float, float, float, float]] = {
    "A-105-13": (1411746.0, 93278.0, 221239.0, 230698.0, 167229.0, 1.08224),
    "A-105-10": (1411746.0, 82417.0, 213267.0, 232151.0, 166430.0, 1.08286),
    "A-105-7": (1411746.0, 72590.0, 208220.0, 233091.0, 165118.0, 1.08350),
    "A-100-10": (1264897.0, 82417.0, 191283.0, 233636.0, 170651.0, 1.08286),
    "A-100-7": (1264897.0, 72590.0, 190282.0, 233836.0, 168599.0, 1.08350),
    "B-105-13": (1411746.0, 93278.0, 195260.0, 235583.0, 172114.0, 1.08224),
    "B-105-10": (1411746.0, 82417.0, 194243.0, 235784.0, 170063.0, 1.08286),
    "B-105-7": (1411746.0, 72590.0, 190282.0, 236572.0, 168599.0, 1.08350),
    "B-110-13": (1571100.0, 93278.0, 229210.0, 231825.0, 165813.0, 1.08224),
    "B-110-10": (1571100.0, 82417.0, 226224.0, 232351.0, 164087.0, 1.08286),
    "design": (1442593.0, 63722.0, 242504.0, 227523.0, 156774.0, 1.08415),
}


def omega_from_cop(point_id: str, cop: float) -> float:
    """Entrainment ratio from a COP, pump work ignored (module docstring)."""
    _pg, _pe, _pc, dh_g, dh_e, _g = BAKED[point_id]
    return cop * dh_g / dh_e
