"""Tables 3-4 of Mazzelli & Milazzo (2015): the Florence R245fa ejector's
third design, measured mass flows with 95 % bands.

F. Mazzelli, A. Milazzo, "Performance analysis of a supersonic ejector cycle
working with R245fa", International Journal of Refrigeration 49 (2015)
79-92, DOI 10.1016/j.ijrefrig.2014.09.020. Read from the accepted manuscript
in the University of Florence repository (FLORE, hdl 2158/924733). Not
committed.

The same 40 kW Frigel chiller as Eames et al. 2013 (eames2013_fig8.py), with
its THIRD ejector. Geometry, Table 1 of the paper (also Table 5.3 of
Grazzini, Milazzo & Mazzelli, *Ejectors for Efficient Refrigeration*,
Springer 2018):

* primary nozzle throat 10.2 mm, exit 20.2 mm, length 66.4 mm;
* diffuser (mixing) throat 31.8 mm, exit 108.3 mm, length 950 mm (86 mm inlet
  region, 272 mm mixing region, 678 mm subsonic diffuser; Fig. 8).

So A_3/A_t = (31.8/10.2)^2 = 9.72, on the SAME 10.2 mm primary throat as
Eames 2013's ducts A (8.65) and B (13.2): a third duct on one nozzle.

**Tables 3 and 4** list, per experimental point, the generator, evaporator
and condenser saturation temperatures and the measured generator
(primary), evaporator (secondary) and condenser mass flows with their 95 %
total uncertainty. Mass flows are measured by equivalence with the external
water circuits' heat fluxes (no refrigerant flow meter), and the uncertainty
analysis (Sec. 2.3: bias from the manufacturers, precision from the standard
deviation of the mean, sequential perturbation through the NIST properties)
gives the bands. The entrainment ratio is omega = m_s / m_p, its band the
two flow bands in quadrature: a MEASUREMENT band, not derived from any model.

The tables are a back-pressure sweep at (nearly) fixed T_g and T_e: the
plateau and the knee. Which points are on the plateau is read from the data
alone (``plateau``): within the measurement band of the sweep's largest
omega.

**The 2014 original-nozzle point** (Milazzo, Rocchetti & Eames, Energy
Procedia 45 (2014) 1245-1254, Sec. 3.4 and Fig. 8): the same diffuser with
the original nozzle (exit 28.8 mm), T_g 100 C, T_e 5 C, critical T_c* about
36.2 C, plateau COP 0.23. Its T_c* is printed in the text; the COP converts
to omega as eames2013_fig8.omega_from_cop does. The authors replaced that
nozzle because it was heavily overexpanded (book Fig. 5.15, CFD).

**Not a reference:** the paper's CFD matches these data only with a fitted
sand-grain wall roughness of 20 um (measured Ka 4-6 um): tuned, so it is not
scored against.

R245fa values are baked from CoolProp by ``generate_mazzelli2015.py``; gamma
follows huang1999_tables.py's convention (ideal-gas, at the entrained-flow
choking plane).
"""

from __future__ import annotations

import math

#: Table 1: diameters [mm].
GEOMETRY = {"d_throat_mm": 10.2, "d_nozzle_exit_mm": 20.2, "d_mix_mm": 31.8}
#: The 2014 nozzle's exit [mm] (book Fig. 5.15).
D_NOZZLE_EXIT_2014_MM = 28.8

#: Tabs. 3-4: (id, T_g, T_e, T_c [C], m_p, u_p, m_s, u_s [kg/s]).
TABLES: list[tuple[str, float, float, float, float, float, float, float]] = [
    ("E5-27.3", 89.0, 5.0, 27.3, 0.35, 0.03, 0.19, 0.01),
    ("E5-27.8", 89.1, 5.1, 27.8, 0.35, 0.03, 0.19, 0.01),
    ("E5-28.3", 89.3, 5.3, 28.3, 0.35, 0.03, 0.18, 0.01),
    ("E5-28.8", 89.6, 5.0, 28.8, 0.35, 0.03, 0.12, 0.01),
    ("E10-28.7", 89.8, 9.9, 28.7, 0.36, 0.03, 0.26, 0.01),
    ("E10-29.4", 89.9, 9.9, 29.4, 0.35, 0.03, 0.24, 0.02),
    ("E10-30.2", 89.9, 10.0, 30.2, 0.35, 0.03, 0.20, 0.02),
    ("E10-30.4", 90.0, 9.8, 30.4, 0.35, 0.03, 0.18, 0.02),
]

#: The 2014 original-nozzle point: T_g, T_e, T_c* [C], plateau COP.
POINT_2014 = {"t_g": 100.0, "t_e": 5.0, "t_c": 36.2, "cop": 0.23}

#: Specific gas constant of R245fa [J/(kg K)].
R_SPEC = 62.0260

#: Per point: P_g, P_e [Pa] (saturation at T_g, T_e), choke-plane gamma, and
#: the real-gas isentropic choked mass flux of saturated vapour at T_g
#: [kg/(m^2 s)] (a physics yardstick, not a shipped model).
BAKED: dict[str, tuple[float, float, float, float]] = {
    "E5-27.3": (982434.0, 66576.0, 1.08393, 4345.7),
    "E5-27.8": (984766.0, 66866.0, 1.08391, 4355.9),
    "E5-28.3": (989441.0, 67451.0, 1.08386, 4376.3),
    "E5-28.8": (996486.0, 66576.0, 1.08393, 4407.0),
    "E10-28.7": (1001203.0, 82074.0, 1.08288, 4427.6),
    "E10-29.4": (1003568.0, 82074.0, 1.08288, 4437.9),
    "E10-30.2": (1003568.0, 82417.0, 1.08286, 4437.9),
    "E10-30.4": (1005937.0, 81731.0, 1.08290, 4448.3),
    "2014": (1264897.0, 66576.0, 1.08393, 5583.3),
}

#: Saturation pressure [Pa] of each condenser temperature [C] used.
P_SAT_COND: dict[float, float] = {
    27.3: 161632.0,
    27.8: 164584.0,
    28.3: 167579.0,
    28.8: 170615.0,
    28.7: 170005.0,
    29.4: 174316.0,
    30.2: 179347.0,
    30.4: 180623.0,
    36.2: 220790.0,
}

#: 2014 point: h_g - h_l, h_e - h_l [J/kg] (saturated vapour at T_g, T_e;
#: saturated liquid at T_c*).
DH_2014 = (228043.0, 161305.0)


def area_ratio_mix() -> float:
    """A_3/A_t of the third design's diffuser on the 10.2 mm throat."""
    return (GEOMETRY["d_mix_mm"] / GEOMETRY["d_throat_mm"]) ** 2


def throat_area() -> float:
    """Primary nozzle throat area [m^2]."""
    return math.pi / 4.0 * (GEOMETRY["d_throat_mm"] * 1e-3) ** 2


def omega(row: tuple) -> tuple[float, float]:
    """Measured omega = m_s / m_p and its 95 % band (the flow bands in
    quadrature)."""
    _pid, _tg, _te, _tc, m_p, u_p, m_s, u_s = row
    w = m_s / m_p
    return w, w * math.hypot(u_p / m_p, u_s / m_s)


def sweeps() -> dict[str, list[tuple]]:
    """The two back-pressure sweeps, in T_c order."""
    out: dict[str, list[tuple]] = {}
    for row in TABLES:
        out.setdefault(row[0].split("-")[0], []).append(row)
    return {k: sorted(v, key=lambda r: r[3]) for k, v in out.items()}


def plateau(sweep: list[tuple]) -> list[tuple]:
    """The sweep's points whose omega is within the measurement band of its
    largest omega: on the critical-mode plateau, by the data alone."""
    w_max, u_max = max((omega(r) for r in sweep), key=lambda t: t[0])
    return [r for r in sweep if w_max - omega(r)[0] <= max(u_max, omega(r)[1])]


def critical_bracket(sweep: list[tuple]) -> tuple[float, float]:
    """(last plateau T_c, first T_c below the plateau) [C]: T_c* lies between."""
    on = plateau(sweep)
    last = max(r[3] for r in on)
    after = [r[3] for r in sweep if r[3] > last]
    return last, min(after)


def omega_from_cop_2014() -> float:
    """The 2014 point's omega from its COP, pump work ignored."""
    dh_g, dh_e = DH_2014
    return POINT_2014["cop"] * dh_g / dh_e
