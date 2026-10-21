"""Regenerate the baked R245fa values in mazzelli2015_tables.py from CoolProp.

Per point: the saturation pressures of T_g and T_e, gamma by the convention of
huang1999_tables.py (ideal-gas cp0 / (cp0 - R) at the entrained-flow choking
plane, as generate_eames2013.py), and the real-gas isentropic choked mass flux
of saturated vapour at T_g (a physics yardstick for the primary flow, not a
model we ship). Per critical-temperature bracket end: its saturation pressure.
For the 2014 original-nozzle point also the enthalpies that turn its COP into
an entrainment ratio. Rounded to the printed digits.

Run with CoolProp available (not a runtime/test dependency):

    uv run --no-project --with CoolProp python validation/ejector/data/generate_mazzelli2015.py
"""

from __future__ import annotations

import math

from CoolProp.CoolProp import PropsSI

FLUID = "R245fa"
R_SPEC = 8.314462618 / PropsSI("molar_mass", FLUID)

# (id, T_g, T_e) [C] as printed in Tabs. 3-4, plus the 2014 point.
POINTS = [
    ("E5-27.3", 89.0, 5.0),
    ("E5-27.8", 89.1, 5.1),
    ("E5-28.3", 89.3, 5.3),
    ("E5-28.8", 89.6, 5.0),
    ("E10-28.7", 89.8, 9.9),
    ("E10-29.4", 89.9, 9.9),
    ("E10-30.2", 89.9, 10.0),
    ("E10-30.4", 90.0, 9.8),
    ("2014", 100.0, 5.0),
]

#: Condenser saturation temperatures [C] whose pressures are needed: the
#: tabulated ones (bracket ends) and the 2014 T_c*.
T_COND = [27.3, 27.8, 28.3, 28.8, 28.7, 29.4, 30.2, 30.4, 36.2]


def k(c: float) -> float:
    return c + 273.15


def choke_plane_gamma(t_e: float) -> float:
    p_ref = 0.5 * PropsSI("P", "T", t_e, "Q", 1, FLUID)
    gamma = 1.1
    for _ in range(30):
        cp0 = PropsSI("CP0MASS", "T", t_e / ((gamma + 1.0) / 2.0), "P", p_ref, FLUID)
        gamma = cp0 / (cp0 - R_SPEC)
    return gamma


def real_gas_choke_flux(t_g: float) -> float:
    """Max of rho*u along the isentrope from saturated vapour at rest [kg/(m^2 s)]."""
    p0 = PropsSI("P", "T", t_g, "Q", 1, FLUID)
    h0 = PropsSI("H", "T", t_g, "Q", 1, FLUID)
    s0 = PropsSI("S", "T", t_g, "Q", 1, FLUID)
    best = 0.0
    for i in range(1, 4000):
        p = p0 * (1.0 - i / 4000.0 * 0.8)
        h = PropsSI("H", "P", p, "S", s0, FLUID)
        rho = PropsSI("D", "P", p, "S", s0, FLUID)
        best = max(best, rho * math.sqrt(2.0 * (h0 - h)))
    return best


def main() -> None:
    print(f"R_SPEC = {R_SPEC:.4f}")
    for pid, tg, te in POINTS:
        pg = PropsSI("P", "T", k(tg), "Q", 1, FLUID)
        pe = PropsSI("P", "T", k(te), "Q", 1, FLUID)
        g = choke_plane_gamma(k(te))
        flux = real_gas_choke_flux(k(tg))
        print(f'    "{pid}": ({pg:.0f}.0, {pe:.0f}.0, {g:.5f}, {flux:.1f}),')
    for tc in T_COND:
        print(f"    {tc}: {PropsSI('P', 'T', k(tc), 'Q', 0, FLUID):.0f}.0,")
    hg = PropsSI("H", "T", k(100.0), "Q", 1, FLUID)
    he = PropsSI("H", "T", k(5.0), "Q", 1, FLUID)
    hl = PropsSI("H", "T", k(36.2), "Q", 0, FLUID)
    print(f"2014 dh_g {hg - hl:.0f}.0 dh_e {he - hl:.0f}.0")


if __name__ == "__main__":
    main()
