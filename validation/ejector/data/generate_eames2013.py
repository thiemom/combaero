"""Regenerate the baked R245fa values in eames2013_fig8.py from CoolProp.

Per point: the saturation pressures of T_g, T_e and the measured T_c*, the
enthalpies that turn the plotted COP into an entrainment ratio, and gamma by the
convention of huang1999_tables.py (ideal-gas cp0 / (cp0 - R) at the
entrained-flow choking plane, T_sy = T_e / ((gamma + 1) / 2), solved as a fixed
point). Rounded to the printed digits.

Run with CoolProp available (not a runtime/test dependency):

    uv run --no-project --with CoolProp python validation/ejector/data/generate_eames2013.py
"""

from __future__ import annotations

from CoolProp.CoolProp import PropsSI

FLUID = "R245fa"
R_SPEC = 8.314462618 / PropsSI("molar_mass", FLUID)

# (id, T_g, T_e, T_c*) [C], as digitised in eames2013_fig8.py, plus the design point.
POINTS = [
    ("A-105-13", 105.0, 13.0, 36.26),
    ("A-105-10", 105.0, 10.0, 35.18),
    ("A-105-7", 105.0, 7.0, 34.48),
    ("A-100-10", 100.0, 10.0, 32.03),
    ("A-100-7", 100.0, 7.0, 31.88),
    ("B-105-13", 105.0, 13.0, 32.62),
    ("B-105-10", 105.0, 10.0, 32.47),
    ("B-105-7", 105.0, 7.0, 31.88),
    ("B-110-13", 110.0, 13.0, 37.31),
    ("B-110-10", 110.0, 10.0, 36.92),
    ("design", 106.0, 4.0, 39.0),
]


def k(c: float) -> float:
    return c + 273.15


def choke_plane_gamma(t_e: float) -> float:
    p_ref = 0.5 * PropsSI("P", "T", t_e, "Q", 1, FLUID)
    gamma = 1.1
    for _ in range(30):
        cp0 = PropsSI("CP0MASS", "T", t_e / ((gamma + 1.0) / 2.0), "P", p_ref, FLUID)
        gamma = cp0 / (cp0 - R_SPEC)
    return gamma


def main() -> None:
    print(f"R_SPEC = {R_SPEC:.4f}")
    for pid, tg, te, tc in POINTS:
        pg = PropsSI("P", "T", k(tg), "Q", 1, FLUID)
        pe = PropsSI("P", "T", k(te), "Q", 1, FLUID)
        pc = PropsSI("P", "T", k(tc), "Q", 1, FLUID)
        hg = PropsSI("H", "T", k(tg), "Q", 1, FLUID)
        he = PropsSI("H", "T", k(te), "Q", 1, FLUID)
        hl = PropsSI("H", "T", k(tc), "Q", 0, FLUID)
        g = choke_plane_gamma(k(te))
        print(
            f'    "{pid}": ({pg:.0f}.0, {pe:.0f}.0, {pc:.0f}.0, '
            f"{hg - hl:.0f}.0, {he - hl:.0f}.0, {g:.5f}),"
        )


if __name__ == "__main__":
    main()
