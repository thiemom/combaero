"""Lucas & Koehler 2012 (Int J Refrig 35:1595): is a transcritical CO2 ejector
inside the element's ideal-gas scope? Its driving flow says no.

The rig's motive flow is supercritical CO2 at 73-103 bar and about
500-700 kg/m^3, flashing to two phases in a converging nozzle with a 0.62 mm
throat (Table 1). The element's primary is an ideal-gas choke. Fig. 5's
measured driving mass flow is set against that choke and against a
homogeneous-equilibrium (HEM) choke on the CoolProp CO2 EOS: the maximum of
rho * sqrt(2 (h0 - h)) along the stagnation isentrope.

The measured points were READ BY EYE off Fig. 5, not digitised: about +-0.3 g/s
and +-0.5 bar. That is ample for a factor-of-two scope verdict, and nothing
else is scored from them.

The rig cannot test the anchors either: its outlet pressure is set by the
refrigeration cycle, not swept, so no plateau end or breakdown is measured.

Not a test dependency. Run with CoolProp available:

    uv run --with CoolProp python validation/ejector/data/lucas2012_choke_check.py
"""

from __future__ import annotations

import math

import CoolProp.CoolProp as CP

FLUID = "CO2"
D_THROAT = 0.62e-3  # m, Table 1
A_THROAT = math.pi * D_THROAT * D_THROAT / 4.0

# (gas cooler outlet T_d [C], high-side pressure [bar], driving flow [g/s]),
# Fig. 5, read by eye
FIG5_BY_EYE = [
    (30, 73, 10.0),
    (30, 79, 12.3),
    (35, 81, 10.2),
    (35, 85, 11.5),
    (35, 94, 16.1),
    (40, 95, 13.0),
    (40, 103, 16.8),
]


def ideal_gas_choke(p0: float, t0: float) -> float:
    """Ideal-gas choked mass flow [kg/s] with the zero-pressure cp of CO2 at t0."""
    r = CP.PropsSI("GAS_CONSTANT", FLUID) / CP.PropsSI("M", FLUID)
    cp0 = CP.PropsSI("Cp0mass", "T", t0, "P", 1e5, FLUID)
    g = cp0 / (cp0 - r)
    return (
        A_THROAT
        * p0
        / math.sqrt(r * t0)
        * math.sqrt(g)
        * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    )


def hem_choke(p0: float, t0: float, n: int = 400) -> float:
    """Homogeneous-equilibrium choked mass flow [kg/s] on the CoolProp EOS."""
    h0 = CP.PropsSI("H", "T", t0, "P", p0, FLUID)
    s0 = CP.PropsSI("S", "T", t0, "P", p0, FLUID)
    best = 0.0
    for i in range(1, n):
        p = p0 * (1.0 - 0.85 * i / n)
        h = CP.PropsSI("H", "P", p, "S", s0, FLUID)
        if h0 > h:
            rho = CP.PropsSI("D", "P", p, "S", s0, FLUID)
            best = max(best, rho * math.sqrt(2.0 * (h0 - h)))
    return A_THROAT * best


def main() -> None:
    print("T_d C  p_h bar  rho0 kg/m3  measured g/s  ideal-gas  HEM")
    for td, ph, m in FIG5_BY_EYE:
        p0, t0 = ph * 1e5, td + 273.15
        rho0 = CP.PropsSI("D", "T", t0, "P", p0, FLUID)
        mi, mh = ideal_gas_choke(p0, t0) * 1e3, hem_choke(p0, t0) * 1e3
        print(
            f"{td:5d}  {ph:7d}  {rho0:10.0f}  {m:12.1f}  "
            f"{(mi / m - 1) * 100:+8.0f}%  {(mh / m - 1) * 100:+4.0f}%"
        )


if __name__ == "__main__":
    main()
