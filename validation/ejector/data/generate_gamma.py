"""Regenerate the per-suction-condition gamma values in huang1999_tables.py, and
the per-row P_c* validation targets from the digitized T_c* column.

Huang's 1-D model uses ideal-gas constant-gamma isentropic relations. The
entrainment ratio is governed by the entrained flow choking at the hypothetical
throat y-y, so the physically-motivated place to evaluate gamma is that plane,
where the entrained-flow static temperature is

    T_sy = T_e / ((gamma + 1) / 2)      (M_sy = 1, Eq. 10 of the paper).

We evaluate the ideal-gas gamma(T) = cp0(T) / (cp0(T) - R) of R141b from the
CoolProp equation of state and solve the small fixed point (gamma appears on
both sides through T_sy). This yields a clearly-defined, reproducible value per
suction condition rather than a fitted constant.

Tables 3 and 4 also print T_c*, the saturated-vapour temperature corresponding
to each row's critical back pressure P_c* -- the target the reference model's
P_c* chain (Eqs. 9-18) is checked against. P_c* itself is a real-fluid
saturation pressure (not part of the ideal-gas model), so it is looked up from
the same CoolProp EOS rather than approximated.

Run with CoolProp available (not a runtime/test dependency):

    uv run --with CoolProp python validation/ejector/data/generate_gamma.py

Paste the printed values into GAMMA_BY_SUCTION_C and TC_STAR_TO_PC_STAR_PA in
huang1999_tables.py.
"""

from __future__ import annotations

from CoolProp.CoolProp import PropsSI

FLUID = "R141b"
R_UNIV = 8.314462618
R_SPEC = R_UNIV / PropsSI("molar_mass", FLUID)  # J/(kg K)

# T_c* values digitized from Tables 3 & 4 (degrees C), in the row order of
# ALL_ROWS in huang1999_tables.py.
TC_STAR_C: list[float] = [
    # Table 3 (P_e = 0.040 MPa, 8 C)
    31.3,
    33.0,
    33.6,
    34.2,
    36.3,
    37.1,
    38.6,
    41.0,
    42.1,
    33.8,
    37.5,
    38.9,
    28.0,
    30.5,
    33.6,
    35.5,
    26.9,
    29.5,
    32.5,
    # Table 4 (P_e = 0.0473 MPa, 12 C)
    33.1,
    34.2,
    34.5,
    39.3,
    42.5,
    32.0,
    36.0,
    39.5,
    28.9,
    36.0,
    25.7,
    29.2,
]


def ideal_gas_gamma(temperature_k: float) -> float:
    """Ideal-gas gamma = cp0 / (cp0 - R) for R141b at temperature_k."""
    cp0 = PropsSI("Cp0mass", "T", temperature_k, "P", 101325.0, FLUID)
    return cp0 / (cp0 - R_SPEC)


def choke_plane_gamma(t_e_kelvin: float, *, iterations: int = 12) -> float:
    """Ideal-gas gamma at the entrained-flow choking plane y-y."""
    gamma = ideal_gas_gamma(t_e_kelvin)
    for _ in range(iterations):
        t_sy = t_e_kelvin / ((gamma + 1.0) / 2.0)
        gamma = ideal_gas_gamma(t_sy)
    return gamma


def saturation_pressure_pa(t_celsius: float) -> float:
    """R141b saturated-vapour pressure at t_celsius, from the CoolProp EOS."""
    return PropsSI("P", "T", t_celsius + 273.15, "Q", 1.0, FLUID)


if __name__ == "__main__":
    print(f"R141b R_specific = {R_SPEC:.4f} J/(kg K)")
    print("gamma (GAMMA_BY_SUCTION_C):")
    for t_e_c in (8.0, 12.0):
        gamma = choke_plane_gamma(t_e_c + 273.15)
        print(f"    {t_e_c:g}: {gamma:.5f},")
    print("P_c* targets (TC_STAR_TO_PC_STAR_PA, same order as TC_STAR_C):")
    for t_c in TC_STAR_C:
        print(f"    {saturation_pressure_pa(t_c):.1f},  # T_c* = {t_c:g} C")
