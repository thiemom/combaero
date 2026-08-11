"""Digitized Tables 3 and 4 of Huang et al. (1999).

Comparison of the 1-D model against 11 R141b ejectors. Each row records the
operating condition, the ejector assembly (nozzle letter + constant-area
section letter, e.g. "EH" = nozzle E + section H), the paper's *theory* column
for the required area ratio A_3/A_t and the entrainment ratio omega, and T_c*
(the saturated-vapour temperature of the row's critical back pressure P_c*).
The experiment columns are not needed to reproduce the model and are omitted.

Nozzle area ratios A_p1/A_t are from Table 1:
    nozzle A: d_t = 2.64 mm, d_p1 = 4.50 mm -> A_p1/A_t = 2.905
    nozzle E: d_t = 2.82 mm, d_p1 = 5.10 mm -> A_p1/A_t = 3.271

Table 3 suction condition: P_e = 0.040 MPa (T_e = 8 C).
Table 4 suction condition: P_e = 0.0473 MPa (T_e = 12 C).

Primary is saturated vapour at P_g (no superheat, Section 4), so T_g is the
saturated-vapour temperature T_gs printed in parentheses in the tables.
"""

from __future__ import annotations

from dataclasses import dataclass

# Nozzle exit/throat area ratios (Table 1).
NOZZLE_AREA_RATIO: dict[str, float] = {"A": 2.905, "E": 3.271}

# Ratio of specific heats for R141b, evaluated per suction condition from the
# CoolProp EOS -- a defined, reproducible value rather than a fitted constant.
#
# Huang's model uses ideal-gas constant-gamma isentropic relations, and the
# entrainment ratio is set by the entrained flow choking at the hypothetical
# throat y-y. gamma is therefore evaluated at that plane's static temperature
# T_sy = T_e / ((gamma + 1) / 2), which depends only on the suction temperature
# T_e -- hence one value per table (Table 3: 8 C, Table 4: 12 C). See and
# regenerate with validation/ejector/data/generate_gamma.py.
#
# This choke-plane gamma (~1.111) reproduces the paper's theory column to a mean
# of 0.8% (max 2.0%). Note the ideal-gas gamma drifts 1.091 (hot stagnation) ->
# 1.123 (cold deep expansion); evaluating instead at the primary STAGNATION
# temperature fits far worse (mean 2.2%, max 7%), because ejector entrainment is
# governed by the cold entrained-choking plane, not the hot primary inlet.
GAMMA_BY_SUCTION_C: dict[float, float] = {
    8.0: 1.11210,  # Table 3, T_e = 8 C
    12.0: 1.11103,  # Table 4, T_e = 12 C
}

# R141b saturated-vapour pressure [Pa] at each row's digitized T_c* (deg C),
# from the CoolProp EOS -- the P_c* validation target for the reference
# model's critical-back-pressure chain (Eqs. 9-18). P_c* is a real-fluid
# saturation pressure, not part of the ideal-gas model, so it is looked up
# rather than derived from gamma. Regenerate with
# validation/ejector/data/generate_gamma.py; baked here so downstream code
# (generate_reference.py, tests, future C++/GUI) needs no CoolProp.
PC_STAR_PA_BY_TC_STAR_C: dict[float, float] = {
    31.3: 98698.5,
    33.0: 104768.6,
    33.6: 106979.3,
    34.2: 109226.4,
    36.3: 117382.8,
    37.1: 120611.9,
    38.6: 126852.8,
    41.0: 137357.5,
    42.1: 142392.1,
    33.8: 107724.3,
    37.5: 122252.3,
    38.9: 128130.6,
    28.0: 87704.2,
    30.5: 95939.1,
    35.5: 114221.4,
    26.9: 84261.8,
    29.5: 92575.4,
    32.5: 102953.8,
    34.5: 110363.6,
    39.3: 129849.8,
    42.5: 144257.9,
    32.0: 101163.6,
    33.1: 105134.6,
    36.0: 116189.4,
    39.5: 130716.1,
    28.9: 90602.1,
    25.7: 80628.8,
    29.2: 91584.5,
}


@dataclass(frozen=True)
class HuangRow:
    """One row of Huang Table 3 or 4 (theory column only)."""

    table: str  # "T3" or "T4"
    ejector: str  # assembly, e.g. "EH"; ejector[0] is the nozzle letter
    p_g_mpa: float  # primary stagnation pressure [MPa]
    t_g_c: float  # primary stagnation (saturated-vapour) temperature [C]
    p_e_mpa: float  # suction stagnation pressure [MPa]
    t_e_c: float  # suction stagnation temperature [C]
    a3_at_theory: float  # required area ratio A_3/A_t (theory)
    omega_theory: float  # entrainment ratio omega (theory)
    tc_star_c: float  # saturated-vapour temperature of P_c* (theory), deg C

    @property
    def nozzle(self) -> str:
        return self.ejector[0]

    @property
    def area_ratio_nozzle(self) -> float:
        return NOZZLE_AREA_RATIO[self.nozzle]

    @property
    def gamma(self) -> float:
        """Ideal-gas gamma for this row's suction condition (choke-plane)."""
        return GAMMA_BY_SUCTION_C[self.t_e_c]

    @property
    def pc_star_theory_pa(self) -> float:
        """R141b saturation pressure at T_c* (theory) -- the P_c* target."""
        return PC_STAR_PA_BY_TC_STAR_C[self.tc_star_c]

    @property
    def p_g_pa(self) -> float:
        return self.p_g_mpa * 1.0e6

    @property
    def p_e_pa(self) -> float:
        return self.p_e_mpa * 1.0e6

    @property
    def t_g_k(self) -> float:
        return self.t_g_c + 273.15

    @property
    def t_e_k(self) -> float:
        return self.t_e_c + 273.15


# Full digitization of Table 3 (P_e = 0.040 MPa, 8 C).
TABLE_3: list[HuangRow] = [
    HuangRow("T3", "EH", 0.604, 95, 0.040, 8, 10.87, 0.4627, 31.3),
    HuangRow("T3", "EF", 0.604, 95, 0.040, 8, 9.67, 0.3774, 33.0),
    HuangRow("T3", "AD", 0.604, 95, 0.040, 8, 9.29, 0.3476, 33.6),
    HuangRow("T3", "EE", 0.604, 95, 0.040, 8, 8.89, 0.3253, 34.2),
    HuangRow("T3", "AC", 0.604, 95, 0.040, 8, 8.57, 0.2983, 36.3),
    HuangRow("T3", "ED", 0.604, 95, 0.040, 8, 8.12, 0.2658, 37.1),
    HuangRow("T3", "AG", 0.604, 95, 0.040, 8, 7.38, 0.2144, 38.6),
    HuangRow("T3", "EG", 0.604, 95, 0.040, 8, 7.05, 0.1919, 41.0),
    HuangRow("T3", "AA", 0.604, 95, 0.040, 8, 6.55, 0.1554, 42.1),
    HuangRow("T3", "AC", 0.538, 90, 0.040, 8, 8.53, 0.3552, 33.8),
    HuangRow("T3", "AB", 0.538, 90, 0.040, 8, 6.65, 0.2093, 37.5),
    HuangRow("T3", "AA", 0.538, 90, 0.040, 8, 6.74, 0.2156, 38.9),
    HuangRow("T3", "AD", 0.465, 84, 0.040, 8, 9.34, 0.5215, 28.0),
    HuangRow("T3", "AC", 0.465, 84, 0.040, 8, 8.68, 0.4605, 30.5),
    HuangRow("T3", "AB", 0.465, 84, 0.040, 8, 6.99, 0.3042, 33.6),
    HuangRow("T3", "AA", 0.465, 84, 0.040, 8, 6.79, 0.2880, 35.5),
    HuangRow("T3", "AC", 0.400, 78, 0.040, 8, 8.97, 0.5966, 26.9),
    HuangRow("T3", "AB", 0.400, 78, 0.040, 8, 7.48, 0.4422, 29.5),
    HuangRow("T3", "AA", 0.400, 78, 0.040, 8, 6.62, 0.3525, 32.5),
]

# Full digitization of Table 4 (P_e = 0.0473 MPa, 12 C).
TABLE_4: list[HuangRow] = [
    HuangRow("T4", "EF", 0.604, 95, 0.0473, 12, 10.43, 0.5482, 33.1),
    HuangRow("T4", "EE", 0.604, 95, 0.0473, 12, 9.67, 0.4894, 34.2),
    HuangRow("T4", "AD", 0.604, 95, 0.0473, 12, 9.47, 0.4708, 34.5),
    HuangRow("T4", "EC", 0.604, 95, 0.0473, 12, 7.69, 0.3235, 39.3),
    HuangRow("T4", "AA", 0.604, 95, 0.0473, 12, 6.91, 0.2573, 42.5),
    HuangRow("T4", "AD", 0.538, 90, 0.0473, 12, 9.50, 0.5573, 32.0),
    HuangRow("T4", "AG", 0.538, 90, 0.0473, 12, 8.00, 0.4142, 36.0),
    HuangRow("T4", "AA", 0.538, 90, 0.0473, 12, 7.03, 0.3257, 39.5),
    HuangRow("T4", "AD", 0.465, 84, 0.0473, 12, 9.63, 0.6906, 28.9),
    HuangRow("T4", "AA", 0.465, 84, 0.0473, 12, 7.07, 0.4147, 36.0),
    HuangRow("T4", "AD", 0.400, 78, 0.0473, 12, 9.85, 0.8626, 25.7),
    HuangRow("T4", "AG", 0.400, 78, 0.0473, 12, 8.26, 0.6659, 29.2),
]

ALL_ROWS: list[HuangRow] = TABLE_3 + TABLE_4

# A compact stratified sample spanning both tables, both primary pressures,
# both nozzles, and the full omega range (0.16 - 0.86). Used for the default
# fast test; the full ALL_ROWS set is exercised separately.
STRATIFIED_SAMPLE: list[HuangRow] = [
    TABLE_3[0],  # EH  P_g 0.604, nozzle E, high omega
    TABLE_3[8],  # AA  P_g 0.604, nozzle A, low omega
    TABLE_3[12],  # AD  P_g 0.465, nozzle A
    TABLE_3[16],  # AC  P_g 0.400, nozzle A
    TABLE_4[0],  # EF  P_g 0.604, nozzle E, Table 4 suction
    TABLE_4[5],  # AD  P_g 0.538, nozzle A, Table 4 suction
    TABLE_4[9],  # AD  P_g 0.400, nozzle A, highest omega
]
