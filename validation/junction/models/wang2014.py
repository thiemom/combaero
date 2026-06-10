"""
Wang 2014 measured K_13 / K_23 for compressible combining flow at 45 deg T-junctions.

Reference: Wenhui Wang, Zhenhua Lu, Kangyao Deng, Shuan Qu,
"An experimental study of compressible combining flow at 45 degree T-junctions",
Proc IMechE Part C: J Mech Eng Sci 229(9):1600-1610, 2015 (received 2014).

No closed-form analytical correlation is provided in the paper; the
authoritative results are Tables 1 and 2 (tabulated K values measured to
within ~1% repeatability). This module exposes:

  - Table 1: K_13 and K_23 at q=0.5, a=1 vs M_3 in {0.1, 0.2, 0.3, 0.4, 0.5, 0.6}
  - Table 2: K_13 and K_23 at M_3=0.5, a=1 vs q in {0, 0.2, 0.5, 0.8, 1}

Notation:
    K_13 = stagnation pressure loss coefficient, branch 1 (lateral) -> 3 (common)
    K_23 = stagnation pressure loss coefficient, branch 2 (straight) -> 3
    a    = area ratio (S_c / S_b in the paper's notation); a=1 = equal areas
    q    = m_dot_1 / m_dot_3 = lateral / common
    M_3  = Mach number in the common branch
"""

from __future__ import annotations

# Wang 2014, Table 1: K_13 and K_23 at q=0.5, a=1.
# Columns: M_3, L_13, K_13, L_23, K_23.
# L = static pressure loss coefficient (eq 4 in the paper), reported for context.
TABLE_1: list[dict[str, float]] = [
    {"M_3": 0.1, "L_13": 0.007, "K_13": 0.133, "L_23": 0.007, "K_23": 0.147},
    {"M_3": 0.2, "L_13": 0.026, "K_13": 0.168, "L_23": 0.026, "K_23": 0.182},
    {"M_3": 0.3, "L_13": 0.061, "K_13": 0.213, "L_23": 0.062, "K_23": 0.226},
    {"M_3": 0.4, "L_13": 0.111, "K_13": 0.215, "L_23": 0.112, "K_23": 0.224},
    {"M_3": 0.5, "L_13": 0.176, "K_13": 0.223, "L_23": 0.178, "K_23": 0.233},
    {"M_3": 0.6, "L_13": 0.278, "K_13": 0.294, "L_23": 0.281, "K_23": 0.302},
]


# Wang 2014, Table 2: K_13, K_23, and K_t (overall) at M_3=0.5, a=1.
# K_t = q*K_13 + (1-q)*K_23  is the overall coefficient.
TABLE_2: list[dict[str, float]] = [
    {"q": 0.0, "K_13": -0.764, "K_23": 0.236, "K_t": 0.236},
    {"q": 0.2, "K_13": -0.229, "K_23": 0.298, "K_t": 0.195},
    {"q": 0.5, "K_13": 0.223, "K_23": 0.233, "K_t": 0.228},
    {"q": 0.8, "K_13": 0.428, "K_23": -0.063, "K_t": 0.331},
    {"q": 1.0, "K_13": 0.438, "K_23": -0.455, "K_t": 0.438},
]


def lookup_table1(M_3: float, K_id: str) -> float:
    """Exact-match lookup in Table 1 (q=0.5, a=1). K_id in {K_13, K_23}.

    Raises ValueError if M_3 is not one of the tabulated values. Use the
    digitized figures (Fig 10 for a=1) for interpolation.
    """
    for row in TABLE_1:
        if abs(row["M_3"] - M_3) < 1e-6:
            return row[K_id]
    raise ValueError(f"M_3={M_3} not in Wang Table 1 (tabulated: {[r['M_3'] for r in TABLE_1]})")


def lookup_table2(q: float, K_id: str) -> float:
    """Exact-match lookup in Table 2 (M_3=0.5, a=1). K_id in {K_13, K_23, K_t}.

    Raises ValueError if q is not one of the tabulated values.
    """
    for row in TABLE_2:
        if abs(row["q"] - q) < 1e-6:
            return row[K_id]
    raise ValueError(f"q={q} not in Wang Table 2 (tabulated: {[r['q'] for r in TABLE_2]})")
