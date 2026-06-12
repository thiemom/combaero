"""
Symbolic derivation of the Mynard Unified0D Jacobian for a 3-port T-junction
in separating flow (1 supplier + 2 collectors).

Generates dK_i / d(port_mdot_j) for i, j in {0, 1, 2} and prints them.
Run this offline to inspect the analytical form; the MPCE-v2 element pastes
the resulting expressions back in (after manual simplification).

Run:
    uv run python scripts/derive_mynard_jacobian.py

Assumptions baked into the symbolic derivation (match the canonical T setup
the validation runner uses; the MPCE-v2 element guards each at runtime):
- Port 0 is the supplier (canonical inlet); ports 1, 2 are collectors.
- Angle convention: theta_0 = pi (Mynard inlet), theta_1 = 0 (straight
  collector), theta_2 in (0, pi] (lateral collector).
- pseudo_col_angle = (theta_1 + theta_2) / 2; no quadrant flip needed.
- pseudo_sup_angle = pi - theta_2/2 (computed from the reoriented frame).
- sign(theta_str_reoriented) = -1, sign(theta_bra_reoriented) = +1.
"""

from __future__ import annotations

import sympy as sp

# Port mdots in junction convention (positive = out of junction).
# Inlet: mdot_0 < 0 (flow IN); outlets: mdot_1, mdot_2 > 0 (flow OUT).
m0, m1, m2 = sp.symbols("m0 m1 m2", real=True)

# Densities and areas (constants for the Jacobian; held at current port state).
rho0, rho1, rho2 = sp.symbols("rho0 rho1 rho2", positive=True)
A0, A1, A2 = sp.symbols("A0 A1 A2", positive=True)

# Branch angle (the only variable angle; inlet at pi and straight at 0).
th2 = sp.symbols("th2", positive=True)  # lateral branch angle (radians)

# Velocities (Mynard convention: positive U = supplier into junction).
U0 = -m0 / (rho0 * A0)
U1 = -m1 / (rho1 * A1)
U2 = -m2 / (rho2 * A2)

# Volumetric flows.
Q0 = U0 * A0
Q1 = U1 * A1  # negative (collector)
Q2 = U2 * A2  # negative (collector)
Qtot = Q0  # single supplier

# Flow ratios (positive for collectors).
FR1 = -Q1 / Qtot
FR2 = -Q2 / Qtot

# Reorientation: theta -> theta - (theta_1 + theta_2) / 2.
# theta_1 starts at 0, theta_2 starts at th2.
# Reoriented angles: theta_str = -th2/2, theta_bra = +th2/2.
# Pseudosupplier angle in the reoriented frame: pi - th2/2 (always positive).
psi_sup = sp.pi - th2 / 2

# sign(theta_collector) in reoriented frame (constants for this geometry).
sign1 = -1  # straight collector ends up at negative angle
sign2 = +1  # branch collector ends up at positive angle

# Energy transfer factor per collector.
etf1 = (sp.Rational(8, 10) * (sp.pi - psi_sup) * sign1 - sp.Rational(2, 10)) * (1 - FR1)
etf2 = (sp.Rational(8, 10) * (sp.pi - psi_sup) * sign2 - sp.Rational(2, 10)) * (1 - FR2)

# Pseudosupplier velocity (Mynard line 110: sum(U*Q)/Qtot for single supplier
# this equals U0 itself).
u_pseudo = U0  # single supplier
tpa1 = Qtot / ((1 - etf1) * u_pseudo)
tpa2 = Qtot / ((1 - etf2) * u_pseudo)

# Area ratios.
AR1 = tpa1 / A1
AR2 = tpa2 / A2

# phi = pseudo_sup_angle - theta_collector (in reoriented frame, wrapped to 2*pi).
phi1 = psi_sup - sign1 * th2 / 2  # straight: psi_sup - (-th2/2) = psi_sup + th2/2
phi2 = psi_sup - sign2 * th2 / 2  # branch:   psi_sup - (+th2/2) = psi_sup - th2/2

# Damping factor (avoids singular C at FR -> 0).
damp1 = 1 - sp.exp(-FR1 / sp.Rational(2, 100))
damp2 = 1 - sp.exp(-FR2 / sp.Rational(2, 100))

C1 = damp1 * (1 - (1 / (AR1 * FR1)) * sp.cos(sp.Rational(3, 4) * (sp.pi - phi1)))
C2 = damp2 * (1 - (1 / (AR2 * FR2)) * sp.cos(sp.Rational(3, 4) * (sp.pi - phi2)))

# K (Matlab line 73): K_col = (U_col^2 / Ucom^2) * (2*C_col + Ucom^2/U_col^2 - 1)
Ucom = U0
K1 = (U1**2 / Ucom**2) * (2 * C1 + Ucom**2 / U1**2 - 1)
K2 = (U2**2 / Ucom**2) * (2 * C2 + Ucom**2 / U2**2 - 1)

# Reference dynamic head (common-side).
q_dyn = sp.Rational(1, 2) * rho0 * Ucom**2

# Per-port K * q_dyn (only collectors nonzero).
KQ0 = sp.Integer(0)  # common port has no loss contribution
KQ1 = K1 * q_dyn
KQ2 = K2 * q_dyn


def show(name: str, expr: sp.Expr) -> None:
    """Pretty-print a (simplified) derivative expression."""
    print(f"=== {name} ===")
    print(sp.simplify(expr))
    print()


print("=" * 70)
print("Mynard Unified0D Jacobian for 3-port separating T-junction")
print("Convention: m_i positive = flow OUT of junction (MPCE convention)")
print("=" * 70)
print()

for i, KQ in enumerate([KQ0, KQ1, KQ2]):
    if KQ == 0:
        print(f"dKQ_{i}/dm_*: 0 (supplier port)")
        print()
        continue
    for j, m in enumerate([m0, m1, m2]):
        deriv = sp.diff(KQ, m)
        show(f"dKQ_{i}/dm_{j}", deriv)
