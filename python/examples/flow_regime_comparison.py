#!/usr/bin/env python3
"""Flow regime comparison: incompressible vs compressible API.

Demonstrates the symmetric submodule API where swapping the import is the
only change needed to switch between flow regimes.  The same loop body runs
identically against both modules and the results are compared side-by-side.

Topics covered:
- Regime-swap pattern (identical call signatures, different physics)
- pipe_flow: Darcy-Weisbach vs Fanno (adiabatic compressible) friction
- pipe_flow_rough: roughness-based friction in both regimes
- orifice_flow (incompressible) vs nozzle_flow (compressible)
- Validity range: when do the two regimes agree / diverge?
"""

from __future__ import annotations

import math

import combaero as ca
import combaero.compressible as comp
import combaero.incompressible as incomp

# ------------------------------------------------------------------------------
# Shared conditions
# ------------------------------------------------------------------------------
X = ca.standard_dry_air_composition()
T = 400.0  # K
P = 200_000.0  # Pa  (2 bar)
D = 0.05  # m   (50 mm pipe)
L = 2.0  # m
f = 0.02  # Darcy friction factor


def separator(title: str) -> None:
    print(f"\n{'=' * 72}")
    print(f"  {title}")
    print("=" * 72)


def row(label: str, incomp_val: float, comp_val: float, unit: str = "") -> None:
    diff_pct = abs(comp_val - incomp_val) / max(abs(incomp_val), 1e-30) * 100
    flag = "  <- diverges" if diff_pct > 5 else ""
    print(
        f"  {label:<22} {incomp_val:>12.4g} {unit:<6}  "
        f"{comp_val:>12.4g} {unit:<6}  d={diff_pct:5.1f}%{flag}"
    )


# ------------------------------------------------------------------------------
# Example 1: pipe_flow -- identical call, two regimes
# ------------------------------------------------------------------------------
separator("Example 1: pipe_flow -- same call, two regimes")

print(f"\n  Conditions: T={T} K, P={P / 1e5:.1f} bar, D={D * 1e3:.0f} mm, L={L} m, f={f}")
print(f"\n  {'Field':<22} {'incompressible':>18}  {'compressible':>18}  {'diff':>8}")
print(f"  {'-' * 22} {'-' * 18}  {'-' * 18}  {'-' * 8}")

MODULES = [incomp, comp]

# The regime-swap loop: identical call, results keyed by sol.regime
results_pipe = {}
for module in MODULES:
    sol = module.pipe_flow(T, P, X, u=10.0, L=L, D=D, f=f)
    results_pipe[sol.regime] = sol

sol_i = results_pipe["incompressible"]
sol_c = results_pipe["compressible"]

row("mdot", sol_i.mdot, sol_c.mdot, "kg/s")
row("v_out", sol_i.v, sol_c.v, "m/s")
row("dP", sol_i.dP, sol_c.dP, "Pa")
row("Re", sol_i.Re, sol_c.Re, "-")
row("rho_in", sol_i.rho, sol_c.rho, "kg/m3")
row("f_avg", sol_i.f, sol_c.f, "-")

print(f"\n  Incompressible: M=n/a  (regime={sol_i.regime!r})")
print(
    f"  Compressible:   M_out={sol_c.M:.5f}  T_out={sol_c.T_out:.2f} K  "
    f"P_out={sol_c.P_out / 1e5:.4f} bar  (regime={sol_c.regime!r})"
)

print(
    "\n  Interpretation\n"
    "  --------------\n"
    "  At low Mach number (M ~ 0.03) the two models agree to < 0.1 %.\n"
    "  The incompressible model is exact in the limit M -> 0.\n"
)

# ------------------------------------------------------------------------------
# Example 2: pipe_flow_rough -- velocity sweep, both regimes
# ------------------------------------------------------------------------------
separator("Example 2: pipe_flow_rough -- velocity sweep (roughness = 0.1 mm)")

roughness = 1e-4  # 0.1 mm absolute roughness

print(f"\n  {'u (m/s)':<10} {'M_out':>8} {'dP_incomp (Pa)':>18} {'dP_comp (Pa)':>16} {'d (%)':>8}")
print(f"  {'-' * 10} {'-' * 8} {'-' * 18} {'-' * 16} {'-' * 8}")

for u in [5.0, 20.0, 50.0, 100.0, 150.0, 200.0]:
    # regime-swap loop
    results = {}
    for module in MODULES:
        sol = module.pipe_flow_rough(T, P, X, u=u, L=L, D=D, roughness=roughness)
        results[sol.regime] = sol

    si = results["incompressible"]
    sc = results["compressible"]
    diff = abs(sc.dP - si.dP) / max(si.dP, 1e-30) * 100
    choked_str = " [CHOKED]" if sc.choked else ""
    print(f"  {u:<10.0f} {sc.M:>8.4f} {si.dP:>18.1f} {sc.dP:>16.1f} {diff:>7.1f}%{choked_str}")

print(
    "\n  Interpretation\n"
    "  --------------\n"
    "  At low velocities (u < 50 m/s, M < 0.15) the two models agree within ~1 %.\n"
    "  Above u = 100 m/s (M ~ 0.3) compressibility causes the models to diverge\n"
    "  noticeably -- the incompressible model underestimates pressure drop.\n"
    "  At very high velocities the Fanno pipe chokes (M -> 1).\n"
)

# ------------------------------------------------------------------------------
# Example 3: orifice_flow vs nozzle_flow -- pressure ratio sweep
# ------------------------------------------------------------------------------
separator("Example 3: orifice_flow (incomp.) vs nozzle_flow (comp.) -- dP sweep")

A = 1e-4  # 1 cm^2 throat area

print(f"\n  {'P_back/P':>10} {'mdot_incomp':>14} {'mdot_comp':>14} {'d (%)':>8} {'choked':>8}")
print(f"  {'-' * 10} {'-' * 14} {'-' * 14} {'-' * 8} {'-' * 8}")

for ratio in [0.99, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50]:
    P_back = P * ratio

    sol_i = incomp.orifice_flow(T, P, X, P_back=P_back, A=A, Cd=0.65)
    sol_c = comp.nozzle_flow(T, P, X, P_back=P_back, A_eff=A * 0.65)

    diff = abs(sol_c.mdot - sol_i.mdot) / max(sol_i.mdot, 1e-30) * 100
    choked_str = "yes" if sol_c.choked else "no"
    print(
        f"  {ratio:>10.2f} {sol_i.mdot:>14.5f} {sol_c.mdot:>14.5f} {diff:>7.1f}%  {choked_str:>8}"
    )

print(
    "\n  Interpretation\n"
    "  --------------\n"
    "  The incompressible orifice equation (Bernoulli + Cd) overestimates mass\n"
    "  flow at large pressure drops because it ignores the density reduction.\n"
    "  At P_back/P = 0.99 (small dP) both models agree closely.\n"
    "  Once the nozzle chokes (P_back/P ~ 0.53 for air at 400 K) the compressible\n"
    "  mass flow saturates while the incompressible model keeps rising -- a clear\n"
    "  sign that the incompressible assumption has broken down.\n"
)

# ------------------------------------------------------------------------------
# Example 4: FlowSolution field contract
# ------------------------------------------------------------------------------
separator("Example 4: FlowSolution field contract -- nan for inapplicable fields")

sol_i = incomp.pipe_flow(T, P, X, u=10.0, L=1.0, D=D, f=f)
sol_c = comp.pipe_flow(T, P, X, u=10.0, L=1.0, D=D, f=f)

print(f"\n  {'Field':<12} {'incompressible':>18}  {'compressible':>18}")
print(f"  {'-' * 12} {'-' * 18}  {'-' * 18}")


def _fmt(v: object) -> str:
    if isinstance(v, float):
        return "nan" if math.isnan(v) else f"{v:.5g}"
    return str(v)


FIELDS = [
    "regime",
    "mdot",
    "v",
    "dP",
    "Re",
    "rho",
    "f",
    "Cd",
    "M",
    "T_out",
    "P_out",
    "h0",
    "choked",
    "L_choke",
]
for field in FIELDS:
    vi = getattr(sol_i, field)
    vc = getattr(sol_c, field)
    print(f"  {field:<12} {_fmt(vi):>18}  {_fmt(vc):>18}")

print(
    "\n  Fields that are not meaningful for a regime carry math.nan (floats) or\n"
    "  False (bools).  Downstream code can always access every attribute without\n"
    "  branching on the regime -- just check math.isnan(sol.M) if needed.\n"
)

if __name__ == "__main__":
    pass  # all output produced at module level above
