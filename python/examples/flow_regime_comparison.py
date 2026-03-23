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

import matplotlib.pyplot as plt
import numpy as np
from plot_utils import show_or_save
from scipy.optimize import brentq as _brentq

import combaero as cb
import combaero.compressible as comp
import combaero.incompressible as incomp
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)

# ------------------------------------------------------------------------------
# Shared conditions
# ------------------------------------------------------------------------------
X = cb.species.dry_air()
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

# ------------------------------------------------------------------------------
# Example 5: Network solver -- validate against direct API
# ------------------------------------------------------------------------------
separator("Example 5: NetworkSolver pipe+orifice -- validate vs direct API")

# Same geometry as Examples 2/3.
A_orif = 1e-4  # 1 cm² throat (same as Example 3)
Cd_orif = 0.65
P_low = 101_325.0  # atmospheric outlet
Y = list(cb.mole_to_mass(X))
area_pipe = 0.25 * math.pi * D**2

velocities = [5.0, 20.0, 50.0, 100.0, 150.0, 200.0]

print(
    f"\n  Network: MassFlowBoundary -> Pipe(L={L}m, D={D * 1e3:.0f}mm, ε={roughness * 1e3:.1f}mm)"
    f" -> Plenum -> Orifice(A={A_orif * 1e4:.0f}cm², Cd={Cd_orif})"
    f" -> PressureBoundary({P_low / 1e5:.1f}bar)"
)
print(f"  Inlet T={T}K,  Outlet P={P_low / 1e5:.2f} bar")
print(f"\n  {'u (m/s)':<10} {'regime':<14} {'P_in direct':>14} {'P_in network':>14} {'err (%)':>9}")
print(f"  {'-' * 10} {'-' * 14} {'-' * 14} {'-' * 14} {'-' * 9}")


def _direct_incompressible(mdot: float) -> float:
    """Self-consistent P_in: iterate on P_junction via brentq."""

    def _residual(P_junc: float) -> float:
        # Orifice: junction → outlet
        rho_j = cb.density(T, P_junc, X)
        mdot_orif = cb.orifice_mdot(P_junc, P_low, A_orif, Cd_orif, rho_j)
        return mdot - mdot_orif

    P_junc = _brentq(_residual, P_low + 1.0, 1e8)

    # Pipe: inlet → junction
    rho_in_est = cb.density(T, P_junc + 1.0, X)
    v_pipe = mdot / (rho_in_est * area_pipe)
    dP_pipe, _, _ = cb.pressure_drop_pipe(
        T,
        P_junc,
        X,
        v_pipe,
        D,
        L,
        roughness,
        "haaland",
    )
    return P_junc + dP_pipe


def _direct_compressible(mdot: float) -> float:
    """Self-consistent P_in using compressible pipe + nozzle."""

    def _residual(P_junc: float) -> float:
        sol_n = comp.nozzle_flow(T, P_junc, X, P_back=P_low, A_eff=A_orif * Cd_orif)
        return mdot - sol_n.mdot

    P_junc = _brentq(_residual, P_low + 1.0, 1e8)

    rho_in_est = cb.density(T, P_junc, X)
    v_pipe = mdot / (rho_in_est * area_pipe)
    sol_p = comp.pipe_flow_rough(T, P_junc, X, u=v_pipe, L=L, D=D, roughness=roughness)
    return P_junc + sol_p.dP


# Collect data for the plot.
plot_u: list[float] = []
plot_Pin_direct_i: list[float] = []
plot_Pin_direct_c: list[float] = []
plot_Pin_net_i: list[float] = []
plot_Pin_net_c: list[float] = []

for u_in in velocities:
    rho_in = cb.density(T, P, X)
    mdot = rho_in * area_pipe * u_in
    plot_u.append(u_in)

    for regime_label, pipe_regime, orif_regime, direct_fn in [
        ("incompressible", "incompressible", "incompressible", _direct_incompressible),
        ("compressible", "compressible", "compressible", _direct_compressible),
    ]:
        # --- Direct API (self-consistent via brentq) ---
        Pin_direct = direct_fn(mdot)

        # --- Network solver ---
        net = FlowNetwork()
        net.add_node(MassFlowBoundary("inlet", m_dot=mdot, T_total=T, Y=Y))
        net.add_node(PlenumNode("junction"))
        net.add_node(PressureBoundary("outlet", P_total=P_low, T_total=T, Y=Y))
        net.add_element(
            PipeElement(
                id="pipe",
                from_node="inlet",
                to_node="junction",
                length=L,
                diameter=D,
                roughness=roughness,
                regime=pipe_regime,
            )
        )
        net.add_element(
            OrificeElement(
                id="orifice",
                from_node="junction",
                to_node="outlet",
                Cd=Cd_orif,
                area=A_orif,
                regime=orif_regime,
            )
        )
        solver = NetworkSolver(net)
        sol_net = solver.solve()
        ok = sol_net["__success__"]
        Pin_net = float(sol_net["inlet.P_total"]) if ok else float("nan")

        err = abs(Pin_net - Pin_direct) / Pin_direct * 100 if ok else float("nan")
        tag = "" if ok else "  [FAILED]"
        print(
            f"  {u_in:<10.0f} {regime_label:<14} "
            f"{Pin_direct:>14.0f} {Pin_net:>14.0f} {err:>8.2f}%{tag}"
        )

        if regime_label == "incompressible":
            plot_Pin_direct_i.append(Pin_direct)
            plot_Pin_net_i.append(Pin_net)
        else:
            plot_Pin_direct_c.append(Pin_direct)
            plot_Pin_net_c.append(Pin_net)

print(
    "\n  Interpretation\n"
    "  --------------\n"
    "  The network solver (pipe + orifice coupled through a plenum) reproduces\n"
    "  the self-consistent brentq reference solution.  Differences arise from\n"
    "  the orifice velocity-of-approach correction and the pipe evaluation\n"
    "  point (inlet vs junction conditions).\n"
)

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

u_arr = np.asarray(plot_u)
ax = axes[0]
ax.plot(u_arr, np.asarray(plot_Pin_direct_i) / 1e5, "s-", label="direct (incomp.)")
ax.plot(u_arr, np.asarray(plot_Pin_net_i) / 1e5, "o--", label="network (incomp.)")
ax.plot(u_arr, np.asarray(plot_Pin_direct_c) / 1e5, "s-", label="direct (comp.)")
ax.plot(u_arr, np.asarray(plot_Pin_net_c) / 1e5, "o--", label="network (comp.)")
ax.set_xlabel("Inlet velocity [m/s]")
ax.set_ylabel("Required inlet pressure [bar]")
ax.set_title("Direct API vs Network Solver")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

ax = axes[1]
err_i = (
    100
    * (np.asarray(plot_Pin_net_i) - np.asarray(plot_Pin_direct_i))
    / np.asarray(plot_Pin_direct_i)
)
err_c = (
    100
    * (np.asarray(plot_Pin_net_c) - np.asarray(plot_Pin_direct_c))
    / np.asarray(plot_Pin_direct_c)
)
ax.plot(u_arr, err_i, "o-", label="incompressible")
ax.plot(u_arr, err_c, "s-", label="compressible")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("Inlet velocity [m/s]")
ax.set_ylabel("Network vs Direct error [%]")
ax.set_title("Network solver validation")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

fig.tight_layout()
show_or_save(fig, "flow_regime_comparison.png")

if __name__ == "__main__":
    pass  # all output produced at module level above
