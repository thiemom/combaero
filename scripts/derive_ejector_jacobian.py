"""
Symbolic derivation of the Huang (1999) / Kracik & Dvorak (2016) ejector
Jacobian: d(omega)/d(p_g, t_g, p_e, t_e) and d(p03)/d(p_g, t_g, p_e, t_e),
the two closed-form chains ported to C++ in include/ejector.h / src/ejector.cpp.

Generates the step-by-step local derivatives used to hand-write the C++
forward-mode accumulation (each C++ intermediate carries its value plus the
4 running partials, mirroring solver_interface.cpp's
orifice_compressible_mdot_and_jacobian style); the C++ code pastes in the
resulting expressions (after manual simplification), same workflow as
scripts/derive_mynard_jacobian.py.

Run:
    uv run python scripts/derive_ejector_jacobian.py

Key simplification (see docs/ejector C++ port plan): `mach_p1` (the primary
nozzle-exit Mach number) is the supersonic root of an isentropic area-Mach
relation that depends ONLY on `area_ratio_nozzle` and `gamma` -- neither is
a Newton-solved unknown for EjectorElement, so mach_p1 has an EXACT zero
derivative w.r.t. p_g, t_g, p_e, t_e. It is therefore treated here as a
precomputed constant symbol (`M_p1`), not differentiated -- no implicit/
root-finding differentiation is needed anywhere in this derivation. Same
treatment for gamma, r_gas, area_ratio_mix, recovery_efficiency, and the
paper-calibrated coefficients (eta_p, eta_s, phi_p): all constants here.

References: see python/combaero/network/_ejector_huang1999.py's docstring
for full paper citations (Huang 1999 Eqs. 1-8; Kracik & Dvorak 2016
Eqs. 7-13).
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

# Differentiation variables (the 4 Newton-solved thermodynamic inputs).
p_g, t_g, p_e, t_e = sp.symbols("p_g t_g p_e t_e", positive=True)

# Constants for this derivation (precomputed once per residual evaluation;
# NOT functions of p_g/t_g/p_e/t_e within this module's scope).
M_p1, gamma, r_gas = sp.symbols("M_p1 gamma r_gas", positive=True)
ar_mix = sp.symbols("ar_mix", positive=True)  # A_3 / A_t (area_ratio_mix)
eta_p, eta_s, phi_p = sp.symbols("eta_p eta_s phi_p", positive=True)
recovery_efficiency = sp.symbols("recovery_efficiency", positive=True)

gm1 = gamma - 1
gp1 = gamma + 1


def area_over_astar(mach: sp.Expr, g: sp.Expr) -> sp.Expr:
    """Isentropic area ratio A/A* (Huang Eq. 2's forward form)."""
    return (1 / mach) * ((2 / (g + 1)) * (1 + sp.Rational(1, 2) * (g - 1) * mach**2)) ** (
        (g + 1) / (2 * (g - 1))
    )


def q_lambda(lam: sp.Expr, g: sp.Expr) -> sp.Expr:
    """Aerodynamic mass-flow function q(lambda) (Kracik & Dvorak Eq. 9)."""
    return (1 - ((g - 1) / (g + 1)) * lam**2) ** (1 / (g - 1)) * ((g + 1) / 2) ** (
        1 / (g - 1)
    ) * lam


# ---------------------------------------------------------------------------
# entrainment_ratio (Huang Eqs. 1-8) -- mirrors _ejector_huang1999.py exactly.
# ---------------------------------------------------------------------------

p_p1 = p_g / (1 + sp.Rational(1, 2) * gm1 * M_p1**2) ** (gamma / gm1)
p_sy = p_e / ((gp1 / 2) ** (gamma / gm1))
p_py = p_sy

core = (1 + sp.Rational(1, 2) * gm1 * M_p1**2) * (p_py / p_p1) ** (-gm1 / gamma)
mach_py_sq = (core - 1) / (sp.Rational(1, 2) * gm1)  # unclamped branch (core > 1 assumed)
mach_py = sp.sqrt(mach_py_sq)

area_py = phi_p * area_over_astar(mach_py, gamma)
area_sy = ar_mix - area_py

omega = (p_e / p_g) * area_sy * sp.sqrt(t_g / t_e) * sp.sqrt(eta_s / eta_p)

# ---------------------------------------------------------------------------
# critical_back_pressure (Kracik & Dvorak Eqs. 7-13) -- continues from the
# entrainment-ratio state above.
# ---------------------------------------------------------------------------

t_py = t_g / (1 + sp.Rational(1, 2) * gm1 * mach_py**2)
t_sy = t_e / (gp1 / 2)
v_py = mach_py * sp.sqrt(gamma * r_gas * t_py)
v_sy = sp.sqrt(gamma * r_gas * t_sy)

c1_star = sp.sqrt(2 * gamma / gp1 * r_gas * t_g)
c2_star = sp.sqrt(2 * gamma / gp1 * r_gas * t_e)
lambda1 = v_py / c1_star
lambda2 = v_sy / c2_star

theta21 = t_e / t_g
gam = omega  # Kracik & Dvorak's Gamma = ms/mp, same definition as omega

z1 = lambda1 + 1 / lambda1
z2 = lambda2 + 1 / lambda2
z3 = (z1 + gam * sp.sqrt(theta21) * z2) / sp.sqrt((1 + gam) * (1 + gam * theta21))
lambda3 = (z3 - sp.sqrt(z3**2 - 4)) / 2  # subsonic root

q1 = q_lambda(lambda1, gamma)
q2 = q_lambda(lambda2, gamma)
q3 = q_lambda(lambda3, gamma)

p03 = (
    p_g
    * sp.sqrt((1 + gam) * (1 + gam * theta21))
    / (1 + (p_g / p_e) * gam * sp.sqrt(theta21) * q1 / q2)
    * q1
    / q3
)
# p_c = recovery_efficiency * p03 -- trivial final scaling, not differentiated here.


def show(name: str, expr: sp.Expr) -> sp.Expr:
    """Print a derivative expression; returns it.

    Deliberately NOT run through sp.simplify(): gamma is a free symbol, so
    every power in this chain has a SYMBOLIC exponent (1/(gamma-1) etc.);
    simplify's pattern matching on symbolic-exponent radicals is
    pathologically slow here (didn't finish in minutes on the full z3/
    lambda3 chain). The unsimplified diff is still exact and instant; the
    numeric cross-check below is the real correctness gate, not the
    printed form.
    """
    print(f"=== {name} ===")
    print(expr)
    print()
    return expr


VARS = [("p_g", p_g), ("t_g", t_g), ("p_e", p_e), ("t_e", t_e)]

print("=" * 70)
print("Ejector Jacobian: d(omega)/d(p_g,t_g,p_e,t_e)")
print("=" * 70)
domega = {}
for name, var in VARS:
    domega[name] = show(f"domega_d{name}", sp.diff(omega, var))

print("=" * 70)
print("Ejector Jacobian: d(p03)/d(p_g,t_g,p_e,t_e)")
print("=" * 70)
dp03 = {}
for name, var in VARS:
    dp03[name] = show(f"dp03_d{name}", sp.diff(p03, var))

# ---------------------------------------------------------------------------
# Numeric cross-check against validation/ejector/data/huang1999_reference.json
# (central-difference Jacobian targets, h = 1e-6*|x|) -- catches derivation
# mistakes here, before any C++ is written.
# ---------------------------------------------------------------------------


def mach_from_area_ratio_supersonic(area_ratio: float, g: float) -> float:
    lo, hi = 1.0 + 1e-12, 50.0

    def a_over_astar(m: float) -> float:
        gm1_ = g - 1.0
        gp1_ = g + 1.0
        return (1.0 / m) * ((2.0 / gp1_) * (1.0 + 0.5 * gm1_ * m * m)) ** (gp1_ / (2.0 * gm1_))

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if a_over_astar(mid) < area_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def check_case(case: dict, omega_fns: dict, p03_fns: dict) -> None:
    inp = case["inputs"]
    jac = case["jacobian"]
    ref = case["reference"]

    subs = {
        p_g: inp["p_g_pa"],
        t_g: inp["t_g_k"],
        p_e: inp["p_e_pa"],
        t_e: inp["t_e_k"],
        gamma: inp["gamma"],
        r_gas: inp["r_gas"],
        ar_mix: inp["area_ratio_mix"],
        eta_p: 0.95,
        eta_s: 0.85,
        phi_p: 0.88,
        M_p1: mach_from_area_ratio_supersonic(inp["area_ratio_nozzle"], inp["gamma"]),
    }

    print(f"--- {case['id']} ---")
    ok = True
    for name, _ in VARS:
        analytic = float(domega[name].subs(subs))
        target = jac[f"domega_d_{name}_pa" if name in ("p_g", "p_e") else f"domega_d_{name}_k"]
        rel = abs(analytic - target) / (abs(target) + 1e-12)
        status = "OK" if rel < 1e-4 else "MISMATCH"
        ok = ok and status == "OK"
        print(f"  domega_d{name}: analytic={analytic:.6e} target={target:.6e} rel={rel:.2e} {status}")

    for name, _ in VARS:
        analytic = float(dp03[name].subs(subs)) * inp["recovery_efficiency"]
        key = f"dp_c_d_{name}_pa" if name in ("p_g", "p_e") else f"dp_c_d_{name}_k"
        target = jac[key]
        rel = abs(analytic - target) / (abs(target) + 1e-12)
        status = "OK" if rel < 1e-4 else "MISMATCH"
        ok = ok and status == "OK"
        print(f"  dp_c_d{name}: analytic={analytic:.6e} target={target:.6e} rel={rel:.2e} {status}")

    omega_val = float(omega.subs(subs))
    p03_val = float(p03.subs(subs))
    rel_omega = abs(omega_val - ref["omega"]) / abs(ref["omega"])
    rel_p03 = abs(p03_val - ref["p_mixed_stagnation_pa"]) / abs(ref["p_mixed_stagnation_pa"])
    print(f"  omega value: {omega_val:.6e} vs {ref['omega']:.6e} rel={rel_omega:.2e}")
    print(f"  p03 value:   {p03_val:.6e} vs {ref['p_mixed_stagnation_pa']:.6e} rel={rel_p03:.2e}")
    ok = ok and rel_omega < 1e-6 and rel_p03 < 1e-6
    print("  ALL OK" if ok else "  *** FAILURES ABOVE ***")
    print()


if __name__ == "__main__":
    ref_path = Path(__file__).resolve().parent.parent / "validation" / "ejector" / "data" / (
        "huang1999_reference.json"
    )
    data = json.loads(ref_path.read_text())
    print("=" * 70)
    print("Numeric cross-check against huang1999_reference.json (CD targets)")
    print("=" * 70)
    for case in data["cases"][:3]:
        check_case(case, domega, dp03)
