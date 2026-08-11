# Ejector validation (Huang 1999 + Kracik & Dvorak 2016)

Groundwork and reference model for a supersonic-ejector network element, plus a
language-neutral golden dataset that the eventual C++ implementation (analytic
Jacobian + ctests), the pybind layer, and the GUI all validate against **without
recreating a CoolProp environment**.

## Source

- `../../docs/ejector/Huang_1d_analysis_ejector.pdf` — B.J. Huang, J.M. Chang,
  C.P. Wang, V.A. Petrenko, "A 1-D analysis of ejector performance",
  *Int. J. Refrigeration* 22 (1999) 354-364. Source of `entrainment_ratio`
  (Eqs. 1-8) and of the digitized Tables 3-4 validation data.
- `../../docs/ejector/Development_of_an_Analytical_Method_for.pdf` — J. Kracik,
  V. Dvorak, "Development of an Analytical Method for Predicting Flow in a
  Supersonic Air Ejector", EPJ Web of Conferences 114, 02059 (2016),
  DOI: 10.1051/epjconf/201611402059. Source of the mixing closure
  `critical_back_pressure` implements (their Eqs. 7-13) -- see Accuracy below
  for why this replaced Huang's original P_c* chain.
- `../../docs/ejector/art_10.1007_s40430-025-05536-7.pdf` — S. Akbarnejad,
  M. Ziabasharhagh, "A novel 1D model for the analysis of double-choked
  ejectors validated by CFD simulations", J. Braz. Soc. Mech. Sci. Eng.
  47:253 (2025), DOI: 10.1007/s40430-025-05536-7. CFD confirmation that a
  naive momentum balance across the mixing region is unreliable (evaluated,
  not adopted -- see Accuracy below).
- `../../docs/ejector/Lienhard_One dimensional.pdf` — R.K. McGovern,
  K.V. Bulusu, M.A. Antar, J.H. Lienhard V, "One-dimensional Model of an
  Optimal Ejector and Parametric Study of Ejector Efficiency" (MIT/GWU/KFUPM).
  Independent corroboration of the same momentum-retaining, shock-free mixing
  pattern as Kracik & Dvorak (evaluated, not adopted -- see Accuracy below).
- `../../docs/ejector/Bartosiewicz_ejector_2004.pdf` — independent air-ejector
  CFD + experiment (validation cross-check target; not a model source).

## Model scope

`models/huang1999.py` implements the **critical (double-choked) mode** 1-D
model: primary flow choked at the nozzle throat, entrained flow choked at the
hypothetical throat y-y, constant-pressure mixing (Huang's Eqs. 1-8,
`entrainment_ratio`). In this regime the entrainment ratio omega is fixed and
independent of back pressure -- the flat plateau of the ejector characteristic.

It also implements the critical back pressure P_c* (`critical_back_pressure`),
using Kracik & Dvorak's mixing closure (their Eqs. 7-13) rather than Huang's
original Eqs. 9-18: mass, momentum and energy are solved together directly
from the y-y state to the fully-mixed, subsonic state -- no separate loss
coefficient on momentum, no explicit shock model; the conjugate subsonic root
of the same conservation laws plays the shock's role. A `recovery_efficiency`
parameter (default 1.0, no artificial loss) scales the resulting lossless
mixed stagnation pressure down to P_c*, exposed for a user to tune against
their own data rather than fitted to this project's validation set. Unlike
omega, this needs the specific gas constant R (each stream's own
critical/sonic velocity is referenced to its own stagnation temperature).

`EjectorGeometry` is pure dimensionless area ratios (A_p1/A_t, A_3/A_t), so
Huang's single-centered-nozzle geometry and a peripheral multi-nozzle geometry
(e.g. Kracik & Dvorak 2016's 12-nozzle wind-tunnel ejector) are the same model
-- only the *inputs* differ. `EjectorGeometry.from_count(count, ...)` builds
the aggregate ratios for `count` identical, independently-choking nozzles
sharing one mixing chamber (their own stated assumption): A_p1/A_t comes out
invariant to count, while A_3/A_t shrinks as count grows, since the aggregate
throat area scales with count but the shared mixing chamber does not.
`count=1` reduces exactly to the direct constructor.

Out of scope for now (roadmap): the subcritical / back-flow branches.

## gamma provenance (the only CoolProp step)

Huang's model uses ideal-gas constant-gamma relations. gamma is evaluated per
suction condition at the **entrained-flow choking plane** y-y (the plane that
sets the entrainment), from the CoolProp R141b EOS:

    T_sy = T_e / ((gamma + 1) / 2),   gamma = cp0(T_sy) / (cp0(T_sy) - R)

giving ~1.111 (Table 3: 8 C -> 1.11210; Table 4: 12 C -> 1.11103). Evaluating at
the primary *stagnation* temperature instead fits far worse (mean 2.2% vs 0.8%),
because entrainment is governed by the cold choking plane, not the hot inlet.

These values are **baked into** `data/huang1999_tables.py`; CoolProp is only
needed to regenerate them, which is rare.

P_c*'s validation target uses CoolProp the same way: each row's T_c* (the
saturated-vapour temperature the paper prints alongside its theory A_3/A_t and
omega) is converted to a saturation pressure P_c* via the CoolProp EOS and
baked into `PC_STAR_PA_BY_TC_STAR_C`.

## Files

| File | What |
|------|------|
| `models/huang1999.py` | Reference model: `entrainment_ratio` + `critical_back_pressure`. |
| `data/huang1999_tables.py` | Digitized Tables 3 & 4 (theory column) + baked gamma and P_c* targets. |
| `data/generate_gamma.py` | Regenerate gamma from CoolProp (**needs CoolProp**). |
| `data/generate_reference.py` | Emit the golden JSON + C++ header (no CoolProp). |
| `data/huang1999_reference.json` | Canonical golden fixture (inputs, outputs, Jacobian, theory). |
| `data/huang1999_reference_data.h` | Same data as a `constexpr` array for ctests. |

The Python regression test lives at `python/tests/test_ejector_huang.py` and
consumes the JSON fixture directly.

## Golden fixture contract

Each case in `huang1999_reference.json` (and the `HuangCase` struct in the
header) carries:

- `inputs` — everything the model needs, including the baked `gamma`.
- `reference` — this model's `omega` and intermediates; the C++ **forward solve
  must match these tightly** (same math).
- `jacobian` — central-difference `d(omega)/d(input)`; the target the C++
  **analytic Jacobian must match to finite-difference tolerance**.
- `paper_theory` — Huang's published theory `omega` (physics check, `omega_rel_tol`
  = 3%) and each row's `T_c*`-derived `p_c_pa` target (physics check,
  `pc_rel_tol` = 15%; see Accuracy below).

## Regeneration workflow

```bash
# Rare: only when the fluid/gamma recipe changes (needs CoolProp).
uv run --with CoolProp python validation/ejector/data/generate_gamma.py
#   -> paste values into GAMMA_BY_SUCTION_C in huang1999_tables.py

# Routine: after any change to the model or table (no CoolProp).
uv run python validation/ejector/data/generate_reference.py
#   -> rewrites huang1999_reference.json and huang1999_reference_data.h

uv run pytest python/tests/test_ejector_huang.py
```

## Accuracy

The reference model reproduces Huang's theory column across all 31 rows (both
primary pressures 0.40-0.60 MPa, both suction conditions, both nozzles, omega
0.16-0.86) to a **mean 0.8%, max ~2%**. The residual is the hard floor of a
constant-gamma ideal-gas model (the true ideal-gas gamma drifts 1.091 -> 1.123
across the expansion) plus table round-off -- not an implementation error. That
floor only comes down by dropping the ideal-gas assumption, which is the job of
the production element on combaero's real-fluid thermo.

### P_c*: why the model changed, and what was tried

Huang's original P_c* chain (Eqs. 9-18: a `phi_m`-weighted momentum balance,
then a separate Rankine-Hugoniot shock, then isentropic diffuser recovery)
was systematically **25% low (max 35%)** against each row's `T_c*`-derived
target. Investigation ruled out gamma by direct test (P_c* is monotonic in
gamma over its entire valid range; the gamma -> 1 limit is still ~20-30%
short; re-anchoring gamma's reference temperature moves P_c* under 1%) and
ruled out an equation error (hand-verified against the printed paper; Eqs.
16-17 reproduce classical normal-shock-table values exactly). That left the
momentum treatment itself as the leading suspect: reproducing the target
needed `phi_m` ~ 0.90-0.95, above Huang's own stated calibrated range
(0.80-0.84), even for rows where `phi_m` sat well inside that range.

Four independent models/reformulations were then implemented against the
same y-y starting state and tested against all 31 rows:

| Approach | mean \|err\| | max \|err\| | sign |
|---|---|---|---|
| Huang (original: `phi_m` momentum + shock + isentropic diffuser) | 25% | 35% | low |
| Lienhard/McGovern (zero-loss momentum, static-T energy, no shock) | 9.0% | 16.0% | high |
| **Kracik & Dvorak (zero-loss, stagnation-referenced, no shock) -- adopted** | **6.2%** | **13.2%** | high |
| Akbarnejad (drops momentum, imposes exact-sonic condition) | 32-60% | -- | low |

Lienhard/McGovern and Kracik & Dvorak both retain momentum in full (mass +
momentum + energy solved together) but drop `phi_m` and the separate shock
step -- both beat Huang's chain by 3-4x and, tellingly, **flip sign**
(over-predict where Huang under-predicts): a fully lossless treatment is a
theoretical upper bound, so landing above target while Huang's damped
treatment lands below is exactly the expected bracketing. Kracik & Dvorak's
`lambda`/`q(lambda)`/`z(lambda)` formalism edges out Lienhard/McGovern
because it preserves each stream's stagnation enthalpy along its own
isentropic path before combining, while Lienhard/McGovern's energy equation
explicitly drops the kinetic-energy term (a mass-weighted *static*
temperature average only, confirmed by re-reading their printed equation
directly).

Akbarnejad & Ziabasharhagh's CFD independently measured the same defect from
a different angle: their `Psi_5` metric (the relative error of assuming
momentum is simply conserved across the mixing region) reaches **up to 30%**
at critical back pressure in their own water-vapor ejector -- close to
Huang's ~25-35% gap here, and real, CFD-quantified confirmation that the
momentum treatment specifically is the weak link. Their own proposed fix
(drop momentum entirely; close on an *externally imposed* exact-sonic
condition) was implemented and tested too, but performed worse than Huang's
own chain (32-60% depending on the reading of their ambiguous mixing-pressure
convention) -- their model is explicitly a geometry-free **design** tool (it
never consumes a given area ratio, only inlet conditions), not a
given-geometry forward analysis like this one, and the paper says so
directly. Dropping momentum is not the right direction for this use case.

**Adopted: Kracik & Dvorak's closure**, with `recovery_efficiency` defaulting
to **1.0** (no artificial loss) rather than a value fitted to this dataset --
fitting `phi_m` or `recovery_efficiency` against Huang's own 31 rows would be
exactly the "never fit against measurements" trap this project avoids
elsewhere; a principled default, if one is ever adopted, should come from
further validated literature or real test/CFD data, not a curve fit to this
one dataset. `test_ejector_huang.py` checks P_c* by regression, by trend
(monotonicity in `A_3/A_t`, and `P_e < P_c* < P_g`), and now also by a real
(if loose) physics tolerance against `paper_theory.p_c_pa` (`pc_rel_tol` =
15%, comfortably above the observed 13.2% max).

## Roadmap

1. ~~Add the P_c* / critical-back-pressure chain to the reference model.~~
   Done -- see Accuracy above; uses Kracik & Dvorak's closure, not Huang's
   original Eqs. 9-18.
2. Draft the production 3-port `EjectorElement` (primary-in, secondary-in,
   outlet) in `python/combaero/network/`, real-fluid thermo, analytic `(f, J)`.
3. Port to C++ with analytic Jacobians (derived symbolically with sympy where
   useful) and ctests consuming `huang1999_reference_data.h`'s central-
   difference Jacobian columns to tight tolerance.
4. pybind bindings, in a separate PR.
5. GUI work and validation, last.
