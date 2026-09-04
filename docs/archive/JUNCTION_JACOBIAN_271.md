# MPCE Junction Jacobians and Validation: Provenance Record (issue #271)

Why the junction Jacobian work went the way it did. The code owns the formulas;
this owns the reasoning, the measurements that decided things, and the dead
ends. Contract and current state live in [issue #271]; the acceptance tables
are reproduced from `validation.junction.network_runner.run_network`.

## The defect that started it

`MPCEv2Element`'s Jacobian omitted the common-port static-pressure column. The
Mynard loss term `K_i * q_dyn_com` depends on the common port's density, hence
on its static pressure, which is a solver unknown. `ConstantKTeeElement` gained
that term in PR #230 after an FD test caught a `K*q/P` error; the same defect
survived in `MPCEv2Element` because `test_mpce_v2_jacobian.py` FD-checks
`dKQ_dmdot_separating_T` **in isolation** -- the `d/dmdot` block only -- and
nothing pinned the assembled row against the residual it differentiates.

## The derivation, and why the obvious one is wrong

Two paths carry the dependence:

- **(a)** `q_dyn_com = mdot^2/(2*rho*A^2) ~ 1/rho`, so `d(q_dyn)/dP = -q_dyn/P`.
- **(b)** Mynard's `K` is a function of `U = -mdot/(rho*A)`, so it moves with
  `rho` too.

Path (a) alone gives `analytic = -4.015e-04` against `fd = -3.843e-04` -- 4.5%
off, with a plausible-looking sign. `U` carries `P` and `mdot` through the same
`1/rho` factor (`dU/dP = -U/P`, `dU/dmdot = U/mdot`), so
`dK/dP = -(mdot_com/P) * dK/dmdot_com`, which lets the existing `dKQ/dmdot`
column supply (b). With `d(q_dyn)/dmdot_com = 2*q_dyn/mdot_com` the two combine
to

    dR_i/dP_com = sign * [ K_i*q_dyn/P - (mdot_com/P) * dKQ_i/dmdot_com ]

and the leading term ends up **positive**: (b) more than cancels the naive
`-K*q/P`. **This is the trap for the C++ port.** Re-deriving only (a) lands
4.5% off in a way inspection will not catch and FD will. Pinned by
`python/tests/test_mpce_v2_jacobian_rows.py`.

## What the FD-fallback audit found

Instrumenting the branch selection and all three silent-zero paths, over 10,918
residual evaluations of the validation dataset (and 6,135 of the test suite):

| regime | sympy | FD |
|---|---|---|
| separating (K5, K6), all topologies | 100% | 0% |
| joining (K11, K12), all topologies | 0% | 100% |
| overall | 55% | 45% |

**100% of declines are `suppliers != 1`.** The guard's other three conditions
never decide anything on a real topology. What is written as a narrow
special case is in practice a separating/joining switch, and the defect is not
"some cases fall back to FD" but **merge junctions have no analytic Jacobian at
all**.

Of the three silent-zero paths, only `sign_flip` fired in those runs -- 8
times, 0.16% of FD evaluations, all joining.

**The other two are NOT dead code, and reading them as such was a mistake.**
Not firing across a dataset is absence of evidence, not proof of
unreachability; the dataset contains only 3-port junctions with well-behaved
iterates. Checked by construction instead:

| guard | reachable when | verified |
|---|---|---|
| `K_shape` | `_mynard2010` assigns `K` only under `len(U) <= 3`, so **K is None for any junction with more than 3 ports** | 4-port probe returns `K is None` |
| `mynard_raised` | the collector mask is empty (`float(U[Ci][0])` -> `IndexError`) or the supplier mask is (`ValueError`) -- both are plausible transient Newton iterates | all-suppliers raises `IndexError`, all-collectors raises `ValueError` |

The `K_shape` case is the serious one. `MPCEv2Element` does not restrict `N`,
so a 4-port junction is constructible through the public API, and it produces
**finite residuals with an entirely zero `dKQ` Jacobian block** -- every
mass-flow derivative of the loss term absent, silently. That is not a 0.16%
corner; it is a whole class of junction running on a systematically wrong
Jacobian with no error and no warning.

So all three guards are live and must be *handled* in the port, not deleted.
`continue` is the wrong response in every case: it converts an unrepresentable
state into a plausible-looking zero. The right shapes are a one-sided
difference across a sign flip, and a loud failure (or an analytic branch) where
the closure genuinely has no value to give.

## The recurring bug this arc kept surfacing

Bassett indexes each loss coefficient by the mass-flow fraction in **its own
leg** (Table 1; `validation/junction/equivalences.py` states it). K6 and K12
take the lateral fraction; K5/K2 and K11 take the straight one. The `1-q`
transform exists in `equivalences.py` and **nothing on the network path applied
it**, so the same defect appeared three times in three places:

| where | fixed in |
|---|---|
| Tier-1 tests (local copy of the K2 formula) | #276 |
| MPCE-v1 validation adapter | #277 |
| MPCE-v2 joining adapter (K11, and the seeded targets) | #279 |

Each is now pinned by a **symmetry** assertion rather than a stored number:
evaluating the straight coefficient at `q` and the lateral one at `1-q` must
build the same network, so both extracted values must agree exactly. That
formulation cannot rot the way a hard-coded expected value can.

## Measurements that must survive the port

| | separating (K5/K6) | joining (K11/K12) |
|---|---|---|
| MAE | 0.1465 | 0.1270 |
| bias | -0.0138 | -0.0106 |
| converged | 184 / 315 | 368 / 435 |
| median wall time | 5.2 ms | 7.8 ms |

**The FD regime is the better-behaved one.** Joining is more accurate and
converges more often than separating, and pays for it in wall time. The case
for porting is architectural homogeneity and the CLAUDE.md `(f, J)` rule, not a
convergence rescue -- anyone starting the port expecting robustness gains
should recalibrate first.

## Dead end

The Jacobian fix was expected to improve convergence. It does not: roots are
unchanged to three decimals and convergence moves 187 -> 184 of 315, which is
noise on near-threshold cases. The missing entry is small next to the row
scale. Correctness before the port is the whole of the value, and saying so is
better than implying a performance win that the numbers do not support.

## Invariants

- Every analytic partial has an independent FD cross-check:
  `test_mpce_v2_jacobian_rows.py` (both v2 subclasses, 38 cases) and
  `test_multi_port_chamber_jacobian.py` (the v1 C++ base).
- A coefficient is compared only against the leg its source indexes it on.
- A change to a closure is scored against the digitised curves before it lands,
  and reverted with the numbers written down when the data says no.
