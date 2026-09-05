# MPCE Junction: Operating-Point Multiplicity and the Harness Scores (issue #271)

Why the pressure-driven junction solves fail, converge twice, or score
themselves. The code owns the formulas; this owns the reasoning and the
measurements that decided things. Companion to
[JUNCTION_JACOBIAN_271.md](JUNCTION_JACOBIAN_271.md), which owns the Jacobian
and validation-wiring history.

Measured 2026-09-05 on a throwaway branch, scripts under `tmp/` (gitignored),
against `validation.junction.network_runner.run_network` at 2073 records.

## The question that started it

Step 2 left 26 solves the degenerate-iterate fallback could not rescue, and
Bassett's own Fig 7b case was among them: theta=45, psi=3, M ~ 0.03, deep
inside the conditions he measured. A documented case failing inside documented
conditions is the highest-priority kind. The working hypothesis was seeding:
start Newton on the physical branch, and detect and steer away from
non-physical ones.

**The hypothesis was wrong, and two claims already in this repository were
wrong with it.** The failures are not basin failures. Most are boundary-value
problems with no solution, and the rest have two solutions, both physical.

## What the boundary conditions actually impose

The three validation topologies constrain the junction very differently, and
the difference decides everything below.

| topology | imposed | free | what a root must satisfy |
|---|---|---|---|
| `imposed_q` | both mass flows | pressures | nothing: q is a boundary condition, the answer is single-valued |
| `mfb_two_pb` | `m_com`, both outlet `Pt` | the split, `Pt_com` | `K_lat(q) - K_str(q) = D_target` |
| `three_pb` | every `Pt` | the split and the flow level | `K_lat(q) / K_str(q) = R_target`, then the level follows |

`imposed_q` is the only topology that pins the operating point. The other two
impose a *relation between the coefficients* and let the junction find its own
q. That is production-realistic, and it is exactly why they are in the harness.
It also means the converged q is an output, not the target.

## Finding 1: the model's coefficient difference is not monotonic

`D(q) = K_lat_model(q) - K_str_model(q)` is U-shaped. Mapped through the
`imposed_q` topology (where q is a boundary condition, so `D` is well defined):

| geometry | min D | at q |
|---|---|---|
| psi=3, theta=45 | 0.551 | 0.27 |
| psi=1, theta=90 | 0.292 | 0.97 |
| psi=1, theta=120 | 1.036 | 0.97 |

A U-shaped `D` splits the `mfb_two_pb` boundary-value problem into three
regimes, and the curve predicts every observed outcome:

| case | D_target | prediction | observed |
|---|---|---|---|
| psi=3, th=45, q=0.2 | 0.122 | below min: **no root** | "not making good progress" |
| psi=3, th=45, q=0.4 | 0.385 | below min: **no root** | "not making good progress" |
| psi=3, th=45, q=0.6 | 1.287 | one root, q=0.60 | converged at q=0.60 |
| psi=3, th=45, q=0.8 | 2.829 | one root, q=0.86 | converged at q=0.857 |
| psi=1, th=120, q=0.406 | 1.109 | **two roots**, 0.060 / 0.940 | 0.055 (default seed) / 0.939 (better seed) |
| psi=1, th=120, q=0.490 | 1.235 | **two roots**, 0.120 / 0.880 | 0.129 / 0.870 |
| psi=1, th=120, q=0.637 | 1.456 | **two roots**, 0.320 / 0.680 | 0.324 / 0.675 |

So the failing Fig 7b points are **infeasible, not mis-seeded**. No initial
guess reaches a root that does not exist, and "the iteration is not making good
progress" is the solver reporting that honestly. Classified across the whole
`mfb_two_pb` population: of 62 non-converged records, 31 are infeasible by
construction and 22 are feasible failures. (The same classification was run for
`three_pb` with the *difference* criterion, which is the wrong one for that
topology -- its condition is the ratio -- so those numbers are not reported.)

**Why the minimum sits so high is the K_straight gap (#272).** For psi=3 the
model's `K_straight` spans -0.13 to +0.93 where Bassett's K5 spans -0.06 to
+0.46: wrong sign and magnitude at low q, which lifts `D` above where the
measured targets sit. Sec 5 of [MPCE_CPP_PORT_DESIGN.md] establishes that a
negative `K_straight` in dividing flow is real physics; the defect is its
size and where it crosses zero, not its existence.

**Correction (2026-09-05, Finding 6).** This section first read "a monotonic
`D` would leave one operating point", offering that to #272 as a second
motivation. That is wrong: `D` is a quadratic whose vertex lies inside (0,1)
for every unequal-area geometry in the dataset, so **Bassett's own `D` is not
monotone either** and monotonicity is not available to aim at. What fixing the
model does deliver is feasibility everywhere, and the removal of a *spurious*
multiplicity at equal areas. See Finding 6.

## Finding 2: where there are two roots, both are physical

Both branches were checked against the model's own `imposed_q` curve at their
own converged q -- the only single-valued reference available:

| case | seed | q converged | K extracted | model K at that q | consistent |
|---|---|---|---|---|---|
| psi=1, th=120, q=0.406 | default | 0.0553 | 1.0032 | 1.0027 | yes |
| psi=1, th=120, q=0.406 | junction-aware | 0.9393 | 1.8838 | 1.8836 | yes |
| psi=1, th=90, q=0.701 | default | 0.0102 | 0.9986 | 0.9982 | yes |
| psi=1, th=90, q=0.701 | junction-aware | 0.6018 | 0.9806 | 0.9801 | yes |
| psi=3, th=45, q=0.8 | junction-aware | 0.8569 | 3.4375 | 3.4370 | yes |

Every root reproduces the model to four decimals. **There is no non-physical
branch to detect here.** The boundary conditions genuinely admit two operating
points and the initial guess alone decides which is reported. What a detector
can and should catch is the *degenerate* branch: q=0.0102 means a lateral leg
carrying 1% of the flow, which passes `verify_solution_consistent` because its
only criterion is `m_dot > 1e-6`.

## Correction 1: the mirror root was not a mirror root

Recorded on #271 and in [JUNCTION_JACOBIAN_271.md] on 2026-09-05: Bassett
Fig 7b `mfb_two_pb` q=0.8 "converges to a mirror root, K_lat 3.44 against
2.77, reported as a success -- the worst outcome the harness can produce."

**Withdrawn.** The solve converges to q=0.857, not 0.8, and its extracted K of
3.4370 equals the model's own `imposed_q` K at q=0.857 of 3.4370. It is the
model's correct root at a shifted operating point. Against Bassett at the
achieved q (3.3334) the model is 3.1% high, in line with its `imposed_q`
accuracy. The 24% "error" was a comparison against the paper at the target q
while the solve sat somewhere else.

Three statements have now been made about this one case. The first claimed the
verifier demotes it, from a run on a reused adapter instance. The second
claimed a mirror root, from comparing K at the wrong q. Only the third is
measured against the model's own curve. **The lesson is the comparison, not
the care: a pressure-driven root can only be judged against the model at the q
it actually reached.**

## Correction 2: the three_pb K score is a tautology

In `three_pb` every total pressure is imposed by a `PressureBoundary` through a
`LosslessConnectionElement`, and `_extract_K` normalises by the fixed
`m_dot_ref` rather than the converged flow. So

    K_extracted = (Pt_com - Pt_lat)_imposed / q_dyn(m_dot_ref) = K_target

by construction, whatever the model does.

Falsified directly: a deliberately wrong model (`eta_scale=3.0`, whose
`imposed_q` K at q=0.8 moves from 2.88 to 3.07) still reports **2.7677** in
`three_pb` against a target of **2.7689**, a 0.05% "match".

| q | eta_scale=3 imposed_q K | eta_scale=3 three_pb K | target |
|---|---|---|---|
| 0.6 | 1.6623 | 1.2464 | 1.2467 |
| 0.8 | 3.0651 | 2.7677 | 2.7689 |

Once the achieved operating point is instrumented (below) the mechanism is
visible in a single row. Bassett Fig 7b, target q=0.8, target K 2.7689:

| model | topology | K reported | q the solve reached |
|---|---|---|---|
| honest | imposed_q | 2.8842 | 0.800 |
| honest | three_pb | **2.7654** | **0.021** |
| honest | mfb_two_pb | 3.4370 | 0.857 |
| eta_scale=3 | imposed_q | 3.0651 | 0.800 |
| eta_scale=3 | three_pb | **2.7677** | **0.737** |

`three_pb` reports the target to three decimals from a junction passing 2% of
its flow down the lateral leg, and reports it again from a model that is 6%
different at a completely different operating point. It measures convergence
and nothing else. The residual 1e-3 differences are compressible/density
corrections, not model content.

**Acted on:** `test_bassett_fig7b_case.py` asserted this tautology when it was
merged in PR #288 and has been reduced to a convergence test.
`network_runner._K_SCORING_TOPOLOGIES` now excludes `three_pb` from K scoring,
and the scorecard prints `taut` in its RMSE and bias columns rather than a dash,
so the exclusion cannot read as missing data. The perturbation above is pinned
in `python/tests/test_junction_score_at_achieved_q.py` -- if `three_pb` ever
becomes informative, that test fails and says to revisit the exclusion rather
than delete it.

## Finding 3: the operating point drifts, and the score does not notice

Drift between the intended q and the converged q, over every converged
pressure-driven record:

| topology | median drift | 90th pct | drift > 0.25 |
|---|---|---|---|
| `three_pb` | 0.1335 | 0.5961 | 26.2% |
| `mfb_two_pb` | 0.0001 | 0.3899 | 12.6% |

`mfb_two_pb` mostly lands where intended, which is what makes the tail
dangerous: the score compared K against the paper at the target q with no
signal that 12.6% of the rows are at a different operating point altogether.
**Reporting the achieved q is the detector this harness was missing**, and it
subsumes the "wrong root" question: a root 0.4 away in q is a different
operating point, not a failed validation of the intended one.

**Acted on:** `NetworkResult.q_converged` (a field that existed and was never
populated) is now filled by `solve_and_extract`, carried on `NetworkRecord`
with a `q_drift` property, aggregated per cell as `median_q_drift`, and printed
as a `dq` column. It reads back exactly the imposed value under `imposed_q`,
which is the calibration of the instrument. Adapters that re-index a paper's q
axis must invert the same transform on the way out -- Bassett K11 does, and a
test pins it, because getting that wrong is the mirrored-axis defect class of
#276/#277/#279/#283.

The first thing it shows: on `K_lateral_sep` psi=1 theta=45, `three_pb`
reported the *best* RMSE in the row (0.057, against `imposed_q`'s 0.061) while
sitting a median of **0.476** away in q -- half the entire operating range.

**Also acted on (9b): K is now scored on the measured curve at the achieved q**,
by linear interpolation between the curve's own digitised points. Three
decisions, each with a number behind it:

- **Interpolation is bounded, and the bound was measured.** Leave-one-out over
  all 64 curves (drop an interior point, rebuild it from its neighbours):
  median |error| 0.030, p90 0.175, worst 10.75 on a steep Idelchik K12 branch.
  That spans two gaps and so overstates single-gap interpolation, but it is the
  floor below which a rescored K resolves nothing. Against a model RMSE of 0.78
  the median is about 4%, which is why the change is sound; the tail is why
  extrapolation is refused rather than clamped.
- **Outside a curve's digitised range there is no ground truth.** Those records
  are reported as *off-curve* and excluded from RMSE, never extrapolated.
- **Rescoring alone would hide a coverage loss.** A record aimed at q=0.2 that
  lands at q=0.9 is a correct measurement of the model at 0.9 and no evidence
  whatever about 0.2, and rescoring it makes the RMSE look better while the
  requested point goes untested. So *off-point* is counted separately, at a
  threshold of 0.05 -- half the median digitised spacing of the two dominant
  sources (Bassett 0.0997, Idelchik 0.1000).

Before and after, over the K-scoring topologies:

| topology | converged | off-point | off-curve | RMSE at target | RMSE on the curve |
|---|---|---|---|---|---|
| `imposed_q` | 602 | 0 | 0 | 0.7595 | 0.7595 |
| `mfb_two_pb` | 618 | 88 (14%) | 10 | 0.7894 | 0.7783 |

`imposed_q` is untouched to the digit, which is the check that matters: the
split is a boundary condition there, so nothing should move. (An earlier cut
compared floats exactly, and float noise pushed 9 `imposed_q` records off the
end of their own curves and out of the score. A snap at the instrument's
calibration tolerance fixed it, and a test pins it.)

On the 88 records the change actually touches, RMSE falls 0.2270 to 0.1724.
The overall move is small because the median drift is 0.0001: `mfb_two_pb`
usually lands where it was aimed, and it is the tail that was being scored
wrongly.

The coverage numbers are the more interesting output:

| topology | converged | landed elsewhere |
|---|---|---|
| `imposed_q` | 602 | 0 |
| `three_pb` | 480 | **349 (73%)** |
| `mfb_two_pb` | 618 | 88 (14%) |

`three_pb` does not merely report a tautological K. It probes an operating
point other than the one asked for in nearly three quarters of its converged
records.

## Finding 4: the seed defect is real, and worth about 70 solves

`NetworkSolver._propagate_analytical_pt_prop` seeds only `jct.P_jct` for these
networks. Its Bernoulli mass-flow seed is restricted to `ChannelElement`, and
both the harness and `gui/backend/graph_builder.py` wire junctions with
`LosslessConnectionElement`, which gets nothing. Every element therefore falls
back to the flat reference flow (0.1 kg/s), so x0 on a three-port junction
violates continuity outright -- 0.1 in, 0.2 out -- with a 1:1 split whatever
the boundary pressures say.

A junction-aware replacement was prototyped: flow scale from a
`MassFlowBoundary` when present, else Bernoulli on the largest imposed `dPt`;
split proportional to `A_i * sqrt(dPt_i)`; rescaled so continuity holds exactly
at x0 with every port on its declared side.

| topology | default seed | junction-aware |
|---|---|---|
| `imposed_q` | 602 / 691 | 602 / 691 |
| `mfb_two_pb` | 618 / 691 | 635 / 691 |
| `three_pb` | 480 / 691 | 537 / 691 |
| all | 1700 / 2073 | 1774 / 2073 |

`imposed_q` is untouched, as expected: the flows are boundary conditions there.

**It is not answer-neutral.** It moves already-converged roots, because of
Finding 1 -- both roots exist and the seed selects one. That is why it landed
last, after the operating-point reporting and the energy check were in place to
judge the difference. The re-measurement with all of that present is Finding 8.

The baseline also moved between runs (1700 and 1711 on identical inputs).
That is Finding 5 below, not noise in the measurement of the seed: both seeds
were measured in the same process, so the comparison between the two columns
stands.

## Finding 5: the run-to-run variation is PYTHONHASHSEED, and it is one line

Recorded on 2026-09-05 as "process-dependent outcomes, cause unknown". Cause
found and proven the same day.

Bassett Fig 7b `mfb_two_pb` q=0.8, run as ten separate identical processes:

| outcome | wall time | count |
|---|---|---|
| converged at q=0.857, K=3.4370 | 23 ms | 5 |
| degenerate root, rejected by the verifier | 8.6 ms | 5 |

Bimodal, so not a timing-truncated solve (the stall detector's window is 5 s
and these solves take milliseconds). Fixing the interpreter's hash seed makes
it perfectly deterministic, and *which* root you get depends on the seed:

| PYTHONHASHSEED | outcome, 5 runs each |
|---|---|
| 0 | converged, 5/5 |
| 1 | degenerate root, 5/5 |
| 7 | degenerate root, 5/5 |

The dependence enters at `NetworkSolver._propagate_pressure_guess`:

    visited = set(p_guess.keys())
    queue = list(visited)          # <- iteration order of a set of strings

`queue` seeds a breadth-first pressure propagation in which every hop applies
`dp_est`, so the visit order decides the propagated pressures, hence x0, hence
which basin Newton falls into. Set iteration order for strings follows
`hash(str)`, which Python randomises per process by default.

**Fixed, one line:** seed the queue from the dict, which is
insertion-ordered.

    queue = list(p_guess.keys())

Full scorecard under four hash seeds, before and after:

| PYTHONHASHSEED | before: mfb_two_pb | after |
|---|---|---|
| 0 | 618 / 691 | 618 |
| 1 | 629 | 618 |
| 7 | 629 | 618 |
| 13 | 629 | 618 |

(`imposed_q` 602 and `three_pb` 480 throughout; all-topology 1700 after,
1700 or 1711 before.) The fix pins every seed to the seed-0 behaviour, which
is 11 convergences fewer than the lucky seeds. **Fewer is not worse here:**
by Finding 2 both branches satisfy the model, so a convergence count rewards
landing on a degenerate root as much as on the intended one, and seed 0 is
the ordering under which Bassett's own Fig 7b case converges to his operating
point. A reproducible 1700 is worth more than an irreproducible 1711.

It is not junction-specific: any network whose solve is basin-sensitive
inherited it, and a validation harness that counts convergences inherited it
twice. Landed with this record because the Fig 7b tests cannot be
deterministic without it.

## Finding 6: monotonicity is not the right target, and D is a quadratic

Asked whether the coefficients could be reformulated to force monotonicity.
They cannot, and they should not be. Bassett's separating pair has closed
forms (Eq 15 and Eq 27), so their difference does too:

    K5(q)          = q^2 - 1.5 q + 0.5
    K6(q,psi,th)   = q^2 psi^2 + 1 - 2 q psi cos(0.75 th)

    D(q) = K6 - K5 = q^2 (psi^2 - 1) + q (1.5 - 2 psi c) + 0.5,  c = cos(0.75 th)

A quadratic, verified against the functions to 1e-16 over the geometry grid.
Three consequences follow directly from its coefficients.

**(a) At equal areas D is exactly LINEAR.** `psi = 1` kills the quadratic term:
`D = q (1.5 - 2c) + 0.5`. So the operating point is unique for every
equal-area geometry, at every q, with no appeal to numerics.

**(b) At unequal areas D is a cup whose vertex sits inside the range**, so
Bassett's own coefficients are non-monotone there:

| psi | theta | vertex q* | monotone on (0,1)? |
|---|---|---|---|
| 1 | any | linear | yes |
| 2 | 45 | 0.304 | no |
| 3 | 45 | 0.218 | no |
| 4 | 45 | 0.172 | no |
| 3.333 | 90 | 0.052 | no |
| 10 | 90 | 0.031 | no |

Forcing monotonicity would contradict the source. It is not a defect to fix.

**(c) Non-monotone does not mean two operating points.** For a target
`D_t = D(q_t)` the second root is the mirror about the vertex,
`q_2 = 2 q* - q_t`, which is physical only when `0 < 2 q* - q_t < 1`. Counted
over every separating point in the dataset, with Bassett's own coefficients:

| verdict | points | share |
|---|---|---|
| unique -- psi = 1, D linear | 86 | 81.9% |
| unique -- mirror outside [0,1] | 9 | 8.6% |
| genuinely two roots | 10 | 9.5% |

**90.5% of the dataset has a unique physical operating point.** The 10 that do
not are all psi=3, theta=45 at q_t < 0.436.

### The model's multiplicity is in the wrong place

Every two-root case measured in Finding 2 was at **psi = 1**, where the physics
is provably unique. The model manufactures an extremum where Bassett has a
straight line:

| psi=1 | Bassett D | model D |
|---|---|---|
| theta=45 | linear, slope -0.163 | monotone, slope -0.7 to -2.9, goes negative |
| theta=90 | linear, slope +0.735 | hump, turning at q=0.30 |
| theta=120 | linear, slope +1.500 | hump, turning at q=0.50 |

Mynard's energy-transfer factor is implicated but is not the whole story:
`eta_scale=0` removes the hump at theta=90, and at theta=120 it leaves `D`
essentially flat (1.000 to 0.994) where Bassett rises from 0.5 to 1.925. The
equal-area defect is structural, not a mis-tuned constant.

### The reformulation that IS available, and it is exact

Not "force monotonicity" but **"reproduce the analytical identity the source
already gives you"**:

    at psi = 1:   K_lateral(q) - K_straight(q) = q (1.5 - 2 cos(0.75 theta)) + 0.5

Linear in q, with a slope and intercept fixed by geometry alone. That is a
sharper acceptance criterion for #272 than any error metric: a structural
identity the model must satisfy exactly in the equal-area limit, currently
violated in both slope and shape. It also disposes of the spurious multiplicity
by construction, since a linear D cannot have two roots.

### What selects the root where multiplicity is genuine

Not monotonicity. Stability. Quasi-steady momentum on the two legs, with the
split as the only freedom (m_com fixed, both outlet Pt fixed):

    2 L d(m_lat)/dt = (Pt_str - Pt_bra) - D(q) Q

Perturbing about a root gives `d(dq)/dt = -Q D'(q) dq / (2 L m_com)`, so

    a root is dynamically stable iff D'(q) > 0.

Checked on every two-root case measured: exactly one root has `D' > 0`.

| psi | theta | q_t | roots | D' at each | stable |
|---|---|---|---|---|---|
| 1 | 120 | 0.406 | 0.055, 0.935 | +1.87, -2.08 | 0.055 |
| 1 | 120 | 0.490 | 0.135, 0.865 | +1.54, -1.51 | 0.135 |
| 1 | 120 | 0.637 | 0.325, 0.675 | +0.75, -0.70 | 0.325 |

`D'` is already available analytically inside the element -- it is what the
Jacobian carries -- so this costs nothing to evaluate. It is the principled
form of defect 10's "reject roots that are not physical", and it is a stronger
criterion than the present minimum-flow test because it rests on a stability
argument rather than a threshold.

**Caveat, stated because it bounds the claim.** The argument is quasi-steady
and one-dimensional, and it says a pressure-driven tee cannot sit on the
falling branch. Bassett measured those low-q points in a flow-controlled rig,
where the split is imposed and stability does not select, so there is no
conflict with the data. The criterion applies to pressure-driven networks, not
to the validation of the coefficients themselves.

## Finding 7: v2 had no energy check, and the model creates work at low q

Defect 10 was recorded as "`verify_solution_consistent` passes near-degenerate
roots (a lateral leg at 1% of the flow); its only criterion is
`m_dot > 1e-6`". **That framing was wrong.** A tee CAN run at 1% lateral flow
if the pressures say so; landing there when a different point was requested is
a coverage question, and item 9 already measures it. Nothing about a small
flow is unphysical.

The real gap is that **v1 has an energy check and v2 had none**. With the
element's own residual convention the dissipation across the junction is
`q_dyn_com * sum_i |m_i| K_i`, so requiring a passive junction not to create
flow work is

    sum_in |m_i| Pt_i - sum_out |m_i| Pt_i >= 0

computable from Pt and m_dot alone, with no closure re-evaluation.

**Why not v1's form.** v1 bounds every collector's Pt by the single supplier's,
which forbids *any* branch gaining total pressure. That is real physics in
dividing flow -- Bassett's K5 and Hager's xi_t both go negative -- and v1 only
escapes rejecting it because its relative tolerance (1e-4 of ~1e5 Pa, so ~10
Pa) is wider than the effect at the validation dynamic head (~2.6 Pa). It is
also single-supplier only, so it skips joining entirely. The mass-weighted
balance permits one branch to gain at another's expense, forbids a net gain,
and covers both flow directions.

### Audited before being imposed

| source | worst mass-weighted mean K | violations |
|---|---|---|
| Bassett analytical pair, correctly paired | **exactly 0** in the no-diversion limit | 0 |
| **the model** | **-0.06** | violates for q below ~0.2 |

**Corrected 2026-09-05.** The first version of this table read "measured curves
+0.21, Bassett analytical +0.19", computed by pairing K5 and K6 at the same q.
That is the wrong pairing: Bassett indexes each coefficient on its own leg, so
K5's argument is the STRAIGHT fraction (Finding 9). Paired correctly the
quantity is `(1-q) K5(1-q) + q K6(q)`, which is positive throughout and tends
to **exactly zero** as the lateral flow vanishes -- a junction that diverts
nothing dissipates nothing, which is the sharper statement and a good check on
the axis reading. The mirrored pairing gives 0.5 there, a loss with no flow
diverted. The conclusion is unchanged and strengthened; the two margin figures
were wrong.

The data satisfies the condition with margin. The model does not: it is
net-generating at low lateral fraction in four of five geometries checked
(psi=1/th=90, psi=1/th=120, psi=3/th=45, psi=3.333/th=90), with the violating
band running from q=0.02 to about q=0.22. That is the K_straight gap of #272
appearing as a hard physical violation rather than an accuracy shortfall,
independent of the equal-area identity in Finding 6.

### Cost, measured

| topology | converged before | after | demoted |
|---|---|---|---|
| `imposed_q` | 602 | 590 | 12 |
| `three_pb` | 480 | 428 | 52 |
| `mfb_two_pb` | 618 | 607 | 11 |
| all | 1700 | 1625 | 75 |

The audit predicted 74 from the converged states alone, so the outcome is the
one that was measured in advance. Four previously-green assertions in
`test_bassett_fig7b_case.py` now fail and have become strict xfail targets:
the model cannot be evaluated at psi=3, theta=45, q=0.2 because it is
inadmissible there, and `three_pb` at q=0.6 and 0.8 wanders down into the same
band.

**The honest reading of the cost.** These 75 were not solver failures being
created; they were physically impossible answers being reported as successes.
Losing them is the gate working. But it does mean the scorecard can no longer
score the model in its inadmissible region, so the region is now visible as
non-convergence rather than as a bad K -- and the solver's message was widened
to say that a non-dissipative closure state is one of the causes, so it is not
misread as a solver artifact.

### A caveat that bounds the check

This is the model's own dissipation identity, not a general entropy statement.
Where port total temperatures differ, true entropy production also carries a
positive mixing term that is not counted, so on a strongly non-isothermal
junction the check is conservative in the wrong direction. The element is
documented for low Mach and the fixtures are isothermal, so the case does not
arise today; it is recorded at the method.

## Finding 8: what the seed is actually worth, measured with the instruments in place

Re-measured after items 8, 9 and 10 landed, which changed the baseline from
1700 to 1625 and gave the harness three things the first prototype lacked: the
achieved operating point, the off-point count, and the energy gate.

| | seed off | seed on |
|---|---|---|
| converged | 1625 | **1732** |
| rejected as inadmissible or artifact | 227 | **163** |
| on-point evidence | 1243 | 1254 |
| RMSE on the common subset, `imposed_q` | 0.7664 | 0.7664 |
| RMSE on the common subset, `mfb_two_pb` | 0.7845 | 0.7844 |

**Read the middle row first.** The seed's real effect is not "+107 rows": it is
that 89 solves which previously landed on an inadmissible or artifact root now
land on an admissible one, against 20 pushed the other way. That is the
"steer away from non-physical branches" objective, and it is only visible
because the energy check of item 10 exists to name those roots.

**The convergence count on its own would have flattered it.** Of the 153 rows
gained, 21 are on-point and 131 are not; 47 rows are lost, 11 of them on-point.
Net evidence at the requested operating points is +11 on 2073 records. Anyone
quoting the +107 as the benefit is quoting mostly off-point rows.

**Accuracy is unchanged.** On the records scored under both seeds, RMSE and
bias agree to four decimals. A first pass compared full populations and showed
`mfb_two_pb` RMSE rising 0.7806 to 0.7884, which was a population effect from
the 36 extra rows, not a degradation -- the common-subset comparison is the
one that answers the question.

**One earlier claim is withdrawn.** The prototype note said the seed "selects
the intended operating point" on at least one case, generalising from six
moved roots of which five landed nearer the target. Measured across every
record converged under both seeds, it is a wash: `three_pb` median drift 0.0949
to 0.0944 with 167 nearer and 195 further; `mfb_two_pb` 0.0001 to 0.0001 with
231 nearer and 205 further. **The seed does not steer toward the requested q.**
It steers toward admissible roots, which is a different and better thing.

### Two things it broke, and what they turned out to mean

**The ejector.** `EjectorElement` extends `MultiPortChamberElement`, so a
blanket junction seed reached it and overwrote the physics-based warm start it
carries for its own port pressures and `P_jct`. The cold reference network in
`test_gui_ejector.py` stopped converging. The seed is now opt-in through a
class attribute, `seeds_ports_by_pressure_split`, true on the chamber junction
and false on the ejector, because an ejector's port flows follow choking and
entrainment rather than a Bernoulli share of the imposed pressure difference.

**A bistable production fixture.** `_three_port_net` in
`test_multi_port_chamber.py` has three pressure boundaries at 2.10, 2.05 and
2.00 bar with friction in every channel. Two flow arrangements satisfy it:

| mode | common port | straight port | branch port |
|---|---|---|---|
| dividing (seed on) | in, 0.450 | **out**, 0.379 | out, 0.071 |
| ejector (seed off) | in, 0.358 | **in**, 0.211 | out, 0.569 |

Both conserve mass, both are net dissipative, both pass the v1 energy check,
and both are consistent with the boundary pressures -- in the second the 2.05
bar boundary also feeds the junction. The seed reaches the first directly; the
plain cold start reaches the mirror root, is demoted, and the auto-retry
reaches the second. **Nothing available says which is "the" answer**, so the
bistability is now pinned as a test rather than resolved, and the two solver
machinery tests that depended on the old path force it explicitly instead of
being weakened.

The change is deterministic: 1732 under hash seeds 0, 1 and 7.

## Finding 9: the separating straight leg was scored on a mirrored axis

Found while re-reading Mynard and Bassett to answer a question about whether a
negative K implies an energy gain (it does not; see Finding 7's correction).

Bassett indexes each coefficient on the fraction in **its own leg**. K6's
argument is the lateral fraction; K5's and K2's is the **straight** fraction.
Every network adapter builds its network from the lateral fraction, so a
straight-leg record needs `1 - q` on the way in, and the achieved operating
point needs the same transform coming back.

Four independent lines settle the convention:

1. `validation/junction/equivalences.py` states it outright: "Hager q is
   lateral/total -> Bassett K5 q is straight/total = 1 - q_Hager".
2. The algebra: `K5(q) = q^2 - 1.5q + 0.5` is exactly Hager's `xi_t` evaluated
   at `1 - q`, and Hager's argument is the lateral fraction.
3. The digitised points fit Bassett's own K5 curve read directly (mean error
   0.037) far better than mirrored (0.240), across all four measured files.
   Hager's `xi_t` files fit his own curve on the lateral axis (0.08 against
   0.21), as documented.
4. The physical limit. Paired correctly, the flow-weighted mean K tends to
   **exactly zero** as the lateral flow vanishes: a junction that diverts
   nothing dissipates nothing. Paired the other way it tends to 0.5, a loss
   with no flow diverted.

### Audited at every site

| adapter | state |
|---|---|
| `mpce_v1_network` | input transform present (added with #277) |
| `mpce_v2_network` | **mirrored** |
| `tee_junction_element_network` | **mirrored** |
| `mynard_analytical` | **mirrored** |
| `tee_junction_raw` | correct: it calls Bassett's own formula on Bassett's own axis |

**None of the four did the output half**, which only became reachable when
`q_converged` was populated in #290 -- so the drift and the on-curve
interpolation added for item 9b were being measured on a mirrored axis for
every straight-leg record.

### Effect, measured

| topology | converged before | after | RMSE before | after |
|---|---|---|---|---|
| `imposed_q` | 77 | 79 | 0.3770 | **0.3178** |
| `mfb_two_pb` | 70 | 70 | 0.1962 | **0.1399** |

Every other canonical coefficient is bit-identical -- lateral separating,
lateral joining and straight joining all unchanged -- which is the check that
the fix is surgical. Overall convergence 1732 to 1734.

**One honest caveat.** For MPCE-v1 the correct axis makes the straight-leg
score *worse* (0.477 mirrored against 0.551 correct). That is not an argument
against the fix: the axis is fixed by the paper's convention, not by which
reading fits better. It is a statement about v1's straight-leg physics, which
is separately known to be poor because its collinear ports are collapsed by
the sin^2 gate.

## What this changes

- The 26 unrescued solves are no longer a mystery: most are infeasible, and
  the infeasibility measures the `K_straight` gap.
- #272 gains a sharper acceptance criterion than any error metric: the
  equal-area identity in Finding 6, which the model currently violates in both
  slope and shape.
- The step-3 C++ port gate must not be "reproduce the `three_pb` K table" --
  that table is reproducible by any model at all.
- `verify_solution_consistent` needs a minimum-flow criterion, not an
  energy-consistency one, for the branch class actually observed here.
- The non-determinism is not an accepted cost of a hard problem. It is one
  line, and until it is fixed no convergence count in this harness is
  reproducible.

## Dead ends

- **Seeding the physical branch by hand.** A hand-built physical seed
  (continuity-consistent, correct signs, `P_jct` at the common total pressure)
  does not rescue the q=0.2 or q=0.4 Fig 7b cases. It cannot: there is no root.
  This was the experiment that turned the investigation from seeding to
  feasibility.
- **Raising the dynamic head.** The imposed differences at q=0.2 are ~10 Pa on
  a 1e5 Pa level, so ill-conditioning was the natural second hypothesis.
  Sweeping `m_in` from 0.1 to 1.6 kg/s (a 256x change in dynamic head) left
  every failing case failing. Not conditioning.
- **Treating the second root as an artifact.** It reproduces the model exactly.
  Rejecting it would have meant writing a detector against a correct answer.

## Tests and code that pin this

- `python/tests/test_bassett_fig7b_case.py` -- the passing points, the
  infeasible ones as documented xfail targets, and no tautological assertion.
- `validation/junction/models/_network_builder.py` -- `_extract_K` and the
  three skeletons; the tautology lives in the combination of the two.
- `python/combaero/network/solver.py` -- `_propagate_analytical_pt_prop`
  (the seed) and `verify_solution_consistent`'s call site (the detector).

[MPCE_CPP_PORT_DESIGN.md]: ../../validation/junction/MPCE_CPP_PORT_DESIGN.md
[JUNCTION_JACOBIAN_271.md]: JUNCTION_JACOBIAN_271.md
[issue #271]: https://github.com/thiemom/combaero/issues/271
