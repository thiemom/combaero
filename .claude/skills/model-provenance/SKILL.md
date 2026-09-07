---
name: model-provenance
description: Ground solver elements and model changes in their literature/physics provenance, verify them independently, and close them into a durable rationale record. Use when implementing a new element or correlation from a paper, deriving or updating an analytical Jacobian or closed-form, accepting an agent-implemented solver change, or closing out a completed model/feature.
---

# Model Provenance

Working code can pass every test while resting on a constant nobody derived.
Provenance is the discipline of knowing, for every element and every
regression test, where it came from and why — before it's trusted.

## Workflow

1. Before wiring a new element or correlation into the network solver,
   reproduce the source's own reported case in isolation — same inputs, same
   tolerance — before integration. No reported case in the source? Use the
   nearest analytical limit as ground truth instead.
2. If the change introduces an analytical derivative (Jacobian, closed-form
   gradient), verify it against a finite-difference or other independent
   check. Never let two code paths lean on the same derivation without a
   second one confirming it.
3. Before accepting any agent-implemented change, ask directly: which
   constants, tolerances, or fallback behaviors trace to "the test passes"
   rather than to physics or the cited literature — and separately, which
   choices were invented to fill a gap the spec left open, rather than
   forced by it? An unanswered or fudge-shaped constant is unresolved, not
   done.
4. Write regression tests against the physics/literature ground truth, not
   the code's current output. Deterministic invariants (conservation, known
   limiting cases) get exact tolerances; tuning-dependent behavior (damping,
   pseudo-transient continuation parameters) gets a labeled band. Falsify
   each new regression test once — revert the fix, confirm red for the right
   reason, restore. Then measure coverage separately; see below for why one
   does not substitute for the other.
5. On completion, close the change into a short provenance record: the
   formulation chosen and the alternatives considered, the literature basis
   (source plus the specific case or correlation), the invariant it must
   preserve, dead ends tried, and a pointer to the test that pins it. WHY
   only — never restate the code.

## Falsification and coverage are different checks

**Falsification checks that the tests bind to the code you wrote; coverage
checks whether you wrote tests for the code you actually have.** Run both.
Neither finds what the other does.

Falsification perturbs the implementation and asserts a test goes red. It
catches a test that passes for the wrong reason — a tautology, a vacuous
assertion, a fixture where every quantity is imposed. Its most useful outcome
is the perturbation that does **not** go red: that names a behaviour nothing
tests. Porting `MPCEv2Element` to C++, deleting the line that keeps the
snapping out of the mass-conservation row changed no test result. The line was
right and the gap was real, so the semantic got its own test.

Coverage asks a different question, and templates make it lie. `dual_number.h`
reported 100% of regions, lines and branches while five of its thirteen
operator overloads — including `operator+(D, D)` — were never called. An
uninstantiated template is never emitted, so it cannot be counted as missed.
**Check instantiation counts, not percentages.** Across four port PRs, coverage
found a real gap in three of them (an angle below `-pi`, the joining term with
a single supplier, a `1e-9` dead band a comment made a claim about) and seven
to eight falsifications per PR found none of the three.

Record which branches are genuinely unreachable and why, so the next reader
does not hunt for them.

## Reward hacking is failure, not a shortcut

A number that improved because the model got better is a result. A number that
improved because the measurement was bent is a defect that looks like a
result, and it is worse than no change at all.

- **Never tune a constant against the score it is judged by.** Empirical terms
  prove themselves on the digitised data and are labelled as tuned. A CFD
  correction is tuning; it does not reproduce the source CFD.
- **Prefer a saturation plateau to a peak.** If the metric is monotone in a
  parameter and flat above some value, there is no peak to fit and the choice
  is safe. If it peaks, ask what the peak is made of before taking it.
- **Falsify the metric, not only the change** — perturb the model and confirm
  the score moves. A cell whose every quantity is imposed scores nothing.
- **Never falsify by removing one half of a matched pair**, and never widen a
  tolerance or delete a case to make a suite green.
- Improving one case while the aggregate stands still is not an improvement;
  say so and revert. A negative result, reported, is a good outcome.

## Rules

- Comments and names carry the WHY — which correction, which paper, why this
  over the alternative. Narrative of what was tried and abandoned goes in
  the provenance record, not inline.
- No lineage names (`_v2`, `_final`, `_momentum_cv_new`) — git history and
  the provenance record hold the past; identifiers describe the present.
- A doc points at the file that owns a fact (`materials.h`, an element's
  source file) rather than transcribing its current contents.

## Done

Every new element or correlation has a spike verified against its source
before integration. Every analytical derivative has an independent
cross-check. Every new test has been falsified once, and coverage has been
measured separately with instantiation counts checked. Every accepted agent
change has an answered provenance question. Every number that moved can be
traced to a change in the model rather than to a change in the measurement.
Every closed model has a short record you could read in six months without
this conversation.
