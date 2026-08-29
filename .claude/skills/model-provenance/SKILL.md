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
   reason, restore.
5. On completion, close the change into a short provenance record: the
   formulation chosen and the alternatives considered, the literature basis
   (source plus the specific case or correlation), the invariant it must
   preserve, dead ends tried, and a pointer to the test that pins it. WHY
   only — never restate the code.

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
cross-check. Every accepted agent change has an answered provenance
question. Every closed model has a short record you could read in six
months without this conversation.
