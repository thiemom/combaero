# Flow-Area Inference: Provenance Record (issue #262)

Why the area-inference rules are what they are. The contract itself lives in
[API_PYTHON.md](../API_PYTHON.md) ("Flow-Area Inference"); the behaviour is
pinned by `python/tests/test_area_inference.py`.

## The defect

Seven `resolve_topology` implementations searched for a neighbouring
`ChannelElement` and, when none was found, assigned a hard-coded constant
(0.01 / 0.02 / 0.1 m^2). Those constants trace to "something had to be there",
not to any geometry -- the fudge-shaped-constant case the model-provenance
discipline exists to catch. With a non-channel neighbour (a momentum chamber,
a combustor, an ejector outlet) the fallback fabricated an area that the
network does not contain.

## What the audit established

Two measurements decided the design; neither was predictable from reading the
code.

**1. Which sites are genuinely affected.** Probing all seven with a
non-channel, area-bearing neighbour: six resolved to their hard-coded default
and ignored the real neighbour area; `BorderCarnotLossElement` refused to
guess and raised, which is where the "raise, don't default" precedent came
from.

**2. Which defaults are load-bearing.** Counting fallback hits across the test
suite, then perturbing the constant 0.1 -> 0.2 and re-running:

| site | suite hits | of which a neighbour area existed | sensitive to the value |
|---|---|---|---|
| `AreaChangeElement.F0` / `.F1` | 0 | - | - |
| `TeeJunctionElement.F_C` | 0 | - | - |
| `ChannelElement.diameter` | 0 | - | - |
| `MomentumChamberNode.area` | 20 | 2 | no |
| `CombustorNode.area` | 29 | 0 | no |
| `PressureLossElement.area` | 13 | 0 | 3 tests |

The three sites the issue is actually about are never exercised, so making
them raise costs nothing. The three that fire constantly are insensitive to
the value in 59 of 62 cases -- in 60 of 62 no neighbour area exists anywhere,
so inference cannot rescue them and a raise would force a migration on
networks that are working correctly.

## Formulation chosen

Split by what the area *means*, not by node-versus-element:

- **The area is the geometry being modelled** (`AreaChangeElement` F0/F1,
  `ChannelElement` diameter, `TeeJunctionElement` F_C,
  `BorderCarnotLossElement` area) -- raise. A default here invents a
  contraction or expansion, and no value is defensible.
- **The area is a dynamic-head reference** (`MomentumChamberNode`,
  `CombustorNode`, `PressureLossElement`) -- default, but warn when a
  correlation will actually consume it. `PressureLossElement` warns only for
  head-loss correlations (which read `area`), not fraction-based ones (which
  do not). Across 1891 tests this fires on exactly the 3 the perturbation
  audit flagged as sensitive: no false positives, no misses.

Channels keep first priority in the search so that every network which
resolved before resolves to the same area now.

## Invariant

An area may be inferred only from a *known* area -- one the user set, or one
already resolved from a known area. An auto-sized area is never a source;
inheriting a placeholder launders one guess into another. Pinned by
`test_auto_sized_neighbour_is_not_used_as_a_source`.

## Alternatives considered

- **Raise at all seven sites.** Rejected on the measurement above: it breaks
  60 in-suite cases where nothing can be inferred and the value is never read.
- **Warn at all six defaulting sites.** Rejected as crying wolf -- 59 of 62
  hits are inert, and a warning that fires on correct networks trains users to
  ignore it.
- **Mirror the GUI's node-area lookup into the library.** The lookup in
  `gui/backend/graph_builder.py` covered `discrete_loss` only and made GUI and
  headless networks resolve different areas. Superseded rather than copied:
  the shared path now covers it and the GUI-side special case was removed.

## Dead end

The reported symptom was a solver stall, and the stall could not be reproduced
in a minimal fixture -- it needed the full ejector network. Mass-flow-driven
and pressure-driven chains both converge with the fabricated area. That turned
out to matter: in a plain chain the fabricated throat is *worse* than a stall,
converging confidently to a wrong answer (chamber stagnation pressure inflated
by 4% at 1 kg/s, 64% at 5 kg/s, 342% at 20 kg/s). The regression test pins
that pressure against the physical dynamic head rather than asserting
convergence, which the pre-fix code also satisfied.
