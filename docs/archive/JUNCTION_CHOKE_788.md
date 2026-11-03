# Momentum chamber past its sonic point (#788)

Provenance record. The code lives in `src/solver_interface.cpp`
(`momentum_chamber_residual_and_jacobian`) and
`NetworkSolver._choked_bounded_junction` (`python/combaero/network/solver.py`).
The tests are in `python/tests/test_junction_choke_788.py`.

## Problem

At fixed (T0, G = |m|/A) the compressible chamber closure Pt_cl(P) reaches its
minimum Pt* at the sonic static pressure P*. Below Pt* it has no root, and
above Pt* it has two: one subsonic and one supersonic.

- **Two roots:** Newton could land on the supersonic one. The issue's 95 mm
  junction "converged" at Mach 1.46.
- **No root:** a choked channel discharging into a junction of its own area,
  when the element downstream needed a Pt below Pt*. The junction closure had
  no root, so the solve stalled.

## Formulation chosen

Below P*, the closure is replaced by its mirror image about Pt*:
g(P) = 2 Pt* - Pt_cl(P).

- **Monotone.** g rises monotonically through P* and is C1 there, with both
  slopes zero. Every Pt therefore has exactly one root. That root is subsonic
  whenever Pt >= Pt*, so branch selection is a property of the residual, not
  of the seed.
- **No new partials needed.** The partials of Pt* need no dP*/dx, because
  dPt_cl/dP = 0 at P* (the envelope theorem). Subsonic evaluations are
  unchanged bit for bit.

A root on the continuation means Pt < Pt*: the chamber's area cannot carry
its flux at that Pt. The solver accepts such a root as a **choked-bounded
junction** only when both of these hold:

- One choked plain compressible channel of no larger area feeds the
  junction. The junction is then that channel's under-expanded exit, the
  same situation a pressure boundary handles as "choked bounded" (#636).
- No row other than the junction's own reads its static pressure. This is
  checked on the Jacobian column. When it holds, the junction's P is a slave
  variable, and the other rows of the network are exactly those of a
  `PlenumNode` in its place.

The root is then the plenum's. The measured cases match the plenum network to
1e-7. The channel's Fanno flux saturates, as it does at a boundary, so the
combustor pressure equals that of the channel discharging straight into the
outlet.

**Conservation.** Mass and Tt pass the junction unchanged. The expansion loss
Pt* - Pt falls inside the channel's control volume, where sigma =
m R ln(Pt_in/Pt) > 0.

**Refusals.** Any other root past the sonic point is refused (INCONSISTENT),
with the reason stated in the message.

**Reporting.** The junction reports as the plenum it is equivalent to:
P = Pt, at rest. `P_sonic`, `Pt_sonic` and `expansion_loss` are reported
alongside.

## Alternatives considered

- **Linear continuation `Pt* - (P* - P)`, or the sonic state at the
  junction's own Pt (P = Pt P*/Pt*).** Both make the unknown look physical.
  Both are C0 only at P*, and the first goes negative at deep bounding. The
  "sonic at own Pt" report also put the junction's static P below the outlet
  pressure (static P rising across the orifice), the very symptom the issue
  flagged. Rejected.
- **Reporting the channel exit static P*.** Above Pt in deep cases
  (bore = pipe). Rejected.
- **A complementarity row in the downstream element** (orifice-controlled or
  pipe-controlled). Correct, but a redesign of every downstream element.
  Unnecessary, because the network's rows already are the plenum's once the
  junction P is a slave.

## Gaps left (measured)

- **Bore == pipe (pass-through), cold solve.** The default x0 put the junction P past the fold, and the first step overshot to a negative P. This is fixed by `NetworkSolver._seed_chambers_past_sonic`, which moves only a start past the sonic point onto the closure's root. Its probe propagation is undone, because lagged reads would otherwise make the first residual depend on it.
- **In-pipe orifice downstream of a choked-bounded junction.** The orifice
  reads its approach from the junction Pt at the pipe's flux, and past the
  pipe's own sonic limit its row's flow columns have a Jacobian error
  (FD 4.4 vs 7.4), issue #958. The same error appears with a `PlenumNode` junction and
  on main. It belongs to the orifice row, not to this change.
