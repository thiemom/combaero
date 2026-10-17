- **An orifice whose Cd vanishes no longer raises.** A crossflow-fed hole's
  McGreehan-Schotsch Cd underflows to exactly 0 where its pressure drop
  vanishes, and the nozzle kernels refused `A_eff <= 0`, so the solver's
  safety net rejected those iterates (35 caught probes across the stress
  effusion family). Both kernels are exactly linear in Cd, so a Cd below
  1e-12 is now evaluated at Cd = 1 and scaled: a trace flow, bit for bit,
  instead of an exception.
- **A trace species no longer divides by zero in a pressure loss's
  Jacobian.** The composition column took a central step of half the mass
  fraction, which is 0 at a denormal fraction (5e-324). Below 2e-6 it now
  takes a forward step of 1e-6.
