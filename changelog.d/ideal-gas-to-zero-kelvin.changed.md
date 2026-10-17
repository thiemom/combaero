- **Species thermo below 200 K is an ideal gas down to 0 K.** Below the
  NASA-9 data, cp/R now fades over the whole [0, 200 K] (C2 at 200 K) to the
  rigid-rotor value, 5/2 for an atom (and H2, whose rotation freezes near
  85 K), 7/2 linear, 4 nonlinear, instead of levelling off at its 200 K value
  50 K below. h(0 K) is finite; `calc_T_from_h` resolves T down to a few K and
  throws below h(0 K). New constants `thermo::CP_R_RIGID_ATOM/_LINEAR/_NONLINEAR`;
  `thermo::TAIL_FRACTION_LOW` is removed.
- **No [200, 5000] K clamp on network node temperatures.** A node carries its
  T(h) as is: heated past 5000 K it is simply hot, energy conserved (the clamp
  discarded ~190 kW of 200 in a heated plenum). `__T_clamped__` now lists the
  nodes with flow at or below 0 K (more heat out than the stream holds as a
  gas), and such a node reports `{node}.no_gas_state` instead of properties.
  The mixer inverts its enthalpy exactly to 0 K (the 134 K linear
  continuation is gone) and continues through it on h's line.
- **A stagnant node takes up no heat.** Below `MIXER_HEAT_MDOT_FLOOR` (1 mg/s)
  the heat a node takes up fades C1 to none at zero flow
  (`mixer_heat_per_unit_flow`, new, with its derivative), so a dead end holds
  the temperature of what reaches it; the previous floor put a 100 W heater in
  a dead end at 1e5 K. A root inside the fade -- a heated node neither flowing
  nor stagnant -- is refused as INCONSISTENT and listed in `__heat_faded__`.
