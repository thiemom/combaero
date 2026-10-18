- **A failed solve says when a declared-direction junction is the reason.**
  `MultiPortChamberElement` and `ConstantKTeeElement` take a `flow_direction`
  and do not reverse a port. When a solve fails with such a junction holding
  a port at zero or against its declared direction, the message now names
  the junction and the port, and points at the boundary pressures,
  `flow_direction` or `IdealMomentumTeeElement`. Every one of the 8 stress
  `tee3` draws that stalled has a root with the ideal tee, with that port
  reversed; before, the message said only "not making good progress".
