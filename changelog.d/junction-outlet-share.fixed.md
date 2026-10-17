- **A junction with several inlets and several outlets feeds each outlet its
  own flow.** The solver gave every sink of a multi-source element each
  source's own flow. That is right for a merge, which has one outlet, but at a
  `MultiPortChamberElement` or `IdealMomentumTeeElement` with two or more
  inlets AND two or more outlets, each outlet chamber took its velocity head
  at the junction's whole throughput.
  - **Converged results were wrong:** in the stress `ideal_nport` draws, port
    flows moved by up to 2x on fully forward roots.
  - **Solves stalled:** a port flow crossing zero jumped its chamber closure
    by about 10 kPa, which stalled 13 of 40 draws. All 40 now converge.
  - **The GUI's 3-port tees** (one inlet or one outlet) were not affected.
