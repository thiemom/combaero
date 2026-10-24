- `AreaChangeElement` diagnostics now report the loss the network solved:
  `zeta` and `dP_static` no longer apply the kernels' `k1 = 1 + 0.35 M^2`,
  which the residual never uses (they read 3.6-5.0% high at a small-section
  Mach of 0.33-0.39). The reference `velocity` and `Re` are now taken at the
  smaller section, as labelled; a contraction used to report them at its
  larger inlet face (28 m/s against 114 m/s).
