- **The metering orifice's Jacobian carries its Cd's flow dependence.**
  `ReaderHarrisGallagher` (the `OrificeElement` default), `Stolz` and
  `Miller` read the pipe Reynolds number Re_D = 4|m|/(pi D mu), so the row
  `m - Cd(Re_D(m)) K` has slope `1 - K dCd/dm`. The Jacobian gave 1: about
  1% off for a 2 cm orifice at 1600 K, in both regimes. The documented gap is
  closed: the C++ `Cd_sharp_thin_plate_and_dRe` (new, also in Python) gives
  the Reader-Harris/Gallagher derivative analytically, and Stolz and Miller
  get their closed-form derivatives. Results do not change; convergence
  paths can.
