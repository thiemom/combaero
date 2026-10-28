- **Impingement jet holes take a crossflow-dependent Cd by default (#800).** A strong crossflow along the gap bends the jets over and lowers their discharge coefficient. `ImpingementPlateElement` and `ImpingementArray` gain `cd_crossflow`, also exposed in the GUI:
  - `'rogers_hersh_1976'` (the new default): Rogers & Hersh (1976, NASA CR-2681) measured single-hole fits. It applies to any plate and stays capped at the hole's own Cd below Gc/Gj ~1.5.
  - `'florschuetz_metzger_su_1982'`: CR-3630 Table 4.1, the authors' own law for the plates they tested, and the one behind Table D.1's flows.
  - `'none'`: the previous constant Cd.

  C++ `rogers_hersh_1976_cd` and `florschuetz_metzger_su_1982_cd` return the Cd with analytic derivatives. The law is written in terms of the ideal jet velocity, so the hole's flow goes as dP^1.07 under a strong crossflow instead of sqrt(dP).
  - **Unchanged:** every Florschuetz (1981) array is bit-identical, since its Gc/Gj stays at or below 0.77, so the Eq. 8 flow model and the Nu chain (MAE 6.88%) do not move.
  - **CR-3630 Table D.1 (cross-source law):**
    - The stagnant upstream rows of B(5,4,1) at mc/mj 0.99 no longer alternate in sign (#799); they now rise monotonically, as the data do.
    - Worst-row Gj/Gj_mean improves: z/d >= 2 max 0.111 -> 0.085, z/d = 1 max 0.269 -> 0.231.
    - P/P0 mean error improves 0.0075 -> 0.0073.
  - **CR-3630 Table D.2:**
    - eta is unchanged (bias +0.013, RMSE 0.077).
    - Nu is unchanged inside the correlation's range.
    - Starved rows go from 8 to 11, with Nu at -81% (-95% before): there is still no channel-flow floor (#467).
  - **Solver:** a cold solve that fails on such a network is retried once from the constant-Cd root. CR-3630 D.2 B(5,4,1) at mc/mj 1.01 with a heated bypass needs it.
