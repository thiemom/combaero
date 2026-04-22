feat(solver): add sharp-edge and conical area change elements with analytical Jacobians

- Add sharp_area_change() physics function (area_change.h/cpp) with:
  * Borda-Carnot expansion loss: zeta = alpha - 2*ar + ar², alpha(Re)
  * Weisbach contraction loss: zeta = 0.5*(1-ar)
  * Sigmoid direction blending for smooth bidirectional flow
  * smooth_abs(m) regularization for gradient continuity at m=0
  * Analytical Jacobians: dS_dm, dS_drho, dS_dmu (including Re→alpha chain)
  * Optional D_h parameter for non-circular ducts (decouples hydraulic diameter from area)
  * Mach correction k1 = 1 + 0.35*M² with MACH_CLAMP=1.5 and MACH_ACCURACY_LIMIT=0.9
  * mach_clamped diagnostic flag when |M| > MACH_CLAMP

- Add conical_area_change() physics function with:
  * Diffuser loss: zeta = eta(theta)*(1-ar)², eta = 1 - exp(-3.2*theta^0.8)
  * Nozzle loss: zeta = C_con(theta)*(1-ar), C_con = 0.04 + 0.46*sin(theta)^1.5
  * Half-angle theta derived from F0, F1, and axial length
  * Purely geometric coefficients (no Re dependence, dS_dmu = 0)
  * Same Mach correction, mach_clamped flag, and smooth_abs as sharp

- Add solver_interface wrappers:
  * area_change_residuals_and_jacobian() for sharp-edge
  * conical_area_change_residuals_and_jacobian() for gradual
  * Evaluates rho, mu from (T, P, Y) internally
  * Chain-rule derivatives: d_dP_dP_static, d_dP_dT_up
  * FD derivatives for composition Y (no analytical drho/dY)
  * mach_clamped propagated to AreaChangeElementResult

- Wire through pybind:
  * sharp_area_change + conical_area_change
  * area_change_residuals_and_jacobian + conical_area_change_residuals_and_jacobian
  * mach_clamped exposed on both result structs

- Add comprehensive Jacobian tests (34 C++ tests):
  * Sharp low-level: dS_dm, dS_drho (expansion, contraction, reverse, zero flow)
  * Sharp solver-interface: d_dP_d_mdot, d_dP_dP_static, d_dP_dT (all geometries)
  * Conical low-level: dS_dm, dS_drho (diffuser, nozzle, reverse), dS_dmu = 0
  * Conical solver-interface: d_dP_d_mdot, d_dP_dP_static, d_dP_dT (expansion, contraction)
  * Conical consistency: length→0 recovers sharp-edge (within 15%)
  * Max relative error vs FD: 6.5e-5

- Add Python validation tests (24 tests):
  * Sharp: physics sanity, zero-crossing smoothness, Mach correction, mach_clamp flag
  * Conical: diffuser > nozzle loss, length→0 sharp recovery, dS_dmu = 0
  * Both: FD Jacobian accuracy (parametrized), network wrapper sanity

- Register in units_data.h for documentation generation
