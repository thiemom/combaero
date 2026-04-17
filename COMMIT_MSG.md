test: add combustor regression tests for pressure loss and mass balance

Add comprehensive regression tests for CombustorNode solver bugs:
- test_combustor_lossless_inlets_mass_balance: guards the fix for
  unconstrained comb.P causing mass imbalance
- test_combustor_channel_outlet_network_converges: guards the original
  failing topology (2x MassFlowBoundary -> Lossless -> Combustor -> Channel)
- test_combustor_pressure_loss_computes_mix_result: verifies mixer
  computes P_total_mix correctly
- test_combustor_pressure_loss_affects_solution (xfail strict): documents
  that ALL pressure loss models have no effect on the solution
- test_combustor_lossless_inlets_pressure_loss_enforced (xfail strict):
  documents loss ignored with LosslessConnectionElement inlets

The parametric test covers all loss variants:
- ConstantFractionLoss
- LinearThetaFractionLoss
- ConstantHeadLoss
- LinearThetaHeadLoss

All are marked xfail(strict=True) so they become XPASS when fixed.

Also adds a code comment in CombustorNode.residuals noting the known bug.
