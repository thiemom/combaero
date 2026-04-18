fix(solver): replace asymmetric FD with exact analytical dT0/dM and dP0/dM

The previous implementations of T0_from_static_and_jacobian_M and
P0_from_static_and_jacobian_M used an asymmetric one-sided finite
difference (smooth-floored M_minus vs. forward M_plus), which produced
O(h) error (~0.05 absolute difference). This caused platform-specific CI
failures because the test used a central FD reference.

Root cause: the smooth-max floor shifts M_minus by a different amount
than M_plus, making the "central" difference evaluate at an asymmetric
interval and producing O(h) truncation error instead of O(h^2).

Fix: replace internal FD with exact analytical chain-rule derivatives:

  dT0/dM = M * a(T)^2 / cp_mass(T0)
    (from h(T0) = h(T) + 0.5*(M*a)^2, differentiate w.r.t. M)

  dP0/dM = P0 * cp_mass(T0) / (T0 * R_specific) * dT0/dM
    (from P0 = P_REF * exp([s(T0,P_REF) - s(T,P)] / R_specific),
     differentiate w.r.t. T0, apply chain rule)

Tests now pass with relative error < 1e-9 against the external central FD
(maximum absolute difference: 6.8e-5 Pa), well within the 1e-6 tolerance.

Also fixes PressureLossElement Jacobian sign error for linear-theta
correlations (dres/dT_burned was negated), and updates PR_BODY.md to
correctly scope the GUI as the primary PR deliverable.
