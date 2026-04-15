"""
Pressure-loss correlation classes for CombustorNode.

Each class is a callable that accepts a PressureLossContext and returns
``(zeta, dzeta_dtheta)`` -- the fractional total-pressure loss and its
analytical derivative w.r.t. the dimensionless temperature-rise ratio
theta = T_ad/T_in - 1.

The derivative is used by the C++ solver for analytical chain-rule
Jacobians through ``combustor_residuals_and_jacobians``.

Future extensions (not in this module):
  - FuelSplitLoss: zeta = f(theta, fuel_split) via safe AST expression
  - RBFInterpolatorLoss: serialized scipy RBFInterpolator for mapped data
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable


@runtime_checkable
class PressureLossCallable(Protocol):
    """Protocol for pressure-loss correlation objects.

    The context object passed at solve time is ``PressureLossContext``
    (exposed from C++ via pybind).  The duck-typed access on ``ctx.theta``
    is sufficient for the correlations in this module.
    """

    def __call__(self, ctx: object) -> tuple[float, float]: ...


@dataclass
class ConstantLoss:
    """Constant fractional total-pressure loss.

    Parameters
    ----------
    zeta0 : float
        Fractional pressure loss [-].  E.g. 0.03 means 3 % of P_in.

    Notes
    -----
    Returns ``dzeta/dtheta = 0`` so the Jacobian contribution of pressure
    loss is zero (loss does not depend on temperature rise).
    """

    zeta0: float = 0.03

    def __call__(self, ctx: object) -> tuple[float, float]:
        return self.zeta0, 0.0


@dataclass
class LinearThetaLoss:
    """Linear pressure loss in the temperature-rise ratio Theta = T_ad/T_in - 1.

    Correlation:  zeta = k * theta + zeta0

    Parameters
    ----------
    k : float
        Sensitivity to temperature rise ratio [-/-].
    zeta0 : float
        Base loss at Theta = 0 (e.g. cold-flow loss) [-].

    Notes
    -----
    Returns the analytical derivative ``dzeta/dtheta = k``, enabling the
    C++ solver to propagate the exact Jacobian through the pressure-loss
    term without finite-differencing.
    """

    k: float = 0.5
    zeta0: float = 0.02

    def __call__(self, ctx: object) -> tuple[float, float]:
        theta = ctx.theta  # type: ignore[attr-defined]
        return self.k * theta + self.zeta0, self.k
