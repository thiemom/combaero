"""
Pressure-loss correlation classes for CombustorNode.

Two families of correlations are provided:

**Fraction family** -- loss expressed as a fraction of inlet total pressure:
    dP_total = xi * P_in        =>   P_out = P_in * (1 - xi)
  Classes: ConstantFractionLoss, LinearThetaFractionLoss

**Head family** -- loss expressed as a multiple of the inlet dynamic head:
    dP = zeta * 0.5 * rho_in * v_in^2
  Classes: ConstantHeadLoss, LinearThetaHeadLoss  (Phase 1b, needs mdot_total in ctx)

All callables return ``(loss_value, d_loss_value_d_theta)`` where
``loss_value`` is the *fractional* total-pressure loss (dimensionless)
consumed by ``P_out = P_in * (1 - loss_value)`` in the C++ solver.
The derivative enables exact analytical Jacobians through
``combustor_residuals_and_jacobians``.

Theta = T_ad / T_in - 1  (dimensionless temperature-rise ratio).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable


@runtime_checkable
class PressureLossCallable(Protocol):
    """Protocol for pressure-loss correlation objects.

    The context object passed at solve time is ``PressureLossContext``
    (exposed from C++ via pybind).
    """

    def __call__(self, ctx: object) -> tuple[float, float]: ...


# ---------------------------------------------------------------------------
# Fraction family:  dP/P_in = xi   (total-pressure ratio convention)
# ---------------------------------------------------------------------------


@dataclass
class ConstantFractionLoss:
    """Constant fractional total-pressure loss.

    Formula:  dP/P_in = xi   =>   P_out = P_in * (1 - xi)

    Parameters
    ----------
    xi : float
        Fractional total-pressure loss [-].  E.g. 0.03 = 3% of P_in.

    Notes
    -----
    ``d(dP/P)/dTheta = 0`` -- loss is independent of temperature rise,
    so the pressure Jacobian contribution is zero.
    """

    xi: float = 0.03

    def __call__(self, ctx: object) -> tuple[float, float]:
        return self.xi, 0.0


@dataclass
class LinearThetaFractionLoss:
    """Total-pressure loss fraction linear in temperature-rise ratio Theta.

    Formula:  dP/P_in = k * Theta + xi0
              Theta   = T_ad / T_in - 1

    Parameters
    ----------
    k : float
        Sensitivity of fractional loss to temperature rise [-/-].
        Positive k means more heat release = higher pressure loss.
    xi0 : float
        Cold-flow base loss fraction [-] (value at Theta = 0).

    Notes
    -----
    Returns the analytical derivative ``d(dP/P)/dTheta = k``, enabling
    exact chain-rule Jacobians in the C++ solver without finite-differencing.
    """

    k: float = 0.5
    xi0: float = 0.02

    def __call__(self, ctx: object) -> tuple[float, float]:
        theta = ctx.theta  # type: ignore[attr-defined]
        return self.k * theta + self.xi0, self.k


# ---------------------------------------------------------------------------
# Backward-compatible aliases (deprecated names kept for existing scripts)
# ---------------------------------------------------------------------------

ConstantLoss = ConstantFractionLoss
LinearThetaLoss = LinearThetaFractionLoss
