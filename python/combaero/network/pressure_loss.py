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
consumed by ``P_out = P_in * (1 - loss_value)``.

WHERE THE LOSS IS EVALUATED. A ``PressureLossElement`` on the combustor's
inlet or outlet edge calls the correlation in Python. Its Jacobian takes the
returned ``d_loss/d_theta`` analytically. Every other dependence (m_dot, P_in,
T_in, composition) is a central difference over correlation calls. The C++
``combustor_residuals_and_jacobians`` also accepts a loss callable, but
``CombustorNode`` passes it a zero loss, so that chain does not run in a
network.

PYTHON BY DESIGN, FOR NOW (#718). The interface is a prototype for
user-supplied loss maps, e.g. a scipy interpolant built from rig data.
Combustors are described by their temperature rise (theta) and their fuel
splits, and a user-friendly input for the splits is not yet designed. Until
it is, the callable stays an arbitrary Python object, and the differences
are what let it be one.

Theta = T_ad / T_in - 1  (dimensionless temperature-rise ratio).
"""

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


def _low_flow_factor(ctx: object, m_ref: float | None) -> float:
    """1 at and above the design flow, falling to 0 like 3 (m/m_ref)^2 below.

    A fraction loss is a head loss quoted at its design point: below it the
    drop scales with the dynamic head, and at zero flow it vanishes. Without
    this a constant fraction holds a pressure difference at zero flow, and a
    reversing flow must jump through it -- the residual steps by xi Pt where
    the flow changes sign, which the solver cannot cross (#550).
    ``1 - (1 - r^2)^3`` is C2 at r = 1.
    """
    if m_ref is None:
        return 1.0
    if m_ref <= 0.0:
        raise ValueError(f"m_ref must be positive, got {m_ref}")
    r = abs(float(ctx.mdot_total)) / m_ref  # type: ignore[attr-defined]
    if r >= 1.0:
        return 1.0
    return 1.0 - (1.0 - r * r) ** 3


@dataclass
class ConstantFractionLoss:
    """Constant fractional total-pressure loss.

    Formula:  dP/P_in = xi   =>   P_out = P_in * (1 - xi)

    Parameters
    ----------
    xi : float
        Fractional total-pressure loss [-].  E.g. 0.03 = 3% of P_in.
    m_ref : float | None
        Design mass flow [kg/s]. Below it the loss scales with the dynamic
        head and vanishes at zero flow (``_low_flow_factor``); at and above
        it the loss is exactly ``xi``. Set it wherever the flow may run low
        or reverse. ``None`` keeps a constant fraction at every flow.

    Notes
    -----
    ``d(dP/P)/dTheta = 0`` -- loss is independent of temperature rise,
    so the pressure Jacobian contribution is zero.
    """

    xi: float = 0.03
    m_ref: float | None = None

    def __call__(self, ctx: object) -> tuple[float, float]:
        return self.xi * _low_flow_factor(ctx, self.m_ref), 0.0


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
    m_ref : float | None
        Design mass flow [kg/s], as for :class:`ConstantFractionLoss`.

    Notes
    -----
    Returns the analytical derivative ``d(dP/P)/dTheta = k``, which the
    ``PressureLossElement`` Jacobian uses directly for its theta column.
    """

    k: float = 0.5
    xi0: float = 0.02
    m_ref: float | None = None

    def __call__(self, ctx: object) -> tuple[float, float]:
        theta = ctx.theta  # type: ignore[attr-defined]
        g = _low_flow_factor(ctx, self.m_ref)
        return (self.k * theta + self.xi0) * g, self.k * g


# ---------------------------------------------------------------------------
# Head family:  dP = zeta * 0.5 * rho_in * v_in^2  (Euler loss coefficient)
#
# v_in = mdot_total / (rho_in * area)
# rho_in is computed from ctx.state_in (T, P, X) via the ideal-gas mixture.
#
# The return value is still a *fraction* of P_in so the C++ solver can use
# P_out = P_in*(1 - loss_value) uniformly.  The conversion is:
#   xi = dP/P_in = zeta * 0.5*rho_in*v_in^2 / P_in
#
# Jacobian note: d(xi)/d(theta) = 0 is used here (the dynamic-head term
# changes with theta only through rho_in(T_in), which is loop-free and
# treated as fixed for the Jacobian).  This is a valid approximation for
# the chain-rule through the Newton solver.
# ---------------------------------------------------------------------------


def _dynamic_head_fraction(ctx: object, area: float) -> float:
    """Compute 0.5*rho_in*v_in^2 / P_in from context and combustor area."""
    import combaero as cb

    T_in: float = ctx.T_in  # type: ignore[attr-defined]
    P_in: float = ctx.P_in  # type: ignore[attr-defined]
    X_in: list[float] = list(ctx.X_in)  # type: ignore[attr-defined]
    mdot_total: float = ctx.mdot_total  # type: ignore[attr-defined]
    if mdot_total <= 0.0 or area <= 0.0 or P_in <= 0.0:
        return 0.0
    state = cb.State()
    state.set_TPX(T_in, P_in, X_in)
    rho_in = state.rho
    if rho_in <= 0.0:
        return 0.0
    v_in = mdot_total / (rho_in * area)
    q_in = 0.5 * rho_in * v_in * v_in
    return q_in / P_in


@dataclass
class ConstantHeadLoss:
    """Constant Euler head-loss coefficient.

    Formula:  dP = zeta * 0.5 * rho_in * v_in^2
              v_in = mdot_total / (rho_in * area)

    Parameters
    ----------
    zeta : float
        Euler loss coefficient [-].
    area : float
        Combustor inlet cross-sectional area [m^2].

    Notes
    -----
    Internally converts to a pressure fraction xi = dP/P_in so the C++
    solver receives a consistent fractional loss.
    ``d(xi)/dTheta = 0`` is assumed (velocity reference is loop-free).
    """

    zeta: float = 5.0
    area: float = 0.1

    def __call__(self, ctx: object) -> tuple[float, float]:
        q_frac = _dynamic_head_fraction(ctx, self.area)
        return self.zeta * q_frac, 0.0


@dataclass
class LinearThetaHeadLoss:
    """Euler head-loss coefficient linear in temperature-rise ratio Theta.

    Formula:  dP = (k * Theta + zeta0) * 0.5 * rho_in * v_in^2
              Theta = T_ad / T_in - 1
              v_in  = mdot_total / (rho_in * area)

    Parameters
    ----------
    k : float
        Sensitivity of loss coefficient to temperature rise [-/-].
    zeta0 : float
        Cold-flow Euler loss coefficient [-] (value at Theta = 0).
    area : float
        Combustor inlet cross-sectional area [m^2].

    Notes
    -----
    ``d(xi)/dTheta = k * q_frac`` (analytical), where q_frac =
    0.5*rho_in*v_in^2 / P_in is treated as loop-free.
    """

    k: float = 1.0
    zeta0: float = 3.0
    area: float = 0.1

    def __call__(self, ctx: object) -> tuple[float, float]:
        theta = ctx.theta  # type: ignore[attr-defined]
        q_frac = _dynamic_head_fraction(ctx, self.area)
        xi = (self.k * theta + self.zeta0) * q_frac
        dxi_dtheta = self.k * q_frac
        return xi, dxi_dtheta


# ---------------------------------------------------------------------------
# Backward-compatible aliases (deprecated names kept for existing scripts)
# ---------------------------------------------------------------------------

ConstantLoss = ConstantFractionLoss
LinearThetaLoss = LinearThetaFractionLoss
