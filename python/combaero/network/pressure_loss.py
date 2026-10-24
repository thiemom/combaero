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

WHERE THE LOSS IS EVALUATED. The four built-in classes evaluate in C++
(``cb.combustor_fraction_loss`` / ``cb.combustor_head_loss``), and their
``evaluate(ctx, Y_in)`` returns every partial (theta, m_dot, T_in, P_in, each
inlet mass fraction). A ``PressureLossElement`` builds its Jacobian from those.
The C++ ``combustor_residuals_and_jacobians`` also accepts a loss callable,
but ``CombustorNode`` passes it a zero loss, so that chain does not run in a
network.

A USER CALLABLE STAYS PYTHON (#718). Any object returning ``(xi, dxi/dtheta)``
from a context still works. Its other dependences are central differences,
which are what let it be arbitrary Python. That interface is the prototype for
user-supplied loss maps (e.g. a scipy interpolant over theta and fuel splits;
the input design for the splits is open, #725).

Theta = T_ad / T_in - 1  (dimensionless temperature-rise ratio).
"""

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from .._core import (
    CombustorLossResult,
    combustor_fraction_loss,
    combustor_head_loss,
    mole_to_mass,
)


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


def _m_ref_arg(m_ref: float | None) -> float:
    """The design flow for C++ (0 = no fade), validated.

    A fraction loss is a head loss quoted at its design point: below it the
    drop scales with the dynamic head (``1 - (1 - r^2)^3``, C2 at r = 1), and
    at zero flow it vanishes. Without this a constant fraction holds a
    pressure difference at zero flow, and a reversing flow must jump through
    it -- the residual steps by xi Pt where the flow changes sign, which the
    solver cannot cross (#550).
    """
    if m_ref is None:
        return 0.0
    if m_ref <= 0.0:
        raise ValueError(f"m_ref must be positive, got {m_ref}")
    return float(m_ref)


def _inlet_Y(ctx: object, Y_in: list[float] | None) -> list[float]:
    if Y_in is not None:
        return list(Y_in)
    return list(mole_to_mass(list(ctx.X_in)))  # type: ignore[attr-defined]


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
        head and vanishes at zero flow (``_m_ref_arg``); at and above
        it the loss is exactly ``xi``. Set it wherever the flow may run low
        or reverse. ``None`` keeps a constant fraction at every flow.

    Notes
    -----
    ``d(dP/P)/dTheta = 0`` -- loss is independent of temperature rise,
    so the pressure Jacobian contribution is zero.
    """

    xi: float = 0.03
    m_ref: float | None = None

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> CombustorLossResult:
        """xi and its partials, from C++."""
        # m_dot matters only below a design flow; a context may omit it.
        m = float(ctx.mdot_total) if self.m_ref is not None else 0.0  # type: ignore[attr-defined]
        return combustor_fraction_loss(0.0, m, 0.0, self.xi, _m_ref_arg(self.m_ref))

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


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

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> CombustorLossResult:
        """xi and its partials, from C++."""
        theta = float(ctx.theta)  # type: ignore[attr-defined]
        # m_dot matters only below a design flow; a context may omit it.
        m = float(ctx.mdot_total) if self.m_ref is not None else 0.0  # type: ignore[attr-defined]
        return combustor_fraction_loss(theta, m, self.k, self.xi0, _m_ref_arg(self.m_ref))

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


# ---------------------------------------------------------------------------
# Head family:  dP = zeta * 0.5 * rho_in * v_in^2  (Euler loss coefficient)
#
# v_in = mdot_total / (rho_in * area), rho_in the ideal-gas density at
# (T_in, P_in, inlet composition). The return value is still a *fraction* of
# P_in, xi = zeta * 0.5 rho_in v_in^2 / P_in, so P_out = P_in (1 - xi) holds
# for both families. The partials in m_dot, T_in, P_in and composition come
# from C++ with the value (``cb.combustor_head_loss``).
# ---------------------------------------------------------------------------


def _head(
    ctx: object, Y_in: list[float] | None, k: float, zeta0: float, area: float
) -> CombustorLossResult:
    return combustor_head_loss(
        float(ctx.theta),  # type: ignore[attr-defined]
        float(ctx.mdot_total),  # type: ignore[attr-defined]
        float(ctx.T_in),  # type: ignore[attr-defined]
        float(ctx.P_in),  # type: ignore[attr-defined]
        _inlet_Y(ctx, Y_in),
        k,
        zeta0,
        area,
    )


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
    Internally converts to a pressure fraction xi = dP/P_in so both families
    share P_out = P_in (1 - xi).
    """

    zeta: float = 5.0
    area: float = 0.1

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> CombustorLossResult:
        """xi and its partials, from C++."""
        return _head(ctx, Y_in, 0.0, self.zeta, self.area)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


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
    ``d(xi)/dTheta = k * q_frac``, q_frac = 0.5*rho_in*v_in^2 / P_in; the
    other partials come with it from ``evaluate``.
    """

    k: float = 1.0
    zeta0: float = 3.0
    area: float = 0.1

    def evaluate(self, ctx: object, Y_in: list[float] | None = None) -> CombustorLossResult:
        """xi and its partials, from C++."""
        return _head(ctx, Y_in, self.k, self.zeta0, self.area)

    def __call__(self, ctx: object) -> tuple[float, float]:
        r = self.evaluate(ctx)
        return r.xi, r.dxi_dtheta


# ---------------------------------------------------------------------------
# Backward-compatible aliases (deprecated names kept for existing scripts)
# ---------------------------------------------------------------------------

ConstantLoss = ConstantFractionLoss
LinearThetaLoss = LinearThetaFractionLoss
