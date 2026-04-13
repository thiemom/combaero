"""Internal cooling channel submodule.

Provides a high-level API for combined convective heat transfer and pressure
loss in internal cooling channels.  All functions accept ``(T, P, X, ...)``
and return a :class:`ChannelResult`.

Available models
----------------
smooth
    Smooth pipe or duct (Gnielinski / Dittus-Boelter / Sieder-Tate / Petukhov).
ribbed
    Rib-enhanced cooling channel (Han et al. 1988).
dimpled
    Dimpled surface cooling channel (Chyu et al. 1997).
pin_fin
    Pin-fin array cooling channel (Metzger et al. 1982).
impingement
    Impingement jet array cooling (Florschuetz et al. 1981 / Martin 1977).

Design notes
------------
- ``ChannelResult`` is returned directly, not ``FlowSolution``.  Use
  ``combaero.incompressible`` or ``combaero.compressible`` for flow-only
  (no heat transfer) elements.
- ``T_aw`` is always computed continuously for all Mach numbers:
  ``T_aw = T_static + r * v^2 / (2*cp)`` with ``r = Pr^(1/3)`` (turbulent).
  At ``v=0`` this reduces to ``T_static`` exactly - no threshold, no kink in
  the Jacobian.
- ``q = h * (T_aw - T_wall)`` when ``T_wall`` is supplied, else ``nan``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

from ._core import ChannelResult
from ._core import channel_dimpled as _channel_dimpled
from ._core import channel_impingement as _channel_impingement
from ._core import channel_pin_fin as _channel_pin_fin
from ._core import channel_ribbed as _channel_ribbed
from ._core import channel_smooth as _channel_smooth

if TYPE_CHECKING:
    pass


def smooth(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    L: float,
    D: float,
    T_wall: float = math.nan,
    correlation: str = "gnielinski",
    heating: bool = True,
    mu_ratio: float = 1.0,
    roughness: float = 0.0,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Smooth pipe or duct: combined HTC + pressure drop.

    Parameters
    ----------
    T:
        Bulk static temperature [K].
    P:
        Bulk static pressure [Pa].
    X:
        Mole fractions [-].
    u:
        Bulk flow velocity [m/s].
    L:
        Channel length [m].
    D:
        Hydraulic diameter [m].
    T_wall:
        Wall temperature [K].  Supply to obtain ``q``; omit (``nan``) for
        flow-only or when wall temperature is unknown.
    correlation:
        Nusselt correlation: ``"gnielinski"`` (default), ``"dittus_boelter"``,
        ``"sieder_tate"``, ``"petukhov"``.
    heating:
        ``True`` if the fluid is being heated (affects Dittus-Boelter exponent).
    mu_ratio:
        ``mu_bulk / mu_wall`` for Sieder-Tate viscosity correction.
    roughness:
        Absolute wall roughness [m].  ``0.0`` for hydraulically smooth.

    Returns
    -------
    ChannelResult
    """
    return _channel_smooth(
        T,
        P,
        list(X),
        u,
        D,
        L,
        T_wall=T_wall,
        correlation=correlation,
        heating=heating,
        mu_ratio=mu_ratio,
        roughness=roughness,
        Nu_multiplier=Nu_multiplier,
        f_multiplier=f_multiplier,
    )


def ribbed(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    L: float,
    D: float,
    e_D: float,
    pitch_to_height: float,
    alpha_deg: float,
    T_wall: float = math.nan,
    heating: bool = True,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Rib-enhanced cooling channel (Han et al. 1988).

    Applies ``rib_enhancement_factor`` to Nu and ``rib_friction_multiplier``
    to f from a smooth-pipe Gnielinski baseline.

    Parameters
    ----------
    T, P, X:
        Bulk static thermodynamic state.
    u:
        Bulk flow velocity [m/s].
    L:
        Channel length [m].
    D:
        Hydraulic diameter [m].
    e_D:
        Rib height / hydraulic diameter [-].  Valid: 0.02-0.1.
    pitch_to_height:
        Rib pitch / rib height [-].  Valid: 5-20.
    alpha_deg:
        Rib angle [deg].  Valid: 30-90.
    T_wall:
        Wall temperature [K].  ``nan`` to skip ``q``.
    heating:
        ``True`` if the fluid is being heated.

    Returns
    -------
    ChannelResult
    """
    return _channel_ribbed(
        T,
        P,
        list(X),
        u,
        D,
        L,
        e_D,
        pitch_to_height,
        alpha_deg,
        T_wall=T_wall,
        heating=heating,
        Nu_multiplier=Nu_multiplier,
        f_multiplier=f_multiplier,
    )


def dimpled(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    L: float,
    D: float,
    d_Dh: float,
    h_d: float,
    S_d: float,
    T_wall: float = math.nan,
    heating: bool = True,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Dimpled surface cooling channel (Chyu et al. 1997).

    Applies ``dimple_nusselt_enhancement`` to Nu and
    ``dimple_friction_multiplier`` to f from a smooth-pipe Gnielinski baseline.

    Parameters
    ----------
    T, P, X:
        Bulk static thermodynamic state.
    u:
        Bulk flow velocity [m/s].
    L:
        Channel length [m].
    D:
        Hydraulic diameter [m].
    d_Dh:
        Dimple diameter / channel height [-].  Valid: 0.1-0.3.
    h_d:
        Dimple depth / diameter [-].  Valid: 0.1-0.3.
    S_d:
        Dimple spacing / diameter [-].  Valid: 1.5-3.0.
    T_wall:
        Wall temperature [K].  ``nan`` to skip ``q``.
    heating:
        ``True`` if the fluid is being heated.

    Returns
    -------
    ChannelResult
    """
    return _channel_dimpled(
        T,
        P,
        list(X),
        u,
        D,
        L,
        d_Dh,
        h_d,
        S_d,
        T_wall=T_wall,
        heating=heating,
        Nu_multiplier=Nu_multiplier,
        f_multiplier=f_multiplier,
    )


def pin_fin(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    u: float,
    H: float,
    d: float,
    S_D: float,
    X_D: float,
    N_rows: int,
    T_wall: float = math.nan,
    is_staggered: bool = True,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Pin-fin array cooling channel (Metzger et al. 1982).

    Re and Nu are based on pin diameter ``d``.
    ``dP = N_rows * f_pin * (rho * v_max^2 / 2)``
    where ``v_max = u * S_D / (S_D - 1)`` is the minimum cross-section velocity.

    Parameters
    ----------
    T, P, X:
        Bulk static thermodynamic state.
    u:
        Approach (upstream) velocity [m/s].
    H:
        Pin length / channel height [m].
    d:
        Pin diameter [m].
    S_D:
        Spanwise pitch / d [-].  Valid: 1.5-4.0.
    X_D:
        Streamwise pitch / d [-].  Valid: 1.5-4.0.
    N_rows:
        Number of pin rows in the streamwise direction.
    T_wall:
        Wall temperature [K].  ``nan`` to skip ``q``.
    is_staggered:
        ``True`` for staggered array (default), ``False`` for inline.

    Returns
    -------
    ChannelResult
    """
    return _channel_pin_fin(
        T,
        P,
        list(X),
        u,
        H,
        d,
        S_D,
        X_D,
        N_rows,
        T_wall=T_wall,
        is_staggered=is_staggered,
        Nu_multiplier=Nu_multiplier,
        f_multiplier=f_multiplier,
    )


def impingement(
    T: float,
    P: float,
    X: Sequence[float],
    *,
    mdot_jet: float,
    d_jet: float,
    z_D: float,
    x_D: float = 0.0,
    y_D: float = 0.0,
    A_target: float,
    T_wall: float = math.nan,
    Cd_jet: float = 0.65,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Impingement jet array cooling (Florschuetz et al. 1981 / Martin 1977).

    Re and Nu are based on jet diameter ``d_jet``.
    ``dP = (1/Cd_jet^2) * rho * v_jet^2 / 2``  (jet-plate orifice loss).

    Parameters
    ----------
    T, P, X:
        Bulk static thermodynamic state of the coolant.
    mdot_jet:
        Total jet mass flow rate [kg/s].
    d_jet:
        Jet hole diameter [m].
    z_D:
        Jet-to-target distance / ``d_jet`` [-].  Valid: 1-12.
    x_D:
        Streamwise jet spacing / ``d_jet`` [-].  ``0.0`` for single jet.
        Valid for arrays: 4-16.
    y_D:
        Spanwise jet spacing / ``d_jet`` [-].  ``0.0`` for single jet.
        Valid for arrays: 4-16.
    A_target:
        Target surface area [m^2].
    T_wall:
        Wall temperature [K].  ``nan`` to skip ``q``.
    Cd_jet:
        Jet hole discharge coefficient [-].  Default ``0.65``.

    Returns
    -------
    ChannelResult
    """
    return _channel_impingement(
        T,
        P,
        list(X),
        mdot_jet,
        d_jet,
        z_D,
        x_D,
        y_D,
        A_target,
        T_wall=T_wall,
        Cd_jet=Cd_jet,
        Nu_multiplier=Nu_multiplier,
        f_multiplier=f_multiplier,
    )


# ============================================================================
# Backward-compatible re-exports from network.components
# ============================================================================

# Re-export convective heat transfer classes for backward compatibility
# ruff: noqa: E402  (import not at top of file - needed for circular import avoidance)
from combaero.network.components import (
    ChannelModel,
    ConvectiveSurface,
    DimpledModel,
    ImpingementModel,
    PinFinModel,
    RibbedModel,
    SmoothModel,
    WallConnection,
)

# Mark as re-exported to prevent ruff F401 errors
__all__ = [
    "ChannelModel",
    "ConvectiveSurface",
    "DimpledModel",
    "ImpingementModel",
    "PinFinModel",
    "RibbedModel",
    "SmoothModel",
    "WallConnection",
]
