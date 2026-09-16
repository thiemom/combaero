"""Internal cooling channel submodule.

Provides a high-level API for combined convective heat transfer and pressure
loss in internal cooling channels.  All functions accept ``(T, P, X, ...)``
and return a :class:`ChannelResult`.

Available models
----------------
smooth
    Smooth channel or duct (Gnielinski / Dittus-Boelter / Sieder-Tate / Petukhov).
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
- ``q = h * (T_aw - T_hot)`` when ``T_hot`` is supplied, else ``nan``.
"""

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

from ._core import ChannelResult
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
    T_hot: float = math.nan,
    correlation: str = "gnielinski",
    heating: bool = True,
    mu_ratio: float = 1.0,
    roughness: float = 0.0,
    Nu_multiplier: float = 1.0,
    f_multiplier: float = 1.0,
) -> ChannelResult:
    """Smooth channel or duct: combined HTC + pressure drop.

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
    T_hot:
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
        T_hot=T_hot,
        correlation=correlation,
        heating=heating,
        mu_ratio=mu_ratio,
        roughness=roughness,
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
    SmoothModel,
)

# Mark as re-exported to prevent ruff F401 errors
__all__ = [
    "ChannelModel",
    "ConvectiveSurface",
    "SmoothModel",
]
