"""Unified FlowSolution dataclass shared by incompressible and compressible submodules.

Fields that are not applicable to a given regime are set to ``math.nan``
(numeric) or ``False`` (boolean) so that downstream code can always access
every attribute without branching on the flow regime.

Regime mapping
--------------
Incompressible pipe / orifice
    mdot, v, dP, Re, rho, f, Cd, regime="incompressible"
    choked=False, M=nan, T_out=nan, P_out=nan, h0=nan, L_choke=nan

Compressible Fanno pipe
    mdot, v, dP, Re, rho, f, M, T_out, P_out, h0, choked, L_choke
    regime="compressible", Cd=nan

Compressible nozzle
    mdot, v, M, T_out, P_out, h0, choked
    regime="compressible", dP=nan, Re=nan, rho=nan, f=nan, Cd=nan, L_choke=nan
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class FlowSolution:
    """Unified result of a high-level flow calculation.

    All fields are always present.  Fields that are not meaningful for the
    current regime carry ``math.nan`` (floats) or ``False`` (bools).

    Parameters
    ----------
    regime:
        ``"incompressible"`` or ``"compressible"``.
    mdot:
        Mass flow rate [kg/s].
    v:
        Characteristic velocity (pipe bulk velocity or nozzle exit velocity)
        [m/s].
    dP:
        Pressure drop P_in − P_out [Pa].  ``nan`` for nozzle flow.
    Re:
        Reynolds number [−].  ``nan`` when not computed.
    rho:
        Inlet density [kg/m³].  ``nan`` when not computed.
    f:
        Darcy friction factor [−].  ``nan`` for orifice / nozzle.
    Cd:
        Discharge coefficient [−].  ``nan`` for pipe / nozzle.
    M:
        Outlet Mach number [−].  ``nan`` for incompressible.
    T_out:
        Outlet static temperature [K].  ``nan`` for incompressible.
    P_out:
        Outlet static pressure [Pa].  ``nan`` for incompressible.
    h0:
        Stagnation enthalpy [J/kg].  ``nan`` for incompressible.
    choked:
        ``True`` if the flow is choked (Fanno or nozzle).
    L_choke:
        Length to choking [m] (Fanno only).  ``nan`` otherwise.
    """

    regime: str
    mdot: float
    v: float
    dP: float = math.nan  # noqa: N815
    Re: float = math.nan  # noqa: N815
    rho: float = math.nan
    f: float = math.nan
    Cd: float = math.nan
    M: float = math.nan
    T_out: float = math.nan
    P_out: float = math.nan
    h0: float = math.nan
    choked: bool = False
    L_choke: float = math.nan
    # Optional axial profile (list of FannoStation or similar); empty by default.
    profile: list = field(default_factory=list, compare=False)

    def __repr__(self) -> str:
        def _fmt(v: object) -> str:
            if isinstance(v, float):
                return "nan" if math.isnan(v) else f"{v:.6g}"
            return repr(v)

        parts = [f"regime={self.regime!r}", f"mdot={_fmt(self.mdot)}"]
        if not math.isnan(self.dP):
            parts.append(f"dP={_fmt(self.dP)}")
        if not math.isnan(self.Re):
            parts.append(f"Re={_fmt(self.Re)}")
        if not math.isnan(self.M):
            parts.append(f"M={_fmt(self.M)}")
        if self.choked:
            parts.append("choked=True")
        return f"FlowSolution({', '.join(parts)})"
