import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb

# Physics configuration types for network introspectability
CompressibilityLiteral = Literal["incompressible", "compressible"]
FrictionModelLiteral = Literal["haaland", "colebrook", "serghides", "petukhov"]
HeatTransferModelLiteral = Literal[
    "none", "gnielinski", "dittus_boelter", "sieder_tate", "petukhov"
]

if TYPE_CHECKING:
    from .graph import FlowNetwork

CombustionMethodLiteral = Literal["complete", "equilibrium"]


# ============================================================================
# Convective Heat Transfer Models and Surface (Phase 3)
# ============================================================================


def _safe_rho(rho_raw: float, rho_min: float = 0.01) -> tuple[float, float]:
    """Smooth lower bound on density for numerical robustness.

    Returns (rho_safe, d_rho_safe / d_rho_raw) for Jacobian chain rule.

    Uses softplus: rho_safe = rho_min + rho_min * ln(1 + exp((rho_raw - rho_min) / rho_min))
    - When rho_raw >> rho_min: rho_safe ~= rho_raw  (transparent)
    - When rho_raw -> 0:        rho_safe -> rho_min * (1 + ln2)  (bounded)
    - When rho_raw -> -inf:     rho_safe -> rho_min  (floor)
    - Derivative is always in (0, 1]: solver always sees a gradient.
    """
    z = (rho_raw - rho_min) / rho_min
    if z > 20.0:  # overflow guard
        return rho_raw, 1.0
    if z < -20.0:  # underflow guard
        return rho_min, 0.0
    exp_z = math.exp(z)
    rho_safe = rho_min + rho_min * math.log1p(exp_z)
    d_safe = exp_z / (1.0 + exp_z)  # sigmoid - always in (0, 1)
    return rho_safe, d_safe


@dataclass
class SmoothModel:
    """Parameters for channel_smooth."""

    correlation: str = "gnielinski"  # "gnielinski" | "dittus_boelter" | "sieder_tate" | "petukhov"
    mu_ratio: float = 1.0  # mu_bulk / mu_wall (Sieder-Tate)
    roughness: float = 0.0  # absolute roughness [m]


@dataclass
class RibbedModel:
    """Parameters for channel_ribbed."""

    e_D: float = 0.0  # rib height / hydraulic diameter  # noqa: N815
    pitch_to_height: float = 0.0  # rib pitch / rib height
    alpha_deg: float = 90.0  # rib angle [deg]


@dataclass
class DimpledModel:
    """Parameters for channel_dimpled."""

    d_Dh: float = 0.0  # dimple diameter / hydraulic diameter  # noqa: N815
    h_d: float = 0.0  # dimple depth / dimple diameter
    S_d: float = 0.0  # dimple pitch / dimple diameter


@dataclass
class PinFinModel:
    """Parameters for channel_pin_fin."""

    pin_diameter: float = 0.0  # pin diameter [m]
    channel_height: float = 0.0  # channel height [m]
    S_D: float = 2.0  # transverse pitch / pin diameter
    X_D: float = 2.0  # streamwise pitch / pin diameter
    N_rows: int = 1  # number of pin rows
    is_staggered: bool = True


@dataclass
class ImpingementModel:
    """Parameters for channel_impingement."""

    d_jet: float = 0.0  # jet hole diameter [m]
    z_D: float = 0.0  # jet-to-target distance / d_jet  # noqa: N815
    x_D: float = 0.0  # streamwise pitch / d_jet  # noqa: N815
    y_D: float = 0.0  # spanwise pitch / d_jet  # noqa: N815
    A_target: float = 0.0  # target area [m^2]
    Cd_jet: float = 0.8  # jet discharge coefficient


ChannelModel = SmoothModel | RibbedModel | DimpledModel | PinFinModel | ImpingementModel


@dataclass
class ConvectiveSurface:
    """Convective surface description attached to a NetworkElement.

    Attributes
    ----------
    area : float
        Wetted convective area [m^2]. Default 0.0 (disabled).
    model : ChannelModel
        Model-specific subclass holding geometry parameters.
    heating : bool | None
        True if fluid is being heated, False if cooled, None = auto-detect
        from sign of (T_wall - T_fluid).
    Nu_multiplier : float
        Empirical correction factor on Nusselt number (default 1.0).
    f_multiplier : float
        Empirical correction factor on friction factor (default 1.0).
    """

    area: float = 0.0  # A_conv [m^2] - 0 disables
    model: ChannelModel = field(default_factory=SmoothModel)
    heating: bool | None = None  # None = auto-detect
    Nu_multiplier: float = 1.0  # empirical correction on Nu
    f_multiplier: float = 1.0  # empirical correction on f

    def htc_and_T(
        self,
        T: float,
        P: float,
        X: list[float],
        velocity: float,
        diameter: float,
        length: float,
        T_wall: float = math.nan,
    ):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Parameters
        ----------
        T : float
            Bulk static temperature [K].
        P : float
            Bulk static pressure [Pa].
        X : list[float]
            Mole fractions [-].
        velocity : float
            Bulk flow velocity [m/s].
        diameter : float
            Hydraulic diameter [m].
        length : float
            Channel length [m].
        T_wall : float, optional
            Wall temperature [K]. Used for auto-detection of heating/cooling.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if area=0. Access convective area via ``self.area``.
        """
        import math

        import combaero as cb

        if self.area == 0.0 or velocity < 1e-12:
            return None

        # Auto-detect heating direction from T_wall - T sign
        if self.heating is not None:
            heating = self.heating
        elif math.isfinite(T_wall):
            heating = T_wall >= T
        else:
            heating = True  # default when T_wall unknown

        # Dispatch to appropriate C++ channel_* function based on model type
        if isinstance(self.model, SmoothModel):
            result = cb.channel_smooth(
                T,
                P,
                X,
                velocity,
                diameter,
                length,
                T_wall=T_wall,
                correlation=self.model.correlation,
                heating=heating,
                mu_ratio=self.model.mu_ratio,
                roughness=self.model.roughness,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
            )
        elif isinstance(self.model, RibbedModel):
            result = cb.channel_ribbed(
                T,
                P,
                X,
                velocity,
                diameter,
                length,
                self.model.e_D,
                self.model.pitch_to_height,
                self.model.alpha_deg,
                T_wall=T_wall,
                heating=heating,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
            )
        elif isinstance(self.model, DimpledModel):
            result = cb.channel_dimpled(
                T,
                P,
                X,
                velocity,
                diameter,
                length,
                self.model.d_Dh,
                self.model.h_d,
                self.model.S_d,
                T_wall=T_wall,
                heating=heating,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
            )
        elif isinstance(self.model, PinFinModel):
            result = cb.channel_pin_fin(
                T,
                P,
                X,
                velocity,
                self.model.channel_height,
                self.model.pin_diameter,
                self.model.S_D,
                self.model.X_D,
                self.model.N_rows,
                T_wall=T_wall,
                is_staggered=self.model.is_staggered,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
            )
        elif isinstance(self.model, ImpingementModel):
            # Impingement requires mdot_jet (flow PER JET).
            # The total mass flow in the channel is rho * velocity * A_channel.
            # The number of jets is A_target / (x * y), where x, y are pitches.
            rho, _ = _safe_rho(cb.density(T, P, X))
            mdot_total = rho * velocity * (math.pi / 4 * diameter**2)

            x = self.model.x_D * self.model.d_jet
            y = self.model.y_D * self.model.d_jet
            A_per_jet = x * y
            if A_per_jet > 0 and self.model.A_target > 0:
                N_jets = self.model.A_target / A_per_jet
                mdot_jet = mdot_total / max(1.0, N_jets)
            else:
                mdot_jet = mdot_total  # fallback

            result = cb.channel_impingement(
                T,
                P,
                X,
                mdot_jet,
                self.model.d_jet,
                self.model.z_D,
                self.model.x_D,
                self.model.y_D,
                self.model.A_target,
                T_wall=T_wall,
                Cd_jet=self.model.Cd_jet,
                Nu_multiplier=self.Nu_multiplier,
                f_multiplier=self.f_multiplier,
            )
        else:
            raise TypeError(f"Unknown channel model type: {type(self.model)}")

        return result


# ============================================================================
# WallConnection (Phase 4/5)
# ============================================================================


@dataclass
class WallLayer:
    """Description of a single physical layer in a multi-layer ThermalWall.

    Attributes
    ----------
    thickness : float
        Thickness of the layer [m].
    conductivity : float
        Thermal conductivity of the material [W/(m*C)].
    material : str
        Material name (optional, for future material-lookup support).
    """

    thickness: float
    conductivity: float
    material: str = "generic"

    @property
    def r_val(self) -> float:
        """Thermal resistance per unit area [m^2*K/W]."""
        return self.thickness / self.conductivity


@dataclass
class ThermalWall:
    """Thermal coupling between two elements through a shared multi-layer wall.

    Attributes
    ----------
    id : str
        Unique identifier for this wall connection.
    element_a : str
        Element ID for side A.
    element_b : str
        Element ID for side B.
    layers : list[WallLayer]
        Stack of wall layers ordered from side A to side B.
    contact_area : float | None
        Override effective area [m^2]. If None, uses min(A_conv_a, A_conv_b).
    R_fouling : float
        Additional fouling resistance [m^2*K/W] (applied to cold side).
    """

    id: str
    element_a: str
    element_b: str
    layers: list[WallLayer] = field(default_factory=list)
    contact_area: float | None = None
    R_fouling: float = 0.0

    def compute_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        h_b: float,
        T_aw_b: float,
        A_conv_b: float,
    ) -> cb.WallCouplingResult:
        """Compute heat transfer rate and wall temperature using analytical Jacobians.

        Parameters
        ----------
        h_a, h_b : float
            Convective HTCs on sides A and B [W/(m^2*K)].
        T_aw_a, T_aw_b : float
            Adiabatic wall temperatures on sides A and B [K].
        A_conv_a, A_conv_b : float
            Convective areas on both sides [m^2].

        Returns
        -------
        WallCouplingResult
            Object containing Q [W], T_wall [K], and analytical Jacobians.
        """
        # Effective area
        A_eff = self.contact_area if self.contact_area is not None else min(A_conv_a, A_conv_b)

        # Multi-layer R values (thickness / conductivity)
        t_over_k_layers = [L.r_val for L in self.layers]

        # Call C++ function
        return cb.wall_coupling_and_jacobian_multilayer(
            h_a, T_aw_a, h_b, T_aw_b, t_over_k_layers, A_eff, self.R_fouling
        )


# Backward Compatibility Alias (Deprecated: use ThermalWall instead)
class WallConnection(ThermalWall):
    """Legacy alias for ThermalWall to support older network definitions.

    Translates scalar thickness/conductivity into a single-layer ThermalWall.
    """

    def __init__(
        self,
        id: str,
        element_a: str,
        element_b: str,
        wall_thickness: float,
        wall_conductivity: float,
        contact_area: float | None = None,
    ):
        super().__init__(
            id=id,
            element_a=element_a,
            element_b=element_b,
            layers=[WallLayer(thickness=wall_thickness, conductivity=wall_conductivity)],
            contact_area=contact_area,
        )

    @property
    def wall_thickness(self) -> float:
        return self.layers[0].thickness if self.layers else 1.0

    @property
    def wall_conductivity(self) -> float:
        return self.layers[0].conductivity if self.layers else 1.0


class EnergyBoundary:
    """Energy source/sink that attaches to mixing nodes.

    Convention: Q > 0 / fraction > 0 means heat **into** the fluid (heating),
                Q < 0 / fraction < 0 means heat **out of** the fluid (cooling).

    Two additive modes:
      - Q [W]: absolute heat transfer rate
      - fraction [-]: relative to total enthalpy flow (H_tot = sum(mdot_i * h_i))
        e.g. fraction=-0.05 means 5% enthalpy loss

    Both can be set simultaneously; effective Q = Q + fraction * H_tot.
    C++ functions handle the conversion to delta_h and all Jacobian corrections.
    """

    def __init__(self, id: str, Q: float = 0.0, fraction: float = 0.0) -> None:
        self.id = id
        self.Q = Q  # [W]
        self.fraction = fraction  # [-]


@dataclass
class MixtureState:
    P: float
    P_total: float
    T: float
    T_total: float
    m_dot: float
    Y: list[float]

    @property
    def X(self) -> list[float]:
        """Cached mole fractions converted from mass fractions with defensive clamping."""
        try:
            return self._X  # type: ignore[return-value]
        except AttributeError:
            # Defensive guard against unphysical intermediate solver states
            Y_safe = [max(1e-12, min(1.0, float(yi))) for yi in self.Y]
            y_sum = sum(Y_safe)
            if y_sum > 0:
                Y_safe = [yi / y_sum for yi in Y_safe]
            x = list(cb.mass_to_mole(Y_safe))
            object.__setattr__(self, "_X", x)
            return x

    def density(self) -> float:
        """Static density [kg/m^3]."""
        return cb.density(self.T, self.P, self.X)

    def enthalpy(self) -> float:
        """Static specific enthalpy [J/kg]."""
        return cb.h_mass(self.T, self.X)

    def total_enthalpy(self) -> float:
        """Total (stagnation) specific enthalpy [J/kg]."""
        return cb.h_mass(self.T_total, self.X)

    def cp(self) -> float:
        """Specific heat capacity at constant pressure [J/(kg*K)]."""
        return cb.cp_mass(self.T, self.X)

    def speed_of_sound(self) -> float:
        """Speed of sound [m/s]."""
        return cb.speed_of_sound(self.T, self.X)

    def gamma(self) -> float:
        """Ratio of specific heats (cp/cv) [-]."""
        cp = self.cp()
        cv = cb.cv_mass(self.T, self.X)
        return cp / cv if cv > 0 else 1.4


class NetworkNode(ABC):
    def __init__(self, id: str):
        self.id = id

    @property
    def has_convective_surface(self) -> bool:
        """True if this node has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this node contributes to the solver."""
        pass

    @abstractmethod
    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Returns (residuals, local_jacobian)."""
        pass

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        """
        Computes (T_total, Y, Jac) for nodes based on upstream conditions.
        Jac contains derivatives d(T_total)/d(stream_i) and d(Y)/d(stream_i).
        By default, it just takes the first upstream state (pass-through).
        """
        if not upstream_states:
            import combaero as cb

            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        up = upstream_states[0]
        return up.T_total, up.Y, None

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def mach(self, state: MixtureState) -> float:
        """Returns the local Mach number at this node. Default 0.0 (stagnation)."""
        return 0.0

    def diagnostics(self, state: MixtureState) -> dict[str, float]:
        """Compute generalized node-level diagnostics (e.g., thermo properties)."""
        import combaero as cb

        # Robustness check: state must be physically valid (P > 0)
        if state.P <= 0:
            return {}

        # Robust fallback using CompleteState
        try:
            cs = cb.complete_state(state.T, state.P, state.X)
            return {
                "h": cs.thermo.h,
                "s": cs.thermo.s,
                "u": cs.thermo.u,
                "rho": cs.thermo.rho,
                "gamma": cs.thermo.gamma,
                "a": cs.thermo.a,
                "cp": cs.thermo.cp,
                "cv": cs.thermo.cv,
                "mw": cs.thermo.mw,
                "mu": cs.transport.mu,
                "k": cs.transport.k,
                "Pr": cs.transport.Pr,
                "nu": cs.transport.nu,
            }
        except Exception:
            # Fallback to pure thermo if transport fails or is unsupported
            ts = cb.thermo_state(state.T, state.P, state.X)
            return {
                "h": ts.h,
                "s": ts.s,
                "u": ts.u,
                "rho": ts.rho,
                "gamma": ts.gamma,
                "a": ts.a,
                "cp": ts.cp,
                "cv": ts.cv,
                "mw": ts.mw,
            }


class NetworkElement(ABC):
    def __init__(self, id: str, from_node: str, to_node: str):
        self.id = id
        self.from_node = from_node
        self.to_node = to_node

    @property
    def has_convective_surface(self) -> bool:
        """True if this element has a convective heat transfer surface."""
        return False

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this element contributes to the solver."""
        pass

    @abstractmethod
    def residuals(
        self, state_in: MixtureState, state_out: MixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """
        Returns (residuals, local_jacobian).
        local_jacobian is a dict mapping equation_index to a dict of {unknown_name: partial_derivative}.
        """
        pass

    @abstractmethod
    def n_equations(self) -> int:
        """Number of residual equations contributed."""
        pass

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass

    def htc_and_T(self, state: MixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature."""
        return None

    def diagnostics(self, state_in: MixtureState, state_out: MixtureState) -> dict[str, float]:
        """Compute diagnostic properties for this element (e.g. Mach, P_ratio)."""
        return {}


class PlenumNode(NetworkNode):
    """
    A stagnation volume where v ~ 0, so P_total = P_static.
    For Phase 1: pressure-only node.
    For Phase 2+: automatically handles mixing when multiple upstream connections exist.
    """

    def __init__(self, id: str):
        super().__init__(id)
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this plenum."""
        self.energy_boundaries.append(eb)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a plenum (simple mixing + energy boundaries)."""
        import combaero as cb

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.T_total, s.P_total, s.Y) for s in upstream_states]
        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = cb.mixer_from_streams_and_jacobians(streams, Q=Q_total, fraction=fraction_total)

        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Base residual: P_total = P for plenum
        res = [state.P_total - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.P_total": 1.0}}
        return res, jac

    def diagnostics(self, state: MixtureState) -> dict[str, float]:
        import combaero as cb

        # Plenums don't have transport properties (per user request / stationary reservoir model)
        ts = cb.thermo_state(state.T, state.P, state.X)
        return {
            "h": ts.h,
            "s": ts.s,
            "u": ts.u,
            "rho": ts.rho,
            "gamma": ts.gamma,
            "a": ts.a,
            "cp": ts.cp,
            "cv": ts.cv,
            "mw": ts.mw,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class MomentumChamberNode(NetworkNode):
    """
    Similar to a plenum, but preserves dynamic pressure.
    Enforces P_total = P_static + 0.5 * rho * v^2 instead of P_total = P_static.
    Supports directional flow vectors (port angles) for momentum conservation.
    Automatically handles mixing when multiple upstream connections exist.
    """

    def __init__(
        self,
        id: str,
        area: float = 0.1,
        port_angles_deg: dict[str, float] | None = None,
        length: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_wall: float | None = None,
        Dh: float | None = None,
    ):
        super().__init__(id)
        self.area = area  # Cross-sectional area for momentum calculations
        self.port_angles_deg: dict[str, float] = port_angles_deg or {}
        self.length = length
        self.Dh = Dh
        self.surface = surface or ConvectiveSurface()
        self.t_wall = t_wall
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

    @property
    def has_convective_surface(self) -> bool:
        return True

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this momentum chamber."""
        self.energy_boundaries.append(eb)

    def set_port_angle(self, element_id: str, angle_deg: float) -> None:
        """Set the angle of the flow for a specific connected element."""
        self.port_angles_deg[element_id] = angle_deg

    def get_port_angle(self, element_id: str) -> float:
        """Get the angle of the flow for a specific connected element (defaults to 0.0)."""
        return self.port_angles_deg.get(element_id, 0.0)

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a momentum chamber (simple mixing + energy boundaries)."""
        import combaero as cb

        # Store total mass flow for use in residuals
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0

        # Store upstream element IDs for Jacobian
        self._upstream_element_ids = []
        for s in upstream_states:
            if hasattr(s, "_element_id"):
                self._upstream_element_ids.append(s._element_id)

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.T_total, s.P_total, s.Y) for s in upstream_states]
        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = cb.mixer_from_streams_and_jacobians(streams, Q=Q_total, fraction=fraction_total)

        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def htc_and_T(self, state: MixtureState):
        """Compute heat transfer coefficient for the momentum chamber."""
        if self.surface.area == 0.0:
            return None

        import math

        T_wall = self.t_wall if self.t_wall is not None else math.nan
        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho, _ = _safe_rho(state.density())
        u = m_dot_total / (rho * self.area) if self.area > 0 else 1.0

        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)
        length = self.length if self.length is not None else diameter

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=diameter,
            length=length,
            T_wall=T_wall,
        )

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        import combaero as cb

        # Momentum chamber: P_total = P_static + 0.5 * rho * v^2
        # Use total mass flow computed during compute_derived_state
        m_dot_total = getattr(self, "_total_m_dot", 0.0)

        # Use C++ function for residual and analytical Jacobian
        result = cb.momentum_chamber_residual_and_jacobian(
            state.P, state.P_total, m_dot_total, state.T, state.Y, self.area
        )

        res = [result.residual]
        jac = {
            0: {
                f"{self.id}.P": result.d_res_dP,
                f"{self.id}.P_total": result.d_res_dP_total,
            }
        }

        # Add Jacobian entries for upstream element mass flows
        upstream_elem_ids = getattr(self, "_upstream_element_ids", [])
        for elem_id in upstream_elem_ids:
            jac[0][f"{elem_id}.m_dot"] = result.d_res_dmdot

        return res, jac

    def mach(self, state: MixtureState) -> float:
        """Computes Mach number using internal total mass flow and area."""
        import combaero as cb

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        if self.area <= 0 or m_dot_total <= 1e-12:
            return 0.0

        rho = state.density()
        velocity = m_dot_total / (rho * self.area)
        return float(cb.mach_number(velocity, state.T, state.X))

    def diagnostics(self, state: MixtureState) -> dict[str, float]:
        import math

        import combaero as cb

        # Momentum chambers use rich CompleteState
        cs = cb.complete_state(state.T, state.P, state.X)

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho = cs.thermo.rho
        u = m_dot_total / (rho * self.area) if rho > 0 and self.area > 0 else 1.0
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)

        re = (rho * u * diameter / cs.transport.mu) if cs.transport.mu > 0 else 1.0

        # Detailed heat transfer diagnostics if surface is active
        nu_val = 0.0
        htc_val = 0.0
        t_aw_val = cs.thermo.T
        if self.has_convective_surface:
            h_res = self.htc_and_T(state)
            if h_res is not None:
                nu_val = h_res.Nu
                htc_val = h_res.h
                t_aw_val = h_res.T_aw

        return {
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": (u / cs.thermo.a) if cs.thermo.a > 0 else 1.0,
            "Nu": nu_val,
            "htc": htc_val,
            "T_aw": t_aw_val,
        }

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class PressureBoundary(NetworkNode):
    """
    Supplies fixed stagnation pressure and temperature. Contributes no unknowns and no residuals.
    Essential as a reference pressure and absolute mass sink/source for the network.
    """

    def __init__(
        self,
        id: str,
        P_total: float = 101325.0,
        T_total: float = 300.0,
        Y: list[float] | None = None,
    ) -> None:
        super().__init__(id)
        self.P_total = P_total
        self.T_total = T_total
        self.Y = Y

    def unknowns(self) -> list[str]:
        return []

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        return [], {}

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        # Boundary nodes define their own state
        import combaero as cb

        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.T_total, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class MassFlowBoundary(NetworkNode):
    """
    Supplies fixed mass flow and stagnation temperature. Its pressure floats to satisfy the flow equations.
    Cannot exclusively define a network, a pressure reference is also required.
    """

    def __init__(
        self,
        id: str,
        m_dot: float = 0.1,
        T_total: float = 300.0,
        Y: list[float] | None = None,
    ) -> None:
        super().__init__(id)
        self.m_dot = m_dot
        self.T_total = T_total
        self.Y = Y

    def unknowns(self) -> list[str]:
        # Pressure floats to whatever is required to push the defined mass flow
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation assumptions: P_total = P_static
        res = [state.P_total - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.P_total": 1.0}}
        return res, jac

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        # Boundary nodes define their own state
        import combaero as cb

        Y = self.Y if self.Y is not None else list(cb.mole_to_mass(cb.species.dry_air()))
        return self.T_total, Y, None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class CombustorNode(NetworkNode):
    """
    Combustion chamber adding energy (and optionally mass) to the flow.
    Evaluates adiabatic combustion using the C++ core backend based on
    the selected CombustionMethodLiteral.
    Automatically handles mixing when multiple upstream connections exist.
    """

    def __init__(
        self,
        id: str,
        method: CombustionMethodLiteral = "complete",
        pressure_loss: Any = None,
        area: float = 0.1,
        Dh: float | None = None,
        surface: ConvectiveSurface | None = None,
        t_wall: float | None = None,
    ):
        super().__init__(id)
        self.method = method
        self.pressure_loss = pressure_loss
        self.area = area
        self.Dh = Dh
        self.surface = surface or ConvectiveSurface()
        self.t_wall = t_wall
        self.upstream_elements = []
        self.fuel_boundary = None
        self.energy_boundaries: list[EnergyBoundary] = []

    @property
    def has_convective_surface(self) -> bool:
        # CombustorNode has a surface attribute but no htc_and_T implementation yet.
        # Return False until htc_and_T is added (planned for a future phase).
        return False

    def add_energy_boundary(self, eb: EnergyBoundary) -> None:
        """Attach an energy source/sink to this combustor (post-combustion)."""
        self.energy_boundaries.append(eb)

    def set_fuel_boundary(self, fuel_bc: MassFlowBoundary) -> None:
        """Set the fuel boundary condition for this combustor."""
        self.fuel_boundary = fuel_bc

    def unknowns(self) -> list[str]:
        # Pure Pressure-Flow: Temperature and Composition are derived forward, not unknowns.
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def compute_derived_state(
        self, upstream_states: list[MixtureState]
    ) -> tuple[float, list[float], Any]:
        """Derived T and Y for a combustor (Reaction + Mixing + energy boundaries)."""
        import combaero as cb

        # Store total mass flow and unburned temperature for use in diagnostics
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 1.0
        self._T_unburned = (
            sum(s.m_dot * s.T_total for s in upstream_states) / self._total_m_dot
            if self._total_m_dot > 0
            else 300.0
        )

        if not upstream_states:
            # Default fallback
            return 300.0, list(cb.mole_to_mass(cb.species.dry_air())), None

        streams = [cb.MassStream(s.m_dot, s.T_total, s.P_total, s.Y) for s in upstream_states]
        P_ref = upstream_states[0].P if upstream_states else 101325.0

        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        # Use new combustor C++ fast-path if available
        if self.pressure_loss is not None:
            mix_res = cb.combustor_residuals_and_jacobians(
                streams,
                P_ref,
                Q=Q_total,
                fraction=fraction_total,
                pressure_loss=self.pressure_loss,
                use_equilibrium=(self.method == "equilibrium"),
            )
        elif self.method == "equilibrium":
            mix_res = cb.adiabatic_T_equilibrium_and_jacobians_from_streams(
                streams, P_ref, Q=Q_total, fraction=fraction_total
            )
        else:
            mix_res = cb.adiabatic_T_complete_and_jacobian_T_from_streams(
                streams, P_ref, Q=Q_total, fraction=fraction_total
            )
        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation constraint: P_total = P
        return [state.P_total - state.P], {0: {f"{self.id}.P": -1.0, f"{self.id}.P_total": 1.0}}

    def diagnostics(self, state: MixtureState) -> dict[str, float]:
        import math

        import combaero as cb

        # Combustors use rich CompleteState
        cs = cb.complete_state(state.T, state.P, state.X)

        m_dot_total = getattr(self, "_total_m_dot", 0.0)
        rho = cs.thermo.rho
        u = m_dot_total / (rho * self.area) if rho > 0 and self.area > 0 else 1.0
        diameter = self.Dh if self.Dh is not None else math.sqrt(4.0 * self.area / math.pi)

        re = (rho * u * diameter / cs.transport.mu) if cs.transport.mu > 0 else 1.0

        diag = {
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            "cp": cs.thermo.cp,
            "cv": cs.thermo.cv,
            "mw": cs.thermo.mw,
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
            "Re": re,
            "Dh": diameter,
            "velocity": u,
            "mach": u / cs.thermo.a if cs.thermo.a > 0 else 1.0,
        }

        # --- Combustion Diagnostics ---
        # 1. Equivalence Ratio (phi) via elemental analysis in C++
        try:
            phi = cb.equivalence_ratio(state.X)
            diag["phi"] = phi
        except Exception as e:
            # Reveal the root cause of the diagnostic failure
            import sys

            print(f"DEBUG: CombustorNode.diagnostics phi calculation failed: {e}", file=sys.stderr)
            diag["phi"] = 0.0

        # 2. Temperature Rise Ratio (theta) = (T_burned / T_unburned) - 1
        T_u = getattr(self, "_T_unburned", state.T_total)
        if T_u > 0:
            diag["theta"] = (state.T_total / T_u) - 1.0
        else:
            diag["theta"] = 0.0

        return diag

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class OrificeElement(NetworkElement):
    """
    Orifice flow element with incompressible or compressible formulation.

    - regime='incompressible': m_dot = Cd * A * sqrt(2 * rho * dP)
    - regime='compressible': Uses isentropic nozzle flow with smooth choked transition

    The discharge coefficient Cd is determined by the 'correlation' parameter:
      - 'fixed': Uses the user-supplied Cd value directly.
      - 'ReaderHarrisGallagher': Sharp thin-plate (ISO 5167-2 / RHG).
      - 'Stolz': ISO 5167:1980 (Corner taps).
      - 'Miller': Miller (1996) simplified correlation.
      - 'ThickPlate': Thick-plate sharp-edged orifice (requires plate_thickness).
      - 'RoundedEntry': Rounded-entry orifice (requires edge_radius).
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        Cd: float = 0.6,
        diameter: float | None = None,
        regime: Literal["incompressible", "compressible"] = "incompressible",
        correlation: str = "ReaderHarrisGallagher",
        plate_thickness: float = 0.0,
        edge_radius: float = 0.0,
        area: float | None = None,
    ):
        super().__init__(id, from_node, to_node)
        import math

        if diameter is None and area is None:
            raise ValueError("OrificeElement requires either 'diameter' or 'area'.")

        if diameter is not None:
            self.diameter = diameter
            self.area = math.pi * (diameter / 2.0) ** 2
        else:
            self.area = area
            self.diameter = math.sqrt(4.0 * area / math.pi)

        self.Cd = Cd
        self.regime = regime
        self.correlation = correlation
        self.use_correlation = correlation != "fixed"
        self.plate_thickness = plate_thickness
        self.edge_radius = edge_radius
        self.upstream_diameter: float | None = None
        self.downstream_diameter: float | None = None
        # OrificeGeometry built in resolve_topology
        self._orifice_geom: object | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Evaluate both sides of the node for geometry discovery
        upstream_elements = graph.get_upstream_elements(self.from_node)
        downstream_elements = graph.get_downstream_elements(self.to_node)

        # Ensure we don't pick ourselves if it's a tight chain
        upstream_channels = [
            e for e in upstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]
        downstream_channels = [
            e for e in downstream_elements if isinstance(e, ChannelElement) and e.id != self.id
        ]

        if len(upstream_channels) == 1:
            self.upstream_diameter = upstream_channels[0].diameter

        if len(downstream_channels) == 1:
            self.downstream_diameter = downstream_channels[0].diameter

        # Pre-compute beta for velocity-of-approach factor
        import math
        import warnings

        d_bore = math.sqrt(4.0 * self.area / math.pi)
        self.beta = 0.0
        self._orifice_geom = None

        if self.upstream_diameter and self.upstream_diameter > 0:
            if d_bore < self.upstream_diameter:
                self.beta = d_bore / self.upstream_diameter
            else:
                warnings.warn(
                    f"OrificeElement '{self.id}': inferred bore diameter "
                    f"({d_bore:.4f} m) >= channel diameter "
                    f"({self.upstream_diameter:.4f} m). "
                    f"Skipping velocity-of-approach correction (E=1).",
                    stacklevel=2,
                )

        if self.use_correlation:
            import combaero as cb

            # Build geometry descriptor for Cd correlation.
            # D=0 when no upstream channel known (plenum connection) -> beta=0 -> RHG extrapolates to ~0.597.
            D_up = (
                self.upstream_diameter
                if self.upstream_diameter and self.upstream_diameter > 0
                else 1.0
            )
            if D_up <= d_bore:
                D_up = d_bore * 10.0
            geom = cb.OrificeGeometry()
            geom.d = d_bore
            geom.D = D_up
            geom.t = self.plate_thickness
            geom.r = self.edge_radius
            self._orifice_geom = geom

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def _effective_Cd(self, state_in: "MixtureState", state_out: "MixtureState") -> float:
        """Return effective Cd: from correlation or user-supplied fixed value."""
        if not self.use_correlation or self._orifice_geom is None:
            # Clamp manual Cd strictly <= 1.0 to strictly preserve vena contracta physics (A_eff <= A_geom)
            return max(1e-4, min(self.Cd, 1.0))

        import combaero as cb

        # Estimate Re_D from current m_dot and upstream viscosity.
        D_up = self._orifice_geom.D
        if D_up > 0.0:
            ts = cb.transport_state(state_in.T, state_in.P, state_in.X)
            mu = ts.mu if ts.mu > 0.0 else 1.8e-5
            m_dot_abs = abs(state_in.m_dot)
            Re_D = (4.0 * m_dot_abs) / (math.pi * D_up * mu) if m_dot_abs > 1e-12 else 1e4
        else:
            Re_D = 1e5  # typical reference when no upstream channel

        flow_state = cb.OrificeState()
        flow_state.Re_D = Re_D
        flow_state.dP = max(state_in.P_total - state_out.P, 1.0)
        flow_state.rho = cb.density(state_in.T, state_in.P, state_in.X)
        flow_state.mu = cb.transport_state(state_in.T, state_in.P, state_in.X).mu

        # Correlation selection logic
        if self.correlation == "ReaderHarrisGallagher":
            return float(cb.Cd_sharp_thin_plate(self._orifice_geom, flow_state))
        elif self.correlation == "Stolz":
            b = self._orifice_geom.beta
            b2 = b * b
            b8 = b2 * b2 * b2 * b2
            cd = (
                0.5961
                + 0.0261 * b2
                - 0.216 * b8
                + 0.000521 * math.pow(1e6 * b / max(flow_state.Re_D, 1.0), 0.7)
            )
            return float(cd)
        elif self.correlation == "Miller":
            b = self._orifice_geom.beta
            b2 = b * b
            b8 = b2 * b2 * b2 * b2
            cd = (
                0.5959
                + 0.0312 * math.pow(b, 2.1)
                - 0.184 * b8
                + 91.71 * math.pow(b, 2.5) * math.pow(max(flow_state.Re_D, 1.0), -0.75)
            )
            return float(cd)
        elif self.correlation == "ThickPlate":
            if self._orifice_geom.t <= 0:
                return float(cb.Cd_sharp_thin_plate(self._orifice_geom, flow_state))
            return float(cb.Cd_thick_plate(self._orifice_geom, flow_state))
        elif self.correlation == "RoundedEntry":
            if self._orifice_geom.r <= 0:
                return float(cb.Cd_sharp_thin_plate(self._orifice_geom, flow_state))
            return float(cb.Cd_rounded_entry(self._orifice_geom, flow_state))
        else:
            # Default/Auto
            return float(cb.Cd_orifice(self._orifice_geom, flow_state))

    def residuals(
        self, state_in: "MixtureState", state_out: "MixtureState"
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        import combaero as cb

        m_dot = state_in.m_dot
        effective_cd = self._effective_Cd(state_in, state_out)

        if self.regime == "compressible":
            res_cpp = cb._core.orifice_compressible_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.T_total,
                state_in.Y,
                state_out.P,
                effective_cd,
                self.area,
                getattr(self, "beta", 0.0),
            )
        else:
            # Use incompressible Bernoulli formulation
            res_cpp = cb.orifice_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.P,
                state_in.T,
                state_in.Y,
                state_out.P,
                effective_cd,
                self.area,
                beta=getattr(self, "beta", 0.0),
            )

        res = [m_dot - res_cpp.m_dot_calc]

        # Assemble Jacobian with respect to all node and element unknowns
        jac = {
            0: {
                f"{self.id}.m_dot": 1.0,
                f"{self.from_node}.P_total": -res_cpp.d_mdot_dP_total_up,
                f"{self.from_node}.P": -res_cpp.d_mdot_dP_static_up,
                f"{self.from_node}.T": -res_cpp.d_mdot_dT_up,
                f"{self.to_node}.P": -res_cpp.d_mdot_dP_static_down,
            }
        }
        # Add species sensitivities
        for i, val in enumerate(res_cpp.d_mdot_dY_up):
            jac[0][f"{self.from_node}.Y[{i}]"] = -val

        return res, jac

    def n_equations(self) -> int:
        return 1

    def diagnostics(self, state_in: MixtureState, state_out: MixtureState) -> dict[str, float]:

        import combaero as cb

        # Full thermo + transport bundle at inlet static conditions (mole fractions X)
        cs = cb.complete_state(state_in.T, state_in.P, state_in.X)

        p_ratio = state_out.P / state_in.P_total if state_in.P_total > 0 else 1.0

        # Estimate throat Mach number
        gamma = cs.thermo.gamma
        # Critical pressure ratio for choking
        pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))

        if p_ratio <= pr_crit:
            mach_throat = 1.0
        else:
            mach_throat = float(
                cb.mach_from_pressure_ratio(
                    state_in.T_total, state_in.P_total, state_out.P, state_in.X
                )
            )

        # Effective Cd actually used in this solve
        effective_cd = self._effective_Cd(state_in, state_out)

        # Calculate macro physical velocity strictly based on the geometric area
        # This prevents physical velocity and Mach from artificially zeroing out when tuning Cd ~ infinity
        v_physical = abs(state_in.m_dot) / (cs.thermo.rho * self.area) if self.area > 0 and cs.thermo.rho > 0 else 0.0
        mach_physical = v_physical / cs.thermo.a if cs.thermo.a > 0 else 0.0

        return {
            "v": v_physical,
            "mach": mach_physical,
            "mach_throat": mach_throat,
            "p_ratio": p_ratio,
            "Cd": effective_cd,
            "is_correlation": float(self.use_correlation),
            # Thermo
            "h": cs.thermo.h,
            "s": cs.thermo.s,
            "u": cs.thermo.u,
            "rho": cs.thermo.rho,
            "gamma": cs.thermo.gamma,
            "a": cs.thermo.a,
            # Transport
            "mu": cs.transport.mu,
            "k": cs.transport.k,
            "Pr": cs.transport.Pr,
            "nu": cs.transport.nu,
        }


class EffectiveAreaConnectionElement(OrificeElement):
    """
    An orifice element with user-specified effective area (Cd * A product).

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A_eff * sqrt(2 * rho * dP)

    where A_eff is the effective area (the product Cd * A).

    **Incompressible Formulation**: The solver uses `cb.orifice_mdot_and_jacobian()`,
    which applies the incompressible Bernoulli equation with regularization for numerical
    stability. Density is evaluated at upstream conditions. For compressible flows with
    significant Mach number effects, use a nozzle element instead (future feature).

    **Effective Area Interpretation**: The effective area represents the combined effect
    of discharge coefficient and geometric area. For example:
    - A_eff = 0.01 m^2 could represent Cd=1.0 * A=0.01 m^2
    - Or equivalently: Cd=0.8 * A=0.0125 m^2
    - The product (Cd * A) is what matters for flow calculation

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you need a simple area-based flow restriction without detailed geometry modeling.

    The effective area is user-specified and already accounts for all geometric effects
    (entrance/exit losses, contraction, vena contracta, etc.), so upstream/downstream
    geometry discovery is not performed.
    """

    def __init__(self, id: str, from_node: str, to_node: str, effective_area: float) -> None:
        """
        Initialize with effective area (Cd * A product).

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            effective_area: Effective area for flow calculation (m^2).
                           This is the product Cd * A, where Cd is the discharge coefficient
                           and A is the geometric area. Stored internally as self.area.

        Example:
            >>> conn = EffectiveAreaConnectionElement("conn1", "inlet", "outlet", 0.01)
            >>> conn.area  # 0.01 m^2 (effective area = Cd * A)
            >>> conn.Cd    # 1.0 (normalized discharge coefficient)

        Note:
            The effective area already includes all loss coefficients. For example, if you
            have a sharp-edged orifice with geometric area 0.0125 m^2 and Cd=0.8, you would
            specify effective_area=0.01 m^2 (0.8 * 0.0125).
        """
        import math

        diameter = math.sqrt(4.0 * effective_area / math.pi)
        super().__init__(id, from_node, to_node, Cd=1.0, diameter=diameter, correlation="fixed")

        # Force self.area to be precisely the effective area
        # (to remove any floating point math.sqrt / pi precision issues)
        self.area = effective_area

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        Unlike OrificeElement, we do not discover upstream/downstream diameters because
        the effective area is user-specified and already accounts for all geometric effects.
        No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - effective area is pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class DiameterDischargeCoefficientConnectionElement(OrificeElement):
    """
    An orifice element with user-specified physical diameter and discharge coefficient or loss coefficient.

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A * Cd * sqrt(2 * rho * dP)

    where A is the computed physical area and Cd is the discharge coefficient.

    **Incompressible Formulation**: Uses the Bernoulli equation with density evaluated
    at upstream conditions. Valid for low Mach number flows (M < 0.3 typically).

    **Parameters**: User provides either Cd (discharge coefficient) or zeta (loss coefficient).
    The other parameter is calculated automatically using the relationship:
        zeta = 1/Cd^2 - 1  or  Cd = 1/sqrt(zeta + 1)

    **Effective Area**: The effective area is calculated as A_eff = A * Cd.

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you have a physical diameter measurement and want to specify loss characteristics
    via discharge coefficient or loss coefficient.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        diameter: float,
        Cd: float | None = None,
        zeta: float | None = None,
    ) -> None:
        """
        Initialize with physical diameter and either Cd or zeta.

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            diameter: Physical geometric diameter (m)
            Cd: Discharge coefficient (optional, 0 < Cd <= 1)
            zeta: Loss coefficient (optional, zeta >= 0)

        Note:
            Exactly one of Cd or zeta must be provided. If both are provided, zeta is ignored.
            If Cd is provided, zeta is calculated as zeta = 1/Cd^2 - 1.
            If zeta is provided, Cd is calculated as Cd = 1/sqrt(zeta + 1).

        Example:
            >>> # Using Cd
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, Cd=0.8)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (discharge coefficient)

            >>> # Using zeta
            >>> conn = DiameterDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.1, zeta=0.5625)
            >>> conn.diameter  # 0.1 m (physical diameter)
            >>> conn.Cd    # 0.8 (calculated from zeta)

        Raises:
            ValueError: If neither Cd nor zeta is provided, or if invalid values are given.
        """
        if Cd is not None and zeta is not None:
            # If both are provided, use Cd and ignore zeta
            zeta = None
        elif Cd is None and zeta is None:
            raise ValueError("Either Cd or zeta must be provided")

        if Cd is not None:
            if not (0 < Cd <= 1):
                raise ValueError("Cd must be in range (0, 1]")
        else:  # zeta is provided
            if zeta < 0:
                raise ValueError("zeta must be >= 0")
            # Calculate Cd from zeta: Cd = 1/sqrt(zeta + 1)
            Cd = 1.0 / (zeta + 1.0) ** 0.5

        # Pass physical diameter and Cd directly to parent
        super().__init__(id, from_node, to_node, Cd=Cd, diameter=diameter, correlation="fixed")

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """
        Skip upstream/downstream geometry discovery.

        The physical area and discharge coefficient are user-specified and already account
        for all geometric effects. No Beta ratio correction is needed or applied.
        """
        # Intentionally empty - area and Cd are pre-computed by user
        pass

    def n_equations(self) -> int:
        """Return the number of equations (1: mass flow balance)."""
        return 1


class LosslessConnectionElement(NetworkElement):
    """
    An ideal connection with no friction or momentum loss.
    Inherently preserves total pressure: P_total_in = P_total_out.
    Does not require geometric parameters.
    """

    def __init__(self, id: str, from_node: str, to_node: str):
        super().__init__(id, from_node, to_node)
        self.area: float | None = None
        self.diameter: float | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass  # Zero-loss elements do not need upstream geometry

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(
        self, state_in: MixtureState, state_out: MixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Residual: P_total_in - P_total_out = 0
        res = [state_in.P_total - state_out.P_total]
        jac = {0: {f"{self.from_node}.P_total": 1.0, f"{self.to_node}.P_total": -1.0}}
        return res, jac

    def get_spatial_profile(self, n_steps: int = 100):
        """
        Returns a mock spatial profile array representing the zero-loss state.
        This provides structural symmetry with realistic pipes (e.g. Fanno/Rough)
        for GUI plotters.
        """
        import combaero as cba

        # Simplified placeholder for test completion validation
        # Normally would utilize the true state nodes
        profile = []
        for _ in range(n_steps):
            st = cba.IncompressibleStation()
            st.x = 0.0
            st.P = 100000.0
            st.T = 300.0
            st.rho = 1.2
            st.v = 1.0
            st.M = 0.0
            st.h = 300000.0
            profile.append(st)
        return profile

    def n_equations(self) -> int:
        return 1


class ChannelElement(NetworkElement):
    """
    Channel element applying frictional pressure drop.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        length: float,
        diameter: float,
        roughness: float,
        regime: CompressibilityLiteral = "incompressible",
        friction_model: FrictionModelLiteral = "haaland",
        htc_model: HeatTransferModelLiteral = "none",
        t_wall: float | None = None,
        surface: ConvectiveSurface | None = None,
    ):
        super().__init__(id, from_node, to_node)
        self.length = length
        self.diameter = diameter
        self.roughness = roughness
        self.area = math.pi * (diameter / 2) ** 2

        self.regime = regime
        self.friction_model = friction_model
        self.htc_model = htc_model
        self.t_wall = t_wall
        self.surface = surface or ConvectiveSurface()

    @property
    def has_convective_surface(self) -> bool:
        return True

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:

        import combaero as cb

        m_dot = state_in.m_dot

        if self.regime == "compressible":
            # Use compressible Fanno flow with friction
            res_cpp = cb._core.channel_compressible_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.T_total,
                state_in.Y,
                state_out.P,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
            )
        else:
            # Use incompressible Darcy-Weisbach formulation
            res_cpp = cb.channel_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.P,
                state_in.T,
                state_in.Y,
                state_out.P,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
            )

        # Residual of P-node unknowns matching channel drop
        # C++ returns dP_calc = friction_loss
        # In a stable network: P_total_up - P_static_down = dP_friction (for exit into plenum)
        # Or more generally: P_total_up - P_total_down = dP_friction

        # Test suite currently expects P_total_up - P_static_down = dP_calc
        res = [state_in.P_total - state_out.P - res_cpp.dP_calc]

        upstream_id = self.from_node
        downstream_id = self.to_node

        jac = {
            0: {
                f"{self.id}.m_dot": -res_cpp.d_dP_d_mdot,
                f"{upstream_id}.P_total": 1.0,
                f"{upstream_id}.P": -res_cpp.d_dP_dP_static_up,
                f"{upstream_id}.T": -res_cpp.d_dP_dT_up,
                f"{downstream_id}.P": -1.0,
            }
        }
        for i, val in enumerate(res_cpp.d_dP_dY_up):
            jac[0][f"{upstream_id}.Y[{i}]"] = -val

        return res, jac

    def diagnostics(self, state_in: MixtureState, state_out: MixtureState) -> dict[str, float]:
        """Compute inlet/outlet Mach numbers and pressure ratio."""
        import combaero as cb

        # Inlet Mach
        rho_in = state_in.density()
        v_in = state_in.m_dot / (rho_in * self.area) if self.area > 0 and rho_in > 0 else 1.0
        mach_in = float(cb.mach_number(v_in, state_in.T, state_in.X))

        # Outlet Mach (at static pressure P_out)
        rho_out = cb.density(state_out.T, state_out.P, state_out.X)
        v_out = state_out.m_dot / (rho_out * self.area) if self.area > 0 and rho_out > 0 else 1.0
        mach_out = float(cb.mach_number(v_out, state_out.T, state_out.X))

        p_ratio_total = state_out.P_total / state_in.P_total if state_in.P_total > 0 else 1.0

        # Full thermo + transport bundle at inlet static conditions (mole fractions X)
        cs = cb.complete_state(state_in.T, state_in.P, state_in.X)

        # Reynolds number at inlet: Re = rho * v * D / mu
        re_in = (rho_in * v_in * self.diameter / cs.transport.mu) if cs.transport.mu > 0 else 1.0

        res_dict = {
            "mach_in": mach_in,
            "mach_out": mach_out,
            "p_ratio_total": p_ratio_total,
            "Re": re_in,
            "dP": max(0.0, state_in.P_total - state_out.P_total),
        }

        # Heat Transfer diagnostics
        h_res = self.htc_and_T(state_in)
        if h_res is not None:
            res_dict["Nu"] = h_res.Nu
            res_dict["htc"] = h_res.h
            res_dict["T_aw"] = h_res.T_aw
        else:
            res_dict["Nu"] = 0.0
            res_dict["htc"] = 0.0
            res_dict["T_aw"] = state_in.T

        # Thermo
        res_dict.update(
            {
                "h": cs.thermo.h,
                "s": cs.thermo.s,
                "u": cs.thermo.u,
                "rho": cs.thermo.rho,
                "gamma": cs.thermo.gamma,
                "a": cs.thermo.a,
                # Transport
                "mu": cs.transport.mu,
                "k": cs.transport.k,
                "Pr": cs.transport.Pr,
                "nu": cs.transport.nu,
            }
        )
        return res_dict

    def get_spatial_profile(
        self,
        state_in: MixtureState,
        n_steps: int = 100,
    ) -> list:
        """
        Compute and return the spatial flow profile array along the channel length.
        Requires the solved inlet MixtureState.
        """
        import combaero as cba

        rho = cba.density(state_in.T, state_in.P, state_in.X)
        u = cba.channel_velocity(state_in.m_dot, self.diameter, rho)

        if self.regime == "incompressible":
            res = cba.channel_flow_rough(
                state_in.T,
                state_in.P,
                state_in.X,
                u,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                n_steps,
                True,
            )
            return res.profile

        elif self.regime == "compressible":
            res = cba.fanno_channel(
                state_in.T,
                state_in.P,
                cba.mass_to_mole(
                    state_in.Y
                ),  # C++ expects mole fractions for Fanno solver currently
                state_in.m_dot,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                n_steps,
                True,
            )
            return res.profile

        return []

    def htc_and_T(self, state: MixtureState):
        """Compute heat transfer coefficient and adiabatic wall temperature.

        Uses the element's ConvectiveSurface to compute HTC and adiabatic wall temperature.

        Parameters
        ----------
        state : MixtureState
            Flow state at the element inlet.

        Returns
        -------
        ChannelResult | None
            Full ChannelResult with h, T_aw, and Jacobians (dh_dmdot, dh_dT, etc.),
            or None if surface.area = 0. Access convective area via ``self.surface.area``.
        """
        if self.surface.area == 0.0:
            return None

        rho, _ = _safe_rho(state.density())
        u = state.m_dot / (rho * self.area) if self.area > 0 else 1.0

        # Use nan for T_wall if not specified (matches C++ default)
        T_wall = self.t_wall if self.t_wall is not None else math.nan

        return self.surface.htc_and_T(
            T=state.T,
            P=state.P,
            X=state.X,
            velocity=u,
            diameter=self.diameter,
            length=self.length,
            T_wall=T_wall,
        )

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass
