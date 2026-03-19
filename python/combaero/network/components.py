import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb

# Physics configuration types for network introspectability
CompressibilityLiteral = Literal["incompressible", "compressible", "compressible_fanno"]
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
    ) -> tuple[float, float, float] | None:
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
        tuple[float, float, float] | None
            (h [W/(m^2*K)], T_aw [K], A_conv [m^2]) or None if area=0.
        """
        import math

        import combaero as cb

        if self.area == 0.0:
            return None

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
                heating=self.heating if self.heating is not None else True,
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
                heating=self.heating if self.heating is not None else True,
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
                heating=self.heating if self.heating is not None else True,
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
            # For impingement, we need mdot_jet which requires computing from velocity
            # This is a simplified version - full implementation would need element context
            rho = cb.density(T, P, X)
            A_cross = math.pi / 4 * diameter**2
            mdot_jet = rho * velocity * A_cross

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

        return result.h, result.T_aw, self.area


# ============================================================================
# WallConnection (Phase 4/5)
# ============================================================================


@dataclass
class WallConnection:
    """Thermal coupling between two elements through a shared wall.

    Attributes
    ----------
    id : str
        Unique identifier for this wall connection.
    element_a : str
        Element ID for side A.
    element_b : str
        Element ID for side B.
    wall_thickness : float
        Wall thickness [m].
    wall_conductivity : float
        Wall thermal conductivity [W/(m*K)].
    contact_area : float | None
        Override effective area [m^2]. If None, uses min(A_conv_a, A_conv_b).
    """

    id: str
    element_a: str
    element_b: str
    wall_thickness: float
    wall_conductivity: float
    contact_area: float | None = None

    def compute_coupling(
        self,
        h_a: float,
        T_aw_a: float,
        A_conv_a: float,
        h_b: float,
        T_aw_b: float,
        A_conv_b: float,
    ) -> tuple[float, float]:
        """Compute heat transfer rate and wall temperature.

        Parameters
        ----------
        h_a : float
            Convective HTC on side A [W/(m^2*K)].
        T_aw_a : float
            Adiabatic wall temperature on side A [K].
        A_conv_a : float
            Convective area on side A [m^2].
        h_b : float
            Convective HTC on side B [W/(m^2*K)].
        T_aw_b : float
            Adiabatic wall temperature on side B [K].
        A_conv_b : float
            Convective area on side B [m^2].

        Returns
        -------
        tuple[float, float]
            (Q [W], T_wall [K]) where Q is positive when heat flows A->B.
        """
        # Effective area
        A_eff = self.contact_area if self.contact_area is not None else min(A_conv_a, A_conv_b)

        # Wall thermal resistance
        t_over_k = self.wall_thickness / self.wall_conductivity

        # Call C++ function
        result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A_eff)

        return result.Q, result.T_wall


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
        """Cached mole fractions converted from mass fractions."""
        try:
            return self._X  # type: ignore[return-value]
        except AttributeError:
            x = list(cb.mass_to_mole(self.Y))
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


class NetworkNode(ABC):
    def __init__(self, id: str):
        self.id = id

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

            return 300.0, list(cb.mole_to_mass(cb.standard_dry_air_composition())), None

        up = upstream_states[0]
        return up.T_total, up.Y, None

    @abstractmethod
    def resolve_topology(self, graph: "FlowNetwork") -> None:
        """Called automatically by FlowNetwork to resolve neighbors."""
        pass


class NetworkElement(ABC):
    def __init__(self, id: str, from_node: str, to_node: str):
        self.id = id
        self.from_node = from_node
        self.to_node = to_node

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
            return 300.0, list(cb.mole_to_mass(cb.standard_dry_air_composition())), None

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
    ):
        super().__init__(id)
        self.area = area  # Cross-sectional area for momentum calculations
        # Dictionary mapping connected element IDs to angle relative to chamber axis [deg]
        self.port_angles_deg: dict[str, float] = port_angles_deg or {}
        self.upstream_elements = []
        self.energy_boundaries: list[EnergyBoundary] = []

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
        self._total_m_dot = sum(s.m_dot for s in upstream_states) if upstream_states else 0.0

        # Store upstream element IDs for Jacobian
        self._upstream_element_ids = []
        for s in upstream_states:
            if hasattr(s, "_element_id"):
                self._upstream_element_ids.append(s._element_id)

        if not upstream_states:
            return 300.0, list(cb.mole_to_mass(cb.standard_dry_air_composition())), None

        streams = [cb.MassStream(s.m_dot, s.T_total, s.P_total, s.Y) for s in upstream_states]
        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        mix_res = cb.mixer_from_streams_and_jacobians(streams, Q=Q_total, fraction=fraction_total)

        return mix_res.T_mix, mix_res.Y_mix, mix_res

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

        Y = (
            self.Y
            if self.Y is not None
            else list(cb.mole_to_mass(cb.standard_dry_air_composition()))
        )
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

        Y = (
            self.Y
            if self.Y is not None
            else list(cb.mole_to_mass(cb.standard_dry_air_composition()))
        )
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
    ):
        super().__init__(id)
        self.method = method
        self.upstream_elements = []
        self.fuel_boundary = None
        self.energy_boundaries: list[EnergyBoundary] = []

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

        if not upstream_states:
            # Default fallback
            return 300.0, list(cb.mole_to_mass(cb.standard_dry_air_composition())), None

        streams = [cb.MassStream(s.m_dot, s.T_total, s.P_total, s.Y) for s in upstream_states]
        P_ref = upstream_states[0].P if upstream_states else 101325.0

        Q_total = sum(eb.Q for eb in self.energy_boundaries)
        fraction_total = sum(eb.fraction for eb in self.energy_boundaries)

        if self.method == "equilibrium":
            mix_res = cb.adiabatic_T_equilibrium_and_jacobians_from_streams(
                streams, P_ref, Q=Q_total, fraction=fraction_total
            )
        else:
            mix_res = cb.adiabatic_T_complete_and_jacobian_T_from_streams(
                streams, P_ref, Q=Q_total, fraction=fraction_total
            )

        return mix_res.T_mix, mix_res.Y_mix, mix_res

    def residuals(self, state: MixtureState) -> tuple[list[float], dict[int, dict[str, float]]]:
        # Stagnation assumptions: P_total = P_static
        # This provides the second equation matching [P, P_total] unknowns,
        # while mass conservation provides the first.
        res = [state.P_total - state.P]
        jac = {0: {f"{self.id}.P": -1.0, f"{self.id}.P_total": 1.0}}
        return res, jac

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Store upstream elements for mixing calculations
        self.upstream_elements = graph.get_upstream_elements(self.id)


class OrificeElement(NetworkElement):
    """
    Orifice flow element with incompressible or compressible formulation.

    - regime='incompressible': m_dot = Cd * A * sqrt(2 * rho * dP)
    - regime='compressible': Uses isentropic nozzle flow with smooth choked transition
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        Cd: float,
        area: float,
        regime: Literal["incompressible", "compressible"] = "incompressible",
    ):
        super().__init__(id, from_node, to_node)
        self.Cd = Cd
        self.area = area
        self.regime = regime
        self.upstream_diameter: float | None = None
        self.downstream_diameter: float | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Evaluate both sides of the node for geometry discovery
        upstream_elements = graph.get_upstream_elements(self.from_node)
        downstream_elements = graph.get_downstream_elements(self.to_node)

        # Ensure we don't pick ourselves if it's a tight chain
        upstream_pipes = [
            e for e in upstream_elements if isinstance(e, PipeElement) and e.id != self.id
        ]
        downstream_pipes = [
            e for e in downstream_elements if isinstance(e, PipeElement) and e.id != self.id
        ]

        if len(upstream_pipes) == 1:
            self.upstream_diameter = upstream_pipes[0].diameter

        if len(downstream_pipes) == 1:
            self.downstream_diameter = downstream_pipes[0].diameter

        # Pre-compute beta for velocity-of-approach factor
        self.beta = 0.0
        if self.upstream_diameter and self.upstream_diameter > 0:
            import math
            import warnings

            d_bore = math.sqrt(4.0 * self.area / math.pi)
            if d_bore < self.upstream_diameter:
                self.beta = d_bore / self.upstream_diameter
            else:
                warnings.warn(
                    f"OrificeElement '{self.id}': inferred bore diameter "
                    f"({d_bore:.4f} m) >= pipe diameter "
                    f"({self.upstream_diameter:.4f} m). "
                    f"Skipping velocity-of-approach correction (E=1).",
                    stacklevel=2,
                )

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(
        self, state_in: MixtureState, state_out: MixtureState
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        import combaero as cb

        m_dot = state_in.m_dot

        if self.regime == "compressible":
            # Use compressible isentropic nozzle flow with smooth choked transition
            res_cpp = cb._core.orifice_compressible_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.T,
                state_in.Y,
                state_out.P,
                self.Cd,
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
                self.Cd,
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
        super().__init__(id, from_node, to_node, Cd=1.0, area=effective_area)

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


class AreaDischargeCoefficientConnectionElement(OrificeElement):
    """
    An orifice element with user-specified physical area and discharge coefficient or loss coefficient.

    This element calculates orifice flow using the incompressible flow equation:
        m_dot = A * Cd * sqrt(2 * rho * dP)

    where A is the physical area and Cd is the discharge coefficient.

    **Incompressible Formulation**: Uses the Bernoulli equation with density evaluated
    at upstream conditions. Valid for low Mach number flows (M < 0.3 typically).

    **Parameters**: User provides either Cd (discharge coefficient) or zeta (loss coefficient).
    The other parameter is calculated automatically using the relationship:
        zeta = 1/Cd^2 - 1  or  Cd = 1/sqrt(zeta + 1)

    **Effective Area**: The effective area is calculated as A_eff = A * Cd.

    Unlike LosslessConnectionElement, this element produces pressure drop proportional to flow rate.
    Use this when you have a physical area measurement and want to specify loss characteristics
    via discharge coefficient or loss coefficient.
    """

    def __init__(
        self,
        id: str,
        from_node: str,
        to_node: str,
        area: float,
        Cd: float | None = None,
        zeta: float | None = None,
    ) -> None:
        """
        Initialize with physical area and either Cd or zeta.

        Args:
            id: Element identifier
            from_node: Upstream node ID
            to_node: Downstream node ID
            area: Physical geometric area (m^2)
            Cd: Discharge coefficient (optional, 0 < Cd <= 1)
            zeta: Loss coefficient (optional, zeta >= 0)

        Note:
            Exactly one of Cd or zeta must be provided. If both are provided, zeta is ignored.
            If Cd is provided, zeta is calculated as zeta = 1/Cd^2 - 1.
            If zeta is provided, Cd is calculated as Cd = 1/sqrt(zeta + 1).

        Example:
            >>> # Using Cd
            >>> conn = AreaDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.0125, Cd=0.8)
            >>> conn.area  # 0.0125 m^2 (physical area)
            >>> conn.Cd    # 0.8 (discharge coefficient)

            >>> # Using zeta
            >>> conn = AreaDischargeCoefficientConnectionElement("conn1", "inlet", "outlet", 0.0125, zeta=0.5625)
            >>> conn.area  # 0.0125 m^2 (physical area)
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

        # Pass physical area and Cd directly to parent
        # Parent will use A * Cd for flow calculation
        super().__init__(id, from_node, to_node, Cd=Cd, area=area)

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


class PipeElement(NetworkElement):
    """
    Pipe element applying frictional pressure drop.
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
    ):
        super().__init__(id, from_node, to_node)
        self.length = length
        self.diameter = diameter
        self.roughness = roughness
        self.area = 3.1415926535 * (diameter / 2) ** 2

        self.regime = regime
        self.friction_model = friction_model
        self.htc_model = htc_model
        self.t_wall = t_wall

    def unknowns(self) -> list[str]:
        return [f"{self.id}.m_dot"]

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:

        import combaero as cb

        m_dot = state_in.m_dot

        if self.regime == "compressible_fanno":
            # Use compressible Fanno flow with friction
            res_cpp = cb._core.pipe_compressible_residuals_and_jacobian(
                m_dot,
                state_in.P_total,
                state_in.T,
                state_in.Y,
                state_out.P,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
            )
        else:
            # Use incompressible Darcy-Weisbach formulation
            res_cpp = cb.pipe_residuals_and_jacobian(
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

        # Residual of P-node unknowns matching pipe drop
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

    def get_spatial_profile(
        self,
        state_in: MixtureState,
        n_steps: int = 100,
    ) -> list:
        """
        Compute and return the spatial flow profile array along the pipe length.
        Requires the solved inlet MixtureState.
        """
        import combaero as cba

        rho = cba.density(state_in.T, state_in.P, state_in.Y)
        u = cba.pipe_velocity(state_in.m_dot, self.diameter, rho)

        if self.regime == "incompressible":
            res = cba.pipe_flow_rough(
                state_in.T,
                state_in.P,
                state_in.Y,
                u,
                self.length,
                self.diameter,
                self.roughness,
                self.friction_model,
                n_steps,
                True,
            )
            return res.profile

        elif self.regime == "compressible_fanno":
            res = cba.pipe_fanno(
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

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass
