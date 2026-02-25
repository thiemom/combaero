from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import combaero as cb

if TYPE_CHECKING:
    from .graph import FlowNetwork

CombustionMethodLiteral = Literal["complete", "equilibrium"]


@dataclass
class MixtureState:
    P: float  # static pressure [Pa]
    P_total: float  # total pressure [Pa]  (= P for plenum nodes)
    T: float  # static temperature [K]
    T_total: float  # total temperature [K]  (= T for plenum nodes)
    m_dot: float  # mass flow rate [kg/s]
    X: list[float]  # mole fractions [14 species]

    def density(self) -> float:
        return cb.density(self.T, self.P, self.X)

    def enthalpy(self) -> float:
        return cb.h_mass(self.T, self.X)

    def speed_of_sound(self) -> float:
        return cb.speed_of_sound(self.T, self.X)


class NetworkNode(ABC):
    def __init__(self, id: str):
        self.id = id

    @abstractmethod
    def unknowns(self) -> list[str]:
        """Names of unknowns this node contributes to the solver."""
        pass

    @abstractmethod
    def residuals(self, state: MixtureState) -> list[float]:
        """Returns zero when node equations are satisfied."""
        pass

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
    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:
        """Returns zero when element equations are satisfied."""
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
    Enforces energy balance (not currently active for Phase 1 pressure networks).
    """

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P"]

    def residuals(self, state: MixtureState) -> list[float]:
        # Residual equations will be governed dynamically by the solver network sum.
        # But locally, P_total must equal P_static.
        return [state.P_total - state.P]

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class MomentumChamberNode(NetworkNode):
    """
    Similar to a plenum, but preserves dynamic pressure.
    Enforces P_total = P_static + 0.5 * rho * v^2 instead of P_total = P_static.
    Supports directional flow vectors (port angles) for momentum conservation.
    """

    def __init__(self, id: str, port_angles_deg: dict[str, float] | None = None):
        super().__init__(id)
        # Dictionary mapping connected element IDs to angle relative to chamber axis [deg]
        self.port_angles_deg: dict[str, float] = port_angles_deg or {}

    def set_port_angle(self, element_id: str, angle_deg: float) -> None:
        """Set the angle of the flow for a specific connected element."""
        self.port_angles_deg[element_id] = angle_deg

    def get_port_angle(self, element_id: str) -> float:
        """Get the angle of the flow for a specific connected element (defaults to 0.0)."""
        return self.port_angles_deg.get(element_id, 0.0)

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def residuals(self, state: MixtureState) -> list[float]:
        # v = m_dot / (rho * A_chamber)
        # For phase 1, we just return a placeholder. The network solver will dynamically
        # inject the momentum conservation equation.
        return [
            state.P_total - state.P
        ]  # Will be dynamically updated by solver with dynamic pressure

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class BoundaryNode(NetworkNode):
    """
    Supplies fixed state values. Contributes no unknowns and no residuals.
    """

    def unknowns(self) -> list[str]:
        return []

    def residuals(self, state: MixtureState) -> list[float]:
        return []

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class CombustorNode(NetworkNode):
    """
    Combustion chamber adding energy (and optionally mass) to the flow.
    Evaluates adiabatic combustion using the C++ core backend based on
    the selected CombustionMethodLiteral.
    """

    def __init__(
        self,
        id: str,
        method: CombustionMethodLiteral = "complete",
        pressure_loss_frac: float = 0.04,
    ):
        super().__init__(id)
        self.method = method
        self.pressure_loss_frac = pressure_loss_frac

    def unknowns(self) -> list[str]:
        # P, T, and species mass fractions depending on the solver architecture
        return [f"{self.id}.P", f"{self.id}.T"]

    def residuals(self, state: MixtureState) -> list[float]:
        # Energy and Mass constraints evaluated dynamically by the FlowNetwork loop
        return [0.0]

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass


class OrificeElement(NetworkElement):
    """
    Incompressible orifice flow. m_dot = Cd * A * sqrt(2 * rho * dP).
    """

    def __init__(self, id: str, from_node: str, to_node: str, Cd: float, area: float):
        super().__init__(id, from_node, to_node)
        self.Cd = Cd
        self.area = area
        self.upstream_diameter: float | None = None

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        # Ask graph what element fed into our "from_node"
        upstream_elements = graph.get_upstream_elements(self.from_node)

        if len(upstream_elements) == 1:
            upstream_element = upstream_elements[0]
            # If the upstream element is a pipe, extract its diameter
            if isinstance(upstream_element, PipeElement):
                self.upstream_diameter = upstream_element.diameter

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:
        state_in.P_total - state_out.P  # wait, Orifice depends on topology upstream
        # The phase 1 solver operates iteratively using the tuple interface.
        return []

    def n_equations(self) -> int:
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

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:
        # Evaluates the native C++ tuple (f, J) interface for lossless P_total conservation
        import combaero as cb

        return [cb.lossless_pressure_and_jacobian(state_in.P_total, state_out.P_total)[0]]

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
    ):
        super().__init__(id, from_node, to_node)
        self.length = length
        self.diameter = diameter
        self.roughness = roughness
        self.area = 3.1415926535 * (diameter / 2) ** 2

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:
        # Phase 1 simple placeholder for now.
        return []

    def n_equations(self) -> int:
        return 1

    def resolve_topology(self, graph: "FlowNetwork") -> None:
        pass
