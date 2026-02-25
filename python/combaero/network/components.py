from abc import ABC, abstractmethod
from dataclasses import dataclass

import combaero as cb


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


class MomentumChamberNode(NetworkNode):
    """
    Similar to a plenum, but preserves dynamic pressure.
    Enforces P_total = P_static + 0.5 * rho * v^2 instead of P_total = P_static.
    """

    def unknowns(self) -> list[str]:
        return [f"{self.id}.P", f"{self.id}.P_total"]

    def residuals(self, state: MixtureState) -> list[float]:
        # v = m_dot / (rho * A_chamber)
        # For phase 1, we just return a placeholder. The network solver will dynamically
        # inject the momentum conservation equation.
        return [
            state.P_total - state.P
        ]  # Will be dynamically updated by solver with dynamic pressure


class BoundaryNode(NetworkNode):
    """
    Supplies fixed state values. Contributes no unknowns and no residuals.
    """

    def unknowns(self) -> list[str]:
        return []

    def residuals(self, state: MixtureState) -> list[float]:
        return []


class OrificeElement(NetworkElement):
    """
    Incompressible orifice flow. m_dot = Cd * A * sqrt(2 * rho * dP).
    """

    def __init__(self, id: str, from_node: str, to_node: str, Cd: float, area: float):
        super().__init__(id, from_node, to_node)
        self.Cd = Cd
        self.area = area

    def residuals(self, state_in: MixtureState, state_out: MixtureState) -> list[float]:
        state_in.P_total - state_out.P_static  # wait, Orifice depends on topology upstream
        # The phase 1 solver operates iteratively using the tuple interface.
        return []

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
