from typing import Literal
from pydantic import BaseModel, Field, ConfigDict

# --- Node Data Definitions ---


class PlenumData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    pass  # Standard plenum just takes connections


class CompositionData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    mode: Literal["mole", "mass"] = "mole"
    source: Literal["dry_air", "humid_air", "fuel", "custom"] = "dry_air"
    custom_fractions: dict[str, float] | None = None
    relative_humidity: float = 0.6  # [0-1] for humid_air
    ambient_T: float = 288.15
    ambient_P: float = 101325.0


class MassBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    m_dot: float = 1.0
    T_total: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)


class PressureBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    P_total: float = 101325.0
    T_total: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)


class CombustorData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    method: Literal["complete", "equilibrium"] = "complete"
    pressure_loss: float | None = None


class MomentumChamberData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    area: float = 0.1


# --- Element Data Definitions ---


class PipeData(BaseModel):
    L: float = 1.0
    D: float = 0.1
    roughness: float = 1e-5


class OrificeData(BaseModel):
    area: float = 0.01
    Cd: float = 0.6


class LosslessConnectionData(BaseModel):
    pass


# --- React Flow Wrapper Schemas ---


class NodePosition(BaseModel):
    model_config = ConfigDict(extra="ignore")
    x: float
    y: float


class ReactFlowNode(BaseModel):
    id: str
    type: str
    position: NodePosition
    data: dict


class ReactFlowEdge(BaseModel):
    id: str
    source: str
    target: str
    data: dict | None = None


class NetworkGraphSchema(BaseModel):
    nodes: list[ReactFlowNode]
    edges: list[ReactFlowEdge]


# --- Results Schemas ---


class StateResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    T: float
    P: float
    P_total: float | None = None
    rho: float | None = None
    h: float | None = None
    s: float | None = None
    mach: float | None = None
    Y: list[float]
    X: list[float] | None = None


class NodeResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    state: StateResult
    success: bool = True


class ElementResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    m_dot: float
    success: bool = True


class SolveResponse(BaseModel):
    success: bool
    message: str = ""
    node_results: dict[str, NodeResult] = {}
    element_results: dict[str, ElementResult] = {}
