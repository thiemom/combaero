from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

# --- Node Data Definitions ---


class PlenumData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


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
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class PressureBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    P_total: float = 101325.0
    T_total: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class CombustorData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    method: Literal["complete", "equilibrium"] = "complete"
    pressure_loss: float | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class MomentumChamberData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    area: float = 0.1
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


# --- Element Data Definitions ---


class PipeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    L: float = 1.0
    D: float = 0.1
    roughness: float = 1e-5
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class OrificeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    area: float = 0.01
    Cd: float = 0.6
    auto_Cd: bool = True
    plate_thickness: float = 0.0  # t [m]: > 0 enables thick-plate correction
    edge_radius: float = 0.0  # r [m]: > 0 enables rounded-entry correction
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class LosslessConnectionData(BaseModel):
    pass


class WallLayerData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    thickness: float = 0.003
    conductivity: float = 20.0
    material: str = "generic"


class ThermalWallData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["thermal"] = "thermal"
    # Support both single-layer (legacy) and multi-layer
    thickness: float | None = None
    conductivity: float | None = None
    layers: list[WallLayerData] | None = None
    area: float | None = 0.05
    R_fouling: float = 0.0


class SolverSettings(BaseModel):
    model_config = ConfigDict(extra="ignore")
    global_regime: Literal["incompressible", "compressible"] = "incompressible"
    init_strategy: Literal["default", "incompressible_warmstart", "homotopy"] = "default"


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
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)


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
    edge_results: dict[str, dict] = {}
