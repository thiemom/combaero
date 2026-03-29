from pydantic import BaseModel, Field, ConfigDict

# --- Node Data Definitions ---


class PlenumData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    pass  # Standard plenum just takes connections


class MassBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    m_dot: float = 1.0
    T_total: float = 300.0
    Y: list[float] | None = None


class PressureBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    P_total: float = 101325.0


class MomentumChamberData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    pass


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


class NodeResult(BaseModel):
    T: float
    P: float
    Y: list[float]
    success: bool = True


class ElementResult(BaseModel):
    m_dot: float
    success: bool = True


class SolveResponse(BaseModel):
    success: bool
    message: str = ""
    node_results: dict[str, NodeResult] = {}
    element_results: dict[str, ElementResult] = {}
