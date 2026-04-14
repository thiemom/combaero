from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SmoothModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["smooth"] = "smooth"


class RibbedModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["ribbed"] = "ribbed"
    e_D: float = 0.05
    pitch_to_height: float = 10.0
    alpha_deg: float = 90.0


class DimpledModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["dimpled"] = "dimpled"
    d_Dh: float = 0.2
    h_d: float = 0.15
    S_d: float = 2.0


class PinFinModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["pin_fin"] = "pin_fin"
    pin_diameter: float = 0.005
    channel_height: float = 0.01  # [m]
    S_D: float = 2.5
    X_D: float = 2.5
    N_rows: int = 10
    is_staggered: bool = True


class ImpingementModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["impingement"] = "impingement"
    d_jet: float = 0.002
    z_D: float = 4.0
    x_D: float = 6.0
    y_D: float = 6.0
    A_target: float = 0.01
    Cd_jet: float = 0.8


SurfaceModelData = (
    SmoothModelData | RibbedModelData | DimpledModelData | PinFinModelData | ImpingementModelData
)


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
    Dh: float = 0.1
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


# --- Element Data Definitions ---


class ChannelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    L: float = 1.0
    D: float = 0.1
    roughness: float = 1e-5
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class OrificeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    diameter: float = 0.08
    Cd: float = 0.6
    auto_Cd: bool = True
    plate_thickness: float = 0.0  # t [m]: > 0 enables thick-plate correction
    edge_radius: float = 0.0  # r [m]: > 0 enables rounded-entry correction
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None

    @model_validator(mode="before")
    @classmethod
    def migrate_area_to_diameter(cls, data: dict) -> dict:
        if isinstance(data, dict) and "area" in data and "diameter" not in data:
            import math

            area_val = data.pop("area")
            if area_val > 0:
                data["diameter"] = math.sqrt(4.0 * area_val / math.pi)
        return data


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

    @model_validator(mode="before")
    @classmethod
    def migrate_legacy_fields(cls, data: dict) -> dict:
        if isinstance(data, dict) and ("layers" not in data or data["layers"] is None):
            # If we have legacy fields but no layers, migrate them
            t = data.get("thickness", 0.003)
            k = data.get("conductivity", 20.0)
            # Ensure they are not None
            t = t if t is not None else 0.003
            k = k if k is not None else 20.0
            data["layers"] = [{"thickness": t, "conductivity": k}]
        return data


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
    type: str | None = None
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
    final_norm: float | None = None
    node_results: dict[str, NodeResult] = {}
    element_results: dict[str, ElementResult] = {}
    edge_results: dict[str, dict] = {}
