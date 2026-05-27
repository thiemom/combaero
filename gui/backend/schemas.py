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
    source: Literal["dry_air", "humid_air", "fuel", "custom"] = "humid_air"
    custom_fractions: dict[str, float] | None = None
    relative_humidity: float | None = None  # None = use global or 0.6
    ambient_T: float | None = None  # None = use global or 288.15
    ambient_P: float | None = None  # None = use global or 101325.0


class MassBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    m_dot: float = 1.0
    Tt: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class PressureBoundaryData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    Pt: float = 101325.0
    Tt: float = 300.0
    composition: CompositionData = Field(default_factory=CompositionData)
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class WallData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)


class ConstantFractionLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["constant_fraction"] = "constant_fraction"
    xi: float = 0.03


class LinearThetaFractionLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["linear_theta_fraction"] = "linear_theta_fraction"
    k: float = 0.5
    xi0: float = 0.02


class ConstantHeadLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["constant_head"] = "constant_head"
    zeta: float = 5.0


class LinearThetaHeadLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["linear_theta_head"] = "linear_theta_head"
    k: float = 1.0
    zeta0: float = 3.0


PressureLossData = (
    ConstantFractionLossData
    | LinearThetaFractionLossData
    | ConstantHeadLossData
    | LinearThetaHeadLossData
)


class CombustorData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    method: Literal["complete", "equilibrium"] = "complete"
    area: float | None = None  # None = derive from Dh (pi/4 * Dh^2)
    Dh: float | None = None
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class MomentumChamberData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    area: float | None = None  # None = derive from Dh (pi/4 * Dh^2)
    Dh: float | None = None  # None = inherit from upstream channel (Dh = D)
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


# --- Element Data Definitions ---


class ChannelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    L: float = 1.0
    D: float | None = None  # None = inherit from upstream element geometry
    Dh: float | None = None  # None = circular (Dh = D); override for non-circular ducts
    roughness: float = 1e-5
    friction_model: Literal["haaland", "serghides", "colebrook", "petukhov"] = "haaland"
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    regime: Literal["default", "incompressible", "compressible"] = "default"
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class OrificeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    diameter: float | None = None  # None = inherit from upstream channel
    Cd: float = 0.6
    correlation: Literal[
        "ReaderHarrisGallagher", "Stolz", "Miller", "ThickPlate", "RoundedEntry", "fixed"
    ] = "ReaderHarrisGallagher"
    plate_thickness: float = 0.0  # t [m]
    edge_radius: float = 0.0  # r [m]
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


class AreaChangeData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    model_type: Literal["sharp", "conical"] = "sharp"
    F0: float | None = None  # None = inherit from upstream channel (A = pi/4*D^2)
    F1: float | None = None  # None = inherit from downstream channel
    length: float | None = None
    D_h: float = 0.0
    initial_guess: dict[str, float] = Field(default_factory=dict)
    label: str | None = None


class LosslessConnectionData(BaseModel):
    pass


class TeeJunctionData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    tee_type: str = "merging"  # "merging" | "branching"
    theta_deg: float = 90.0  # branch angle [deg], converted to radians in graph_builder
    F_C: float | None = None  # None = inherit from common/straight-arm channel (A = pi/4*D^2)
    F_branch: float | None = (
        None  # None = inherit from branch-arm channel; psi is computed from F_C/F_branch
    )
    psi: float = 1.0  # fallback ratio used only when F_branch is None and not inherited
    initial_guess: dict[str, float] = Field(default_factory=dict)


class VortexData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    r_c: float = 0.02  # vortex core radius [m]
    r_out: float = 0.10  # outer radius where pressure is evaluated [m]
    r_in: float = 0.0  # inner radius (0 = on-axis)
    omega_rpm: float | None = None  # shaft speed [rpm]; None = use global solver setting
    n: float = 2.0  # Vatistas shape parameter (>= 1, default n=2)
    initial_guess: dict[str, float] = Field(default_factory=dict)


class DiscreteLossData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    correlation_type: Literal[
        "constant_fraction", "constant_head", "linear_theta_fraction", "linear_theta_head"
    ] = "constant_fraction"
    xi: float = 0.03
    k: float = 0.001
    xi0: float = 0.02
    zeta: float = 1.0
    zeta0: float = 1.0
    area: float | None = None
    theta_source: str | None = None
    surface: SurfaceModelData = Field(default_factory=SmoothModelData)
    Nu_multiplier: float = 1.0
    f_multiplier: float = 1.0
    label: str | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)


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
    init_strategy: Literal["default", "incompressible_warmstart", "homotopy", "continuation"] = (
        "default"
    )
    method: Literal[
        "hybr",
        "lm",
        "broyden1",
        "broyden2",
        "anderson",
        "linearmixing",
        "diagbroyden",
        "excitingmixing",
        "krylov",
        "df-sane",
    ] = "hybr"
    timeout: float | None = 180.0
    omega_rpm: float | None = None  # global shaft speed [rpm]; None = disabled
    Nu_multiplier: float | None = (
        None  # global heat-transfer scale; stacks multiplicatively with per-element
    )
    f_multiplier: float | None = (
        None  # global friction scale; stacks multiplicatively with per-element
    )
    ambient_T: float | None = None  # ISO 2314: 288.15 K; overrides per-node value when humid_air
    ambient_P: float | None = None  # ISO 2314: 101325 Pa; overrides per-node value when humid_air
    ambient_RH: float | None = (
        None  # ISO 2314: 0.6; overrides per-node relative_humidity when humid_air
    )


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
    sourceHandle: str | None = None
    targetHandle: str | None = None
    type: str | None = None
    data: dict | None = None


class NetworkGraphSchema(BaseModel):
    nodes: list[ReactFlowNode]
    edges: list[ReactFlowEdge]
    solver_settings: SolverSettings = Field(default_factory=SolverSettings)

    @model_validator(mode="before")
    @classmethod
    def _normalise_solver_settings_key(cls, data: object) -> object:
        if isinstance(data, dict) and "solverSettings" in data and "solver_settings" not in data:
            data = dict(data)
            data["solver_settings"] = data.pop("solverSettings")
        return data


# --- Results Schemas ---


class StateResult(BaseModel):
    model_config = ConfigDict(extra="allow")
    T: float
    P: float
    Pt: float | None = None
    Tt: float | None = None
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
    m_dot: float = 0.0
    success: bool = True


class ConvergencePoint(BaseModel):
    eval: int
    t: float
    norm: float


class WorstResidual(BaseModel):
    name: str
    residual: float


class SolveResponse(BaseModel):
    success: bool
    message: str = ""
    final_norm: float | None = None
    node_results: dict[str, NodeResult] = {}
    element_results: dict[str, ElementResult] = {}
    edge_results: dict[str, dict] = {}
    convergence_history: list[ConvergencePoint] = []
    worst_residuals: list[WorstResidual] = []
    solver_settings_used: dict = {}
    lm_started_at_eval: int | None = None
