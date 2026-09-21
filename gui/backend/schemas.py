from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SmoothModelData(BaseModel):
    model_config = ConfigDict(extra="ignore")
    type: Literal["smooth"] = "smooth"


class RibbedModelData(BaseModel):
    """Rib-roughened walls, on a provenanced correlation set.

    Restored in 0.8.0 after the previous rib correlation was removed for
    failing provenance review -- its friction multiplier was 4-5x below the
    only rib datum in the repository. See issue #334.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["ribbed"] = "ribbed"
    e_D: float = 0.06
    p_e: float = 10.0
    alpha_deg: float = 90.0
    # Aspect ratio. Carried here rather than on the channel because a
    # hydraulic diameter does not determine one -- a 2:1 duct and a square
    # duct can share a Dh.
    W_H: float = 1.0
    # 1, 2 or 4. Two means two OPPOSITE walls, the configuration Han measured.
    n_ribbed_walls: int = 2
    # A user knob on the smooth walls' contribution. Plain smooth walls give
    # h_s/h_r ~ 0.42 where Han's own channel average implies 0.70, because ribs
    # enhance the adjacent smooth wall too. About 1.67 reproduces Han.
    smooth_wall_Nu_multiplier: float = 1.0


class ImpingementModelData(BaseModel):
    """Jet array impingement cooling, one row (Florschuetz, Truman and
    Metzger 1981), on a real crossflow term.

    Restored in 0.9.0 after the previous impingement correlation was removed
    for citing this same paper while taking no crossflow input at all,
    making it structurally unable to be that correlation. See issue #337.

    One element models one spanwise row. A full array is a chain of these,
    each with its own ``row`` -- there is no single "channel Nu" for a whole
    array the way there is for a smooth or ribbed duct, since downstream
    rows see progressively more crossflow than upstream ones.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["impingement"] = "impingement"
    d_jet: float = 0.002
    xn_d: float = 8.0
    yn_d: float = 6.0
    z_d: float = 2.0
    # 1-indexed, counting from upstream. Row 1 sees zero crossflow by
    # definition (Gc/Gj = 0, Florschuetz's own Nu1).
    row: int = 1
    # Jet-plate discharge coefficient. 0.79 is the source's own recommended
    # default absent a measured value; combaero has no discharge-coefficient
    # correlation of its own for a jet-plate array yet (issue #375).
    C_D: float = 0.79


class SingleJetImpingementModelData(BaseModel):
    """A single free round jet (Goldstein, Behbahani and Heppelmann 1986),
    at a representative radial position.

    Unlike ImpingementModelData, there is no array and no crossflow. Nu is
    LOCAL to the radial position R_D, not area-averaged -- this reports one
    representative value for the target patch, not a profile.
    """

    model_config = ConfigDict(extra="ignore")
    type: Literal["single_jet_impingement"] = "single_jet_impingement"
    bc: Literal["constant_heat_flux", "constant_wall_temperature"] = "constant_heat_flux"
    d_jet: float = 0.003
    # Jet-to-target-plate spacing / d_jet. 7.75 is the correlation's own
    # optimum spacing.
    L_D: float = 7.75
    # Radial distance from the jet centerline / d_jet. 5.0 matches the
    # source's own closed-form check point, a worked example, not a default
    # that suits every rig.
    R_D: float = 5.0


# Dimpled and pin-fin surfaces were removed in 0.7.0 -- their correlations
# could not be traced to their cited sources -- and remain deferred (issue
# #339). Saved networks carrying them are rejected with a message naming the
# reason; see graph_builder. Ribbed and impingement were both restored on
# provenanced replacements (#334, #337).
SurfaceModelData = (
    SmoothModelData | RibbedModelData | ImpingementModelData | SingleJetImpingementModelData
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


class MPCETeeData(BaseModel):
    """Mynard-based multi-port chamber tee element with constrained topology.

    3-port junction using Mynard's Unified0D pseudosupplier residual. The
    ``flow_direction`` field declares whether this junction is a merge
    (joining flow, 2 inlets + 1 outlet) or a branch (separating flow,
    1 inlet + 2 outlets); the runtime element raises if the observed
    flow pattern disagrees.
    """

    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    flow_direction: Literal["merge", "branch"] = "branch"
    # Junction closure model. "mynard" (default) is the full Unified0D
    # pseudosupplier physics; "constant_k" is the simplest-model tier:
    # fixed handbook loss coefficients referenced to the common-leg
    # dynamic head (Idelchik convention), K independent of the split.
    junction_model: Literal["mynard", "constant_k"] = "mynard"
    K_straight: float = 0.4  # constant_k only: straight-arm loss coeff [-]
    K_branch: float = 1.0  # constant_k only: branch-arm loss coeff [-]
    theta_deg: float = 90.0  # branch angle [deg], converted to radians in graph_builder
    F_C: float | None = None  # None = inherit from common/straight-arm channel
    F_branch: float | None = None  # None = inherit from branch-arm channel
    psi: float = 1.0  # fallback ratio used only when F_branch is None and not inherited
    # Joining-side etransfer correction scale (combaero extension to faithful
    # Mynard). Default 0.2 calibrated against Bassett K11/K12 analytical +
    # Idelchik 1966 tabulated values. Set to 0.0 to disable the correction
    # (faithful Mynard); set to a custom value for re-calibrated workflows.
    joining_etransfer_alpha: float | None = None
    initial_guess: dict[str, float] = Field(default_factory=dict)


class EjectorData(BaseModel):
    """Critical-mode supersonic ejector (3-port: primary + secondary inlets,
    one outlet). Wraps ``combaero.network.ejector_element.EjectorElement``.

    Geometry is three absolute areas (the ratios are computed by the
    element). Validation (throat > 0, nozzle_exit > throat, mixing >
    nozzle_exit, recovery_efficiency > 0) is enforced by the element
    constructor at build time; the frontend mirrors it for early feedback.
    """

    model_config = ConfigDict(extra="ignore")
    label: str | None = None
    throat_area: float = 3.14e-5  # primary nozzle throat area A_t [m^2]
    nozzle_exit_area: float = 1.0e-4  # primary nozzle exit area A_p1 [m^2], > throat
    mixing_area: float = 8.0e-4  # constant-area mixing section A_3 [m^2], > nozzle_exit
    # Multiplies the lossless mixed stagnation pressure to give the critical
    # back pressure P_c*. 1.0 = no artificial loss (not fitted to data).
    recovery_efficiency: float = 1.0
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
    init_strategy: Literal[
        "default",
        "analytical_pt_prop",
        "incompressible_warmstart",
        "outletref_warmstart",
        "homotopy",
        "continuation",
    ] = "default"
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
