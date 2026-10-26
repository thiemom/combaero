from .components import (
    CrossflowSegmentElement,
    PressureBoundary,
    MassFlowBoundary,
    CombustorNode,
    EffectiveAreaConnectionElement,
    EffusionPlateElement,
    FilmCoolingElement,
    SlotFilmElement,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    EnergyBoundary,
    LosslessConnectionElement,
    PortDischargeElement,
    NetworkMixtureState,
    MomentumChamberNode,
    OrificeElement,
    ChannelElement,
    PlenumNode,
    PressureLossElement,
    NetworkNode,
    NetworkElement,
    DiameterDischargeCoefficientConnectionElement,
    ConvectiveSurface,
    ImpingementModel,
    PinFinModel,
    RibbedModel,
    SingleJetImpingementModel,
    SmoothModel,
    WallLayer,
    ThermalWall,
    AreaChangeElement,
    TeeJunctionElement,
    VortexElement,
    MultiPortChamberBase,
    BorderCarnotLossElement,
    WallNode,
)
from .mpce_element import ConstantKTeeElement, MultiPortChamberElement
from .ideal_tee_element import IdealMomentumTeeElement
from .combustion import (
    CombustionResult,
    combustion_from_streams,
    combustion_from_phi,
    mix_streams,
)
from .graph import FlowNetwork
from .impingement_array import ImpingementArray
from .effusion_liner import EffusionLiner
from .pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
)
from .schema import DiagnosticSpec, spec_for, validate_diagnostics, ALL_KNOWN_FIELDS
from .solver import NetworkSolver, SolveOutcome

__all__: list[str] = [
    "PlenumNode",
    "MomentumChamberNode",
    "PressureBoundary",
    "MassFlowBoundary",
    "WallNode",
    "CombustorNode",
    "ChannelElement",
    "OrificeElement",
    "EffectiveAreaConnectionElement",
    "EffusionPlateElement",
    "FilmCoolingElement",
    "SlotFilmElement",
    "CrossflowSegmentElement",
    "ImpingementCrossflowElement",
    "ImpingementPlateElement",
    "ImpingementArray",
    "EffusionLiner",
    "LosslessConnectionElement",
    "PortDischargeElement",
    "DiameterDischargeCoefficientConnectionElement",
    "PressureLossElement",
    "AreaChangeElement",
    "TeeJunctionElement",
    "VortexElement",
    "MultiPortChamberBase",
    "MultiPortChamberElement",
    "ConstantKTeeElement",
    "IdealMomentumTeeElement",
    "BorderCarnotLossElement",
    "EnergyBoundary",
    "NetworkMixtureState",
    "CombustionResult",
    "combustion_from_streams",
    "combustion_from_phi",
    "mix_streams",
    "FlowNetwork",
    "NetworkSolver",
    "SolveOutcome",
    "NetworkNode",
    "NetworkElement",
    "DiagnosticSpec",
    "spec_for",
    "validate_diagnostics",
    "ALL_KNOWN_FIELDS",
    "ConvectiveSurface",
    "ImpingementModel",
    "PinFinModel",
    "RibbedModel",
    "SingleJetImpingementModel",
    "SmoothModel",
    "WallLayer",
    "ThermalWall",
    "ConstantFractionLoss",
    "ConstantHeadLoss",
    "LinearThetaFractionLoss",
    "LinearThetaHeadLoss",
]


def __getattr__(name: str) -> type:
    """Deprecated aliases (ConstantLoss, LinearThetaLoss) warn on access."""
    from . import pressure_loss as _pressure_loss

    if name in _pressure_loss._DEPRECATED_ALIASES:
        return getattr(_pressure_loss, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
