"""Canonical diagnostic field schema for network elements and nodes.

This is the single source of truth for what every diagnostics() method must
return after Phase 2.  Use validate_diagnostics() in tests to catch drift.

Field naming conventions
------------------------
- Element pressure / temperature fields use _in / _out suffixes.
- Element thermo/transport properties use a _ref suffix (evaluated at the
  inlet node state) so it is unambiguous which state the value refers to.
- Node stagnation quantities use Tt / Pt (not Tt / Pt).
- Mach at the smaller cross-section of an area change: mach_loss_ref.
- Area ratio of an area change: area_ratio.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Shared field groups
# ---------------------------------------------------------------------------

# Pressure/temperature/Mach block - every element with geometry emits all of these.
ELEMENT_PRESSURE_FIELDS: tuple[str, ...] = (
    "P_in",
    "P_out",
    "T_in",
    "T_out",
    "Pt_in",
    "Pt_out",
    "Tt_in",
    "Tt_out",
    "mach_in",
    "mach_out",
    "dPt",
    "dP",
    "pr_total",
    "pr_static",
)

# Thermodynamic properties at the inlet (reference) state.
ELEMENT_THERMO_REF_FIELDS: tuple[str, ...] = (
    "rho_ref",
    "h_ref",
    "s_ref",
    "u_ref",
    "cp_ref",
    "cv_ref",
    "gamma_ref",
    "a_ref",
    "mw_ref",
)

# Transport properties at the inlet (reference) state.
ELEMENT_TRANSPORT_REF_FIELDS: tuple[str, ...] = (
    "mu_ref",
    "k_ref",
    "nu_ref",
    "Pr_ref",
)

# Convective heat transfer block - present only when an element computes HTC.
HEAT_TRANSFER_FIELDS: tuple[str, ...] = ("Nu", "htc", "T_aw")

# Minimum pressure block for lossless / geometry-free elements.
LOSSLESS_ELEMENT_FIELDS: tuple[str, ...] = (
    "m_dot",
    *ELEMENT_PRESSURE_FIELDS,
    *ELEMENT_THERMO_REF_FIELDS,
    *ELEMENT_TRANSPORT_REF_FIELDS,
)

# Full universal block for elements that have a defined cross-sectional area.
ELEMENT_UNIVERSAL_FIELDS: tuple[str, ...] = (
    "m_dot",
    "velocity",
    "Re",
    *ELEMENT_PRESSURE_FIELDS,
    *ELEMENT_THERMO_REF_FIELDS,
    *ELEMENT_TRANSPORT_REF_FIELDS,
)

# ---------------------------------------------------------------------------
# Node field groups
# ---------------------------------------------------------------------------

NODE_UNIVERSAL_FIELDS: tuple[str, ...] = (
    "T",
    "P",
    "Tt",
    "Pt",
    "mach",
    "rho",
    "h",
    "s",
    "u",
    "gamma",
    "a",
    "cp",
    "cv",
    "mw",
)

NODE_TRANSPORT_FIELDS: tuple[str, ...] = ("mu", "k", "nu", "Pr")

# Extra fields present on chamber nodes (momentum chamber, combustor).
NODE_CHAMBER_FIELDS: tuple[str, ...] = ("Re", "Dh", "velocity", *HEAT_TRANSFER_FIELDS)

# Combustor-only fields.
NODE_COMBUSTOR_FIELDS: tuple[str, ...] = ("phi",)

# ---------------------------------------------------------------------------
# Per-kind diagnostic spec
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DiagnosticSpec:
    """Required and optional field sets for one element or node kind."""

    required: tuple[str, ...]
    optional: tuple[str, ...] = field(default_factory=tuple)


_SPECS: dict[str, DiagnosticSpec] = {
    # -- Elements --
    "OrificeElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "Cd", "mach_throat", "is_correlation"),
    ),
    "EffectiveAreaConnectionElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "Cd", "mach_throat", "is_correlation"),
    ),
    "DiameterDischargeCoefficientConnectionElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "Cd", "mach_throat", "is_correlation"),
    ),
    "PressureLossElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "xi", "theta"),
        optional=(*HEAT_TRANSFER_FIELDS, "f"),
    ),
    "ChannelElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "f", "Dh", *HEAT_TRANSFER_FIELDS),
    ),
    "AreaChangeElement": DiagnosticSpec(
        required=(*ELEMENT_UNIVERSAL_FIELDS, "zeta", "area_ratio", "mach_loss_ref"),
    ),
    "LosslessConnectionElement": DiagnosticSpec(
        required=LOSSLESS_ELEMENT_FIELDS,
    ),
    # -- Nodes --
    "PlenumNode": DiagnosticSpec(
        required=NODE_UNIVERSAL_FIELDS,
    ),
    "MomentumChamberNode": DiagnosticSpec(
        required=(
            *NODE_UNIVERSAL_FIELDS,
            *NODE_TRANSPORT_FIELDS,
            *NODE_CHAMBER_FIELDS,
        ),
    ),
    "CombustorNode": DiagnosticSpec(
        required=(
            *NODE_UNIVERSAL_FIELDS,
            *NODE_TRANSPORT_FIELDS,
            *NODE_CHAMBER_FIELDS,
            *NODE_COMBUSTOR_FIELDS,
        ),
    ),
    "PressureBoundary": DiagnosticSpec(
        required=(*NODE_UNIVERSAL_FIELDS, *NODE_TRANSPORT_FIELDS),
    ),
    "MassFlowBoundary": DiagnosticSpec(
        required=(*NODE_UNIVERSAL_FIELDS, *NODE_TRANSPORT_FIELDS),
    ),
    "EnergyBoundary": DiagnosticSpec(
        required=(*NODE_UNIVERSAL_FIELDS, *NODE_TRANSPORT_FIELDS),
    ),
}


def spec_for(kind: str) -> DiagnosticSpec | None:
    """Return the DiagnosticSpec for the given class name, or None if unknown."""
    return _SPECS.get(kind)


def validate_diagnostics(kind: str, d: dict[str, object]) -> list[str]:
    """Return required fields missing from diagnostic dict *d* for *kind*.

    Returns an empty list when all required fields are present or when *kind*
    is not registered (unknown kinds are not validated).
    """
    spec = _SPECS.get(kind)
    if spec is None:
        return []
    return [f for f in spec.required if f not in d]


# Convenience: all known field names across all specs (for UNIT_MAP coverage checks).
ALL_KNOWN_FIELDS: frozenset[str] = frozenset(
    f for spec in _SPECS.values() for f in (*spec.required, *spec.optional)
)
