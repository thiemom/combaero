"""Unit catalogue for result-export column headers.

Mirrors the labels/units of the frontend ``QUANTITY_CATALOGUE`` so the CSV
produced by ``/export`` is self-describing (``P [Pa]``, ``T [K]`` …).

Keep this module dependency-free; it is imported by ``main.py`` during the
CSV-streaming response.
"""

from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Scalar quantity → unit mapping
# ---------------------------------------------------------------------------

UNIT_MAP: dict[str, str] = {
    # --- Thermodynamic state ---
    "T": "K",
    "T_total": "K",
    "T_in": "K",
    "T_out": "K",
    "Tt_in": "K",
    "Tt_out": "K",
    "T_aw": "K",
    "T_hot": "K",
    "T_interface": "K",
    "P": "Pa",
    "P_total": "Pa",
    "P_in": "Pa",
    "P_out": "Pa",
    "Pt_in": "Pa",
    "Pt_out": "Pa",
    "dP_total": "Pa",
    "dP": "Pa",
    # --- Flow ---
    "m_dot": "kg/s",
    "velocity": "m/s",
    "v": "m/s",
    "a": "m/s",
    "mach": "-",
    "mach_in": "-",
    "mach_out": "-",
    "mach_throat": "-",
    # --- Transport ---
    "mu": "Pa·s",
    "k": "W/m/K",
    "nu": "m²/s",
    "Pr": "-",
    "Re": "-",
    # --- Thermo per-mass ---
    "rho": "kg/m³",
    "h": "J/kg",
    "s": "J/kg/K",
    "u": "J/kg",
    "cp": "J/kg/K",
    "cv": "J/kg/K",
    "gamma": "-",
    "mw": "kg/mol",
    # --- Heat transfer ---
    "Nu": "-",
    "htc": "W/m²/K",
    "h_a": "W/m²/K",
    "h_b": "W/m²/K",
    "Q": "W",
    # --- Element diagnostics ---
    "xi": "-",
    "theta": "-",
    "p_ratio": "-",
    "p_ratio_total": "-",
    "f": "-",
    "Cd": "-",
    "phi": "-",
    # --- Geometry ---
    "area": "m²",
    "diameter": "m",
    "length": "m",
    "roughness": "m",
}

# Composition keys render as ``Y[<label>] [-]``. ``<label>`` may be a species
# name (e.g. ``N2``, ``H2O``, ``CH4``) or a numeric index for forward
# compatibility with older exports.
_COMPOSITION_RE = re.compile(r"^[YX]\[[A-Za-z0-9_+\-()]+\]$")


def label_with_unit(column: str) -> str:
    """Return ``'<column> [<unit>]'`` when a unit is known, else ``column``.

    Meta columns such as ``id``, ``type``, ``label``, ``kind``, ``success``
    and ``is_boundary`` are passed through unchanged.
    """
    if column in _META_COLUMNS:
        return column
    if column in UNIT_MAP:
        return f"{column} [{UNIT_MAP[column]}]"
    if _COMPOSITION_RE.match(column):
        return f"{column} [-]"
    return column


_META_COLUMNS: frozenset[str] = frozenset(
    {
        "id",
        "label",
        "type",
        "kind",
        "success",
        "is_boundary",
        "combustion_method",
    }
)
