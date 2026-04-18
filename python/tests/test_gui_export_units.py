"""Unit tests for :mod:`gui.backend.units` used by the CSV export pipeline."""

from __future__ import annotations

from gui.backend.units import UNIT_MAP, label_with_unit


def test_label_with_unit_known_scalar_pressure():
    assert label_with_unit("P") == "P [Pa]"
    assert label_with_unit("P_total") == "P_total [Pa]"
    assert label_with_unit("dP_total") == "dP_total [Pa]"


def test_label_with_unit_known_scalar_temperature_and_flow():
    assert label_with_unit("T") == "T [K]"
    assert label_with_unit("T_aw") == "T_aw [K]"
    assert label_with_unit("m_dot") == "m_dot [kg/s]"


def test_label_with_unit_heat_transfer_fields():
    assert label_with_unit("htc") == "htc [W/m²/K]"
    assert label_with_unit("Q") == "Q [W]"
    assert label_with_unit("Nu") == "Nu [-]"


def test_label_with_unit_meta_columns_pass_through():
    for col in (
        "id",
        "label",
        "type",
        "kind",
        "success",
        "is_boundary",
        "combustion_method",
    ):
        assert label_with_unit(col) == col


def test_label_with_unit_composition_numeric_indices_annotated():
    # Legacy numeric-index form still recognised.
    assert label_with_unit("Y[0]") == "Y[0] [-]"
    assert label_with_unit("Y[12]") == "Y[12] [-]"
    assert label_with_unit("X[3]") == "X[3] [-]"


def test_label_with_unit_composition_species_names_annotated():
    # New species-name form (produced by export).
    assert label_with_unit("Y[N2]") == "Y[N2] [-]"
    assert label_with_unit("Y[H2O]") == "Y[H2O] [-]"
    assert label_with_unit("Y[CH4]") == "Y[CH4] [-]"
    assert label_with_unit("X[O2]") == "X[O2] [-]"
    assert label_with_unit("X[CO2]") == "X[CO2] [-]"


def test_label_with_unit_unknown_column_returns_unchanged():
    assert label_with_unit("foo_bar") == "foo_bar"
    assert label_with_unit("wall_thickness") == "wall_thickness"


def test_unit_map_core_coverage():
    """Guard against accidental removal of critical units."""
    core = {
        "T",
        "P",
        "P_total",
        "m_dot",
        "rho",
        "mach",
        "Re",
        "Nu",
        "htc",
        "Q",
        "xi",
        "p_ratio",
    }
    missing = core - UNIT_MAP.keys()
    assert not missing, f"Missing expected unit-map entries: {sorted(missing)}"
