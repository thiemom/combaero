"""Unit tests for :mod:`gui.backend.units` used by the CSV export pipeline."""

import csv
import io
import sys
import warnings
from pathlib import Path

import pytest

from gui.backend.units import _META_COLUMNS, UNIT_MAP, is_text_key, label_with_unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import gui_result_keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402

# One network per palette kind plus the variants that switch on optional
# diagnostics (scripts/gui_result_keys.py, the source of the GUI's key list).
_NETWORKS = list(gui_result_keys.networks())


@pytest.mark.parametrize(
    "graph",
    [g for _, _, g in _NETWORKS],
    ids=[f"{kind}-{variant}" for kind, variant, _ in _NETWORKS],
)
def test_every_export_column_has_a_unit(graph: dict) -> None:
    """Every column of the CSV export carries its unit (#774), unless it is a
    meta column (id, label, status ...) or a text-valued key (is_text_key).
    A new result key with no unit fails here: add it to UNIT_MAP in
    gui/backend/units.py (or to TEXT_KEYS when its value is text)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = http._post("/export", graph)
    assert r.status_code == 200, r.text[:300]
    header = next(csv.reader(io.StringIO(r.text)))
    unlabelled = [
        c
        for c in header
        if not (c in _META_COLUMNS or is_text_key(c) or (" [" in c and c.endswith("]")))
    ]
    assert unlabelled == []


def test_label_with_unit_known_scalar_pressure():
    assert label_with_unit("P") == "P [Pa]"
    assert label_with_unit("Pt") == "Pt [Pa]"
    assert label_with_unit("dPt") == "dPt [Pa]"


def test_label_with_unit_known_scalar_temperature_and_flow():
    assert label_with_unit("T") == "T [K]"
    assert label_with_unit("T_aw") == "T_aw [K]"
    assert label_with_unit("m_dot") == "m_dot [kg/s]"


def test_label_with_unit_heat_transfer_fields():
    assert label_with_unit("htc") == "htc [W/m^2/K]"
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
        "Pt",
        "m_dot",
        "rho",
        "mach",
        "Re",
        "Nu",
        "htc",
        "Q",
        "xi",
        "pr_total",
    }
    missing = core - UNIT_MAP.keys()
    assert not missing, f"Missing expected unit-map entries: {sorted(missing)}"


def test_unit_map_covers_all_schema_fields():
    """Every field declared in schema.py must have a unit in UNIT_MAP.

    This test catches UNIT_MAP drift whenever schema.py gains new fields.
    Composition keys (Y[...], X[...]) are handled by the regex branch in
    label_with_unit and are intentionally excluded from UNIT_MAP.
    """
    from combaero.network.schema import ALL_KNOWN_FIELDS

    missing = ALL_KNOWN_FIELDS - UNIT_MAP.keys()
    assert not missing, (
        f"Schema fields not covered by UNIT_MAP: {sorted(missing)}\n"
        "Add them to gui/backend/units.py."
    )
