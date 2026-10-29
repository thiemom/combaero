import os
import re
import sys

# Add GUI backend to sys.path so we can import units.py
gui_backend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../gui/backend"))
sys.path.insert(0, gui_backend_dir)

import units  # noqa: E402


def test_gui_backend_units_catalogue():
    """Verify the backend unit mapping matches Phase 2/4 roadmap."""
    assert "Tt" in units.UNIT_MAP
    assert "Pt" in units.UNIT_MAP
    assert units.UNIT_MAP["Tt"] == "K"
    assert units.UNIT_MAP["Pt"] == "Pa"

    assert "mach_loss_ref" in units.UNIT_MAP
    assert "area_ratio" in units.UNIT_MAP
    assert "dP" in units.UNIT_MAP
    assert "pr_static" in units.UNIT_MAP

    assert "ref_location" in units._META_COLUMNS


def test_label_with_unit():
    """Verify column header formatting."""
    assert units.label_with_unit("Tt") == "Tt [K]"
    assert units.label_with_unit("Pt") == "Pt [Pa]"
    assert units.label_with_unit("dP") == "dP [Pa]"
    assert units.label_with_unit("ref_location") == "ref_location"
    assert units.label_with_unit("Y[N2]") == "Y[N2] [-]"
    # State::mw is g/mol (include/units_data.h); the export said kg/mol.
    assert units.label_with_unit("mw") == "mw [g/mol]"
    # Every key with a unit is labelled in the export, the junction port and
    # step patterns and the arrays' per-row lists included (#774, 0.8.0).
    assert units.unit_of("port_2_P") == "Pa"
    assert units.unit_of("step_branch_dPt") == "Pa"
    assert units.label_with_unit("port_2_P") == "port_2_P [Pa]"
    assert units.label_with_unit("step_branch_dPt") == "step_branch_dPt [Pa]"
    assert units.label_with_unit("port_0_area") == "port_0_area [m^2]"
    assert units.label_with_unit("entrainment_ratio") == "entrainment_ratio [-]"
    assert units.label_with_unit("rows_T_hot") == "rows_T_hot [K]"
    assert units.label_with_unit("UA") == "UA [W/K]"
    # Text-valued keys carry no unit.
    assert units.is_text_key("port_1_area_source")
    assert units.label_with_unit("port_1_area_source") == "port_1_area_source"
    assert units.label_with_unit("second_law") == "second_law"


def test_frontend_quantities_sync():
    """
    Static analysis test to verify frontend quantities.ts matches backend keys.
    """
    frontend_ts_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../gui/frontend/src/utils/quantities.ts")
    )
    if not os.path.exists(frontend_ts_path):
        return  # Skip if frontend is not present

    with open(frontend_ts_path, encoding="utf-8") as f:
        content = f.read()

    # Verify Tt and Pt are present
    assert re.search(r'Tt:\s*\{\s*label:\s*"Tt"', content), "Tt not found in quantities.ts"
    assert re.search(r'Pt:\s*\{\s*label:\s*"Pt"', content), "Pt not found in quantities.ts"

    # Verify Phase 2 specific area change keys
    assert "mach_loss_ref" in content, "mach_loss_ref missing from quantities.ts"
    assert "area_ratio" in content, "area_ratio missing from quantities.ts"
    assert "dP" in content, "dP missing from quantities.ts"
    assert "pr_static" in content, "pr_static missing from quantities.ts"

    # ref_location is a string (a meta column of the CSV export, see
    # test_gui_backend_units_catalogue); the catalogue holds numeric keys only.
    # The full key-by-key sync is test_gui_result_keys_sync.py plus
    # gui/frontend/src/utils/quantities.test.ts.
    assert "ref_location:" not in content
