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

    # Verify ref_location
    assert "ref_location" in content, "ref_location missing from quantities.ts"
