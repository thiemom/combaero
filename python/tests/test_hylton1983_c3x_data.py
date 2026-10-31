"""The C3X transcription (Hylton et al. 1983, NASA CR-168015) checks itself.

validation/cooling/data/hylton1983/ holds tables read off the report's page
images (validation/cooling/extractions/hylton_1983_cr168015_c3x.md). These
tests pin the internal checks the extraction record states, so an edit to a
CSV that breaks one is caught.
"""

from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[2] / "validation" / "cooling" / "data" / "hylton1983"


def _rows(name: str) -> list[dict[str, str]]:
    with open(DATA / name) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_two_unit_columns_agree() -> None:
    """Table III and Fig. 7 print cm and inches; the transcriptions agree."""
    for r in _rows("c3x_coordinates.csv"):
        assert float(r["x_cm"]) / 2.54 == pytest.approx(float(r["x_in"]), abs=1.5e-4), r
        assert float(r["y_cm"]) / 2.54 == pytest.approx(float(r["y_in"]), abs=1.5e-4), r
    for r in _rows("c3x_cooling_holes.csv"):
        for k in ("U", "V", "D"):
            assert float(r[f"{k}_cm"]) / 2.54 == pytest.approx(float(r[f"{k}_in"]), abs=1.5e-3), r


def test_point_29_misprint_is_corrected() -> None:
    p29 = next(r for r in _rows("c3x_coordinates.csv") if r["point"] == "29")
    assert float(p29["y_cm"]) == pytest.approx(0.0411, abs=1e-4)
    assert "0.4115" in p29["note"]


def test_table_ix_pa_column_is_psia_times_133() -> None:
    """The column headed "Pa" is psia x 133.3, not pascals; kPa comes from psia."""
    runs = _rows("c3x_run_conditions.csv")
    assert len(runs) == 18
    for r in runs:
        psia = float(r["PT1_psia"])
        assert float(r["PT1_Pa_printed"]) == pytest.approx(psia * 133.3, rel=2e-3), r
        assert float(r["PT1_kPa_from_psia"]) == pytest.approx(psia * 6.894757, rel=1e-4), r


@pytest.mark.parametrize("run", [149, 158])
def test_coolant_units_agree(run: int) -> None:
    holes = _rows(f"c3x_run{run}_coolant.csv")
    assert [int(h["hole"]) for h in holes] == list(range(1, 11))
    for h in holes:
        assert (float(h["T_avg_F"]) - 32.0) / 1.8 + 273.15 == pytest.approx(
            float(h["T_avg_K"]), abs=0.02
        )
        assert float(h["mdot_lbm_s"]) * 0.45359237 == pytest.approx(float(h["mdot_kg_s"]), rel=6e-3)


@pytest.mark.parametrize("run", [149, 158])
def test_surface_tables_are_complete(run: int) -> None:
    tw = _rows(f"c3x_run{run}_Tw_h.csv")
    assert len(tw) == 75
    for r in tw:
        assert float(r["Tw_over_811K"]) * 811.0 == pytest.approx(float(r["Tw_K"]), abs=0.06)
        assert float(r["h_over_1135"]) * 1135.0 == pytest.approx(float(r["h_W_m2K"]), abs=0.6)
    assert {r["surface"] for r in tw} == {"PS", "SS", "LE"}


def test_our_surface_parametrisation_reproduces_x_over_cax() -> None:
    """All 176 Appendix A rows: our arc length gives the report's x/c_ax."""
    chk = _rows("c3x_station_check.csv")
    assert len(chk) == 176
    assert max(abs(float(r["diff"])) for r in chk) < 0.011


def test_derived_geometry_is_reproducible(tmp_path: Path) -> None:
    """derive_geometry.py rewrites the derived tables byte for byte."""
    work = tmp_path / "hylton1983"
    shutil.copytree(DATA, work)
    subprocess.run(
        [sys.executable, str(work / "derive_geometry.py")], check=True, capture_output=True
    )
    for name in (
        "c3x_cooling_holes_xy.csv",
        "c3x_wall_thickness.csv",
        "c3x_passage_geometry.csv",
        "c3x_station_check.csv",
    ):
        assert (work / name).read_bytes() == (DATA / name).read_bytes(), name


def test_throat_is_within_three_percent_of_table_iv() -> None:
    """The derived throat (3.385 cm) is 2.8% above Table IV's 3.292 cm."""
    header = (DATA / "c3x_passage_geometry.csv").read_text()
    throat = float(
        header.split("Computed throat (min vane-to-neighbour distance) = ")[1].split()[0]
    )
    table_iv = next(
        float(r["value"]) for r in _rows("c3x_cascade_geometry.csv") if r["quantity"] == "throat"
    )
    assert math.isclose(throat / table_iv, 1.028, abs_tol=0.002)
