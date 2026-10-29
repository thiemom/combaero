"""Worked examples from papers (#793): each script and its GUI file agree.

The scripts in python/examples/papers/ run in the examples smoke job
(scripts/run_examples.py). This file pins what the smoke run cannot see:

- the GUI file (gui/examples/<case>.json, served by /examples) solves to the
  same numbers as the script's own network, so the two cannot drift apart;
- the numbers the GUI description panel quotes still come out of the solve;
- every example names its source register entry, and the files ship as
  combaero-gui package data.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tomllib
from pathlib import Path
from types import ModuleType

import pytest

from gui.backend.runner import NetworkRunner

_REPO = Path(__file__).resolve().parents[2]
_GUI_EXAMPLES = _REPO / "gui" / "examples"
_SCRIPTS = _REPO / "python" / "examples" / "papers"
_CASES = sorted(p.stem for p in _GUI_EXAMPLES.glob("*.json"))


def _script(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"paper_example_{name}", _SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _gui(name: str) -> dict:
    data = json.loads((_GUI_EXAMPLES / f"{name}.json").read_text())
    res = NetworkRunner.from_dict(data).solve(timeout=120.0)
    assert res.success, name
    return {k: v.model_dump() for k, v in res._element_results.items()}


def test_the_two_decided_examples_ship() -> None:
    assert _CASES == ["florschuetz1981_jet_array", "iso5167_orifice_metering"]


@pytest.mark.parametrize("name", _CASES)
def test_each_gui_example_has_its_script(name: str) -> None:
    info = json.loads((_GUI_EXAMPLES / f"{name}.json").read_text())["example"]
    assert info["script"] == f"python/examples/papers/{name}.py"
    assert (_REPO / info["script"]).is_file()
    assert info["label"] in ("fidelity", "cross-source")


@pytest.mark.parametrize(
    ("name", "register_key"),
    [
        (
            "florschuetz1981_jet_array",
            "sources/heat_transfer/flohrshuetz/flohrshuetz_truman_metzger_1981.pdf",
        ),
        ("iso5167_orifice_metering", "ISO 5167-2:2003"),
    ],
)
def test_each_script_names_its_source_register_entry(name: str, register_key: str) -> None:
    header = (_SCRIPTS / f"{name}.py").read_text().split('"""')[1]
    assert "validation/sources.yaml" in header
    assert register_key.split("/")[-1] in header.replace("\n    ", "")
    assert register_key in (_REPO / "validation" / "sources.yaml").read_text()
    for section in ("What is reproduced", "What is NOT reproduced", "Validation label"):
        assert section in header, (name, section)


def test_the_examples_ship_as_package_data() -> None:
    cfg = tomllib.loads((_REPO / "gui" / "pyproject.toml").read_text())
    assert "examples/*.json" in cfg["tool"]["setuptools"]["package-data"]["gui"]
    assert "recursive-include examples *.json" in (_REPO / "gui" / "MANIFEST.in").read_text()


def test_florschuetz_gui_file_is_the_scripts_network() -> None:
    """Same rows, same Gc/Gj and Nu, to solver precision; and the panel's
    "Gc/Gj about 0.65 at row 10" and "Eq. 8 to 1%" hold."""
    mod = _script("florschuetz1981_jet_array")
    gui = _gui("florschuetz1981_jet_array")["array"]
    ref = mod.solve_array(5.0, 4.0, 2.0)
    for key in ("rows_Gc_Gj", "rows_Nu", "rows_m_dot"):
        assert gui[key] == pytest.approx(ref[key], rel=1e-6), key
    assert gui["rows_Gc_Gj"][0] == 0.0
    assert gui["rows_Gc_Gj"][-1] == pytest.approx(0.65, abs=0.01)
    for i in range(2, 11):
        eq8 = mod.cb.crossflow_to_jet_ratio_at_row(4.0, 2.0, mod.CD, i)
        assert gui["rows_Gc_Gj"][i - 1] == pytest.approx(eq8, rel=0.01)
    assert min(gui["rows_Re_j"]) > 2500.0 and max(gui["rows_Re_j"]) < 70000.0


def test_florschuetz_fig6_score_at_the_gui_geometry() -> None:
    """The panel quotes 2% mean error at B(5,4)I z/d 2 (measured 2026-10-29:
    9 points, bias -0.9%, MAE 2.1%)."""
    errors = _script("florschuetz1981_jet_array").compare("rc_00", 2.0, verbose=False)
    assert len(errors) == 9
    assert sum(abs(e) for e in errors) / len(errors) < 0.03


def test_iso_gui_file_is_the_scripts_network() -> None:
    """Same flow and differentials; and the panel's "about 6.5 kPa", "about
    4.8 kPa, 73%" hold."""
    mod = _script("iso5167_orifice_metering")
    gui = _gui("iso5167_orifice_metering")["plate"]
    ref = mod.solve_run(0.5, "compressible")
    for key in ("m_dot", "dP_tap", "dP_permanent", "Cd"):
        assert gui[key] == pytest.approx(ref[key], rel=1e-6), key
    assert gui["dP_tap"] == pytest.approx(6.5e3, rel=0.02)
    assert gui["dP_permanent"] == pytest.approx(4.8e3, rel=0.02)
    assert gui["dP_permanent"] / gui["dP_tap"] == pytest.approx(0.73, abs=0.01)


@pytest.mark.parametrize("beta", [0.3, 0.5, 0.7])
def test_iso_incompressible_reads_back_exactly(beta: float) -> None:
    """Run incompressible with eps = 1: the flange differential read through
    ISO's flow equation gives the solved flow back, and the permanent loss is
    5.4.1's long form on the D-D/2 differential, dp_DD2 = dp_fl (C_fl/C_DD2)^2."""
    mod = _script("iso5167_orifice_metering")
    r = mod.solve_run(beta, "incompressible")
    q, C_fl = mod.iso_flow(beta, r["dP_tap"], r["P_approach"], r["T_in"], eps_on=False)
    # 1e-5: the element reads mu for Re_D at the total temperature, the
    # standard at the static one (2e-6 in C at beta 0.7).
    assert q == pytest.approx(r["m_dot"], rel=1e-5)
    dp_dd2 = r["dP_tap"] * (C_fl / r["Cd"]) ** 2
    assert r["dP_permanent"] == pytest.approx(mod.iso_loss_ratio(beta, r["Cd"]) * dp_dd2, rel=1e-6)
