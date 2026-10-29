"""The results table of gui.backend.results (#784): one row per case,
(kind, entity, quantity) columns, units and solver state in attrs."""

import copy
import csv
import inspect
import io
import math
import sys
import typing
import warnings
from pathlib import Path

import pandas as pd
import pytest

from gui.backend.main import _neutralise_formulas
from gui.backend.results import (
    NETWORK_ENTITY,
    SOLVE_ENTITY,
    append_case,
    load_results,
    varied_inputs,
)
from gui.backend.runner import NetworkRunner
from gui.backend.units import input_unit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import gui_result_keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402

import gui.backend.schemas as schemas  # noqa: E402

_SIMPLE: dict = {
    "nodes": [
        {
            "id": "inlet",
            "type": "mass_boundary",
            "position": {"x": 0, "y": 0},
            "data": {
                "label": "inlet",
                "m_dot": 1.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
        {
            "id": "n_chan",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {
                "label": "chan",
                "L": 1.0,
                "D": 0.1,
                "roughness": 1e-5,
                "initial_guess": {"m_dot": 0.9},
            },
        },
        {
            "id": "outlet",
            "type": "pressure_boundary",
            "position": {"x": 400, "y": 0},
            "data": {
                "label": "outlet",
                "Pt": 101325.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
                "initial_guess": {"Pt": 101325.0},
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "inlet", "target": "n_chan"},
        {"id": "e2", "source": "n_chan", "target": "outlet"},
    ],
}


@pytest.fixture(scope="module")
def runner() -> NetworkRunner:
    return NetworkRunner.from_dict(_SIMPLE)


@pytest.fixture(scope="module")
def sweep(runner: NetworkRunner) -> pd.DataFrame:
    return runner.sweep(pd.DataFrame({"inlet.m_dot": [0.5, 1.0, 1.5]}))


def _export_csv(result, path: Path) -> Path:
    # Exactly what the GUI's /export route writes (main.export_results).
    _neutralise_formulas(result.to_dataframe()).to_csv(path, index=False)
    return path


def test_sweep_shape(sweep: pd.DataFrame) -> None:
    assert list(sweep.index) == [0, 1, 2]
    assert sweep.columns.nlevels == 3
    assert list(sweep.columns.names) == ["kind", "entity", "quantity"]
    assert set(sweep.columns.get_level_values(0)) == {"inputs", "outputs"}
    # inputs first, then outputs
    kinds = list(sweep.columns.get_level_values(0))
    assert kinds == sorted(kinds)
    out = sweep["outputs"]
    assert out.columns.nlevels == 2
    for i, m in enumerate([0.5, 1.0, 1.5]):
        assert math.isclose(out["chan", "m_dot"].iloc[i], m, rel_tol=1e-6)
        assert sweep["inputs"]["inlet", "m_dot"].iloc[i] == m
    assert out[SOLVE_ENTITY, "success"].all()


def test_sweep_index_names_the_cases(runner: NetworkRunner) -> None:
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0]}, index=["low", "high"])
    df = runner.sweep(params, metrics=["chan.m_dot"])
    assert list(df.index) == ["low", "high"]
    assert df.attrs["warm_start"] == {"low": None, "high": "low"}


def test_sweep_raising_case_keeps_its_inputs(runner: NetworkRunner) -> None:
    params = pd.DataFrame({"no_such_node.m_dot": [1.0]})
    with pytest.warns(UserWarning, match="raised an exception"):
        df = runner.sweep(params, metrics=["chan.m_dot"])
    assert len(df) == 1
    assert not df["outputs"][SOLVE_ENTITY, "success"].iloc[0]
    assert df["outputs"][SOLVE_ENTITY, "outcome"].iloc[0] == "exception"
    assert math.isnan(df["outputs"]["chan", "m_dot"].iloc[0])
    assert df["inputs"]["chan", "L"].iloc[0] == 1.0


def _assert_same_table(a: pd.DataFrame, b: pd.DataFrame, rel_tol: float = 1e-12) -> None:
    """Same rows, columns (labels), values, units and solver state."""
    assert list(a.index) == list(b.index)
    assert list(a.columns) == list(b.columns)
    for col in a.columns:
        for x, y in zip(a[col], b[col], strict=True):
            if isinstance(x, float) and isinstance(y, float):
                assert (math.isnan(x) and math.isnan(y)) or math.isclose(x, y, rel_tol=rel_tol), col
            else:
                assert x == y, (col, x, y)
    assert a.attrs == b.attrs


def test_export_csv_round_trip(runner: NetworkRunner, tmp_path: Path) -> None:
    # The export carries the inputs (#784): the overridden m_dot comes back
    # from the CSV alone, as do the units, initial guesses and solver options.
    result = runner.solve({"inlet.m_dot": 0.7})
    direct = load_results(result)
    loaded = load_results(_export_csv(result, tmp_path / "x.csv"))
    _assert_same_table(direct, loaded)
    assert loaded["inputs"].loc[0, ("inlet", "m_dot")] == 0.7
    assert loaded.attrs["units"][("inputs", "chan", "L")] == "m"
    assert loaded.attrs["initial_guesses"][0] == {
        "chan": {"m_dot": 0.9},
        "outlet": {"Pt": 101325.0},
    }
    assert loaded.attrs["solver_options"][0]["method"] == "hybr"


def test_export_route_round_trip() -> None:
    # The GUI's own /export response reads back to the runner's table.
    r = http._post("/export", _SIMPLE)
    assert r.status_code == 200, r.text[:300]
    from_route = load_results(pd.read_csv(io.StringIO(r.text), dtype={"id": str, "label": str}))
    from_runner = load_results(NetworkRunner.from_dict(_SIMPLE).solve())
    _assert_same_table(from_runner, from_route, rel_tol=1e-9)


def test_export_keeps_the_774_output_columns(runner: NetworkRunner) -> None:
    # Inputs are appended after every result column, on a new network row:
    # the unit-labelled output headers and the entity rows are as before.
    long = runner.solve().to_dataframe()
    cols = list(long.columns)
    first_input = next(i for i, c in enumerate(cols) if c.startswith("input:"))
    assert all(c.startswith(("input:", "initial_guess:", "solver:")) for c in cols[first_input:])
    assert "T [K]" in cols[:first_input]
    assert list(long["type"])[-1] == "network"
    assert set(list(long["type"])[:-1]) == {"node", "element"}
    net = long[long["type"] == "network"].iloc[0]
    assert net["id"] == "_network"
    assert net["solver:method"] == "hybr"
    assert net["input:global_regime"] in ("compressible", "incompressible")


def test_table_csv_round_trip(sweep: pd.DataFrame, tmp_path: Path) -> None:
    path = tmp_path / "table.csv"
    sweep.to_csv(path)
    back = load_results(path)
    assert list(back.index) == list(sweep.index)
    assert list(back.columns) == list(sweep.columns)
    pd.testing.assert_series_equal(
        back["outputs"]["chan", "m_dot"], sweep["outputs"]["chan", "m_dot"], check_names=False
    )
    assert back.attrs["units"][("outputs", "chan", "m_dot")] == "kg/s"


def test_append_case_reads_the_full_inputs(sweep: pd.DataFrame, runner: NetworkRunner) -> None:
    result = runner.solve({"inlet.m_dot": 2.0})
    df = append_case(sweep, result)
    assert list(df.index) == [0, 1, 2, 3]
    assert len(sweep) == 3  # not changed in place
    row = df["inputs"].loc[3]
    assert row["inlet", "m_dot"] == 2.0  # the override
    assert row["chan", "L"] == 1.0  # geometry, not overridden
    assert row["outlet", "Pt"] == 101325.0
    assert row["inlet", "composition.source"] == "dry_air"
    assert row[NETWORK_ENTITY, "global_regime"] in ("compressible", "incompressible")
    assert math.isclose(df["outputs"].loc[3, ("chan", "m_dot")], 2.0, rel_tol=1e-6)


def test_append_case_name_and_collision(sweep: pd.DataFrame, runner: NetworkRunner) -> None:
    result = runner.solve()
    df = append_case(sweep, result, name="design")
    assert list(df.index) == [0, 1, 2, "design"]
    with pytest.raises(ValueError, match="already has a case named 'design'"):
        append_case(df, result, name="design")


def test_append_case_to_empty_frame(runner: NetworkRunner) -> None:
    df = append_case(pd.DataFrame(), runner.solve())
    df = append_case(df, runner.solve({"inlet.m_dot": 0.4}))
    assert list(df.index) == [0, 1]
    assert list(varied_inputs(df).columns) == [("inlet", "m_dot")]


def test_varied_inputs(sweep: pd.DataFrame) -> None:
    varied = varied_inputs(sweep)
    assert list(varied.columns) == [("inlet", "m_dot")]
    assert list(varied["inlet", "m_dot"]) == [0.5, 1.0, 1.5]


def test_units_in_attrs_not_in_names(sweep: pd.DataFrame) -> None:
    units = sweep.attrs["units"]
    assert units[("outputs", "chan", "m_dot")] == "kg/s"
    assert units[("outputs", "outlet", "T")] == "K"
    assert units[("outputs", "outlet", "rho")] == "kg/m^3"
    assert units[("inputs", "inlet", "m_dot")] == "kg/s"
    assert all(" [" not in q for q in sweep.columns.get_level_values(2))


def test_guesses_and_solver_options_are_not_inputs(sweep: pd.DataFrame) -> None:
    quantities = set(sweep["inputs"].columns.get_level_values(1))
    assert not any(q.startswith("initial_guess") for q in quantities)
    assert not quantities & {"method", "init_strategy", "timeout", "incompressible_mach_limit"}
    assert sweep.attrs["initial_guesses"][0] == {
        "chan": {"m_dot": 0.9},
        "outlet": {"Pt": 101325.0},
    }
    assert sweep.attrs["solver_options"][0]["method"] == "hybr"
    assert sweep.attrs["warm_start"] == {0: None, 1: 0, 2: 1}


def test_duplicate_labels_refused() -> None:
    net = copy.deepcopy(_SIMPLE)
    net["nodes"][2]["data"]["label"] = "chan"
    with pytest.raises(
        ValueError, match=r"named 'chan' \(ids '(n_chan|outlet)' and '(n_chan|outlet)'"
    ):
        NetworkRunner.from_dict(net).sweep(pd.DataFrame({"inlet.m_dot": [1.0]}))


def test_label_falls_back_to_id(runner: NetworkRunner) -> None:
    net = copy.deepcopy(_SIMPLE)
    del net["nodes"][1]["data"]["label"]
    df = load_results(NetworkRunner.from_dict(net).solve())
    assert ("outputs", "n_chan", "m_dot") in df.columns
    assert ("inputs", "n_chan", "L") in df.columns


def test_formula_label_survives_the_export(tmp_path: Path) -> None:
    net = copy.deepcopy(_SIMPLE)
    net["nodes"][1]["data"]["label"] = "=chan"
    result = NetworkRunner.from_dict(net).solve()
    df = load_results(_export_csv(result, tmp_path / "x.csv"))
    assert ("outputs", "=chan", "m_dot") in df.columns


def test_legacy_sweep_frame(runner: NetworkRunner) -> None:
    # A 0.7 full-detail sweep frame: to_dataframe() stacked per case.
    parts = []
    for i, m in enumerate([0.5, 1.0]):
        detail = runner.solve({"inlet.m_dot": m}).to_dataframe()
        detail.insert(0, "inlet.m_dot", m)
        detail.insert(0, "_sweep_index", i)
        parts.append(detail)
    df = load_results(pd.concat(parts, ignore_index=True))
    assert list(df.index) == [0, 1]
    assert list(df["inputs"]["inlet", "m_dot"]) == [0.5, 1.0]
    assert math.isclose(df["outputs"].loc[1, ("chan", "m_dot")], 1.0, rel_tol=1e-6)


# ---------------------------------------------------------------------------
# Units of the inputs
# ---------------------------------------------------------------------------

_NETWORKS = list(gui_result_keys.networks())


@pytest.mark.parametrize(
    "graph",
    [g for _, _, g in _NETWORKS],
    ids=[f"{kind}-{variant}" for kind, variant, _ in _NETWORKS],
)
def test_every_numeric_input_column_has_a_unit(graph: dict) -> None:
    """Every input column of the export with a number in it carries its unit
    (as test_every_export_column_has_a_unit does for the outputs). Text-valued
    inputs stay bare. A new numeric field fails here: add it to
    INPUT_UNIT_MAP in gui/backend/units.py."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = http._post("/export", graph)
    assert r.status_code == 200, r.text[:300]
    rows = list(csv.DictReader(io.StringIO(r.text)))
    unlabelled = set()
    for col in rows[0]:
        if not col.startswith("input:") or col.endswith("]"):
            continue
        for row in rows:
            cell = row[col]
            if cell in ("True", "False"):
                unlabelled.add(col)
                continue
            try:
                float(cell)
            except ValueError:
                continue
            unlabelled.add(col)
    assert unlabelled == set()


_NOT_INPUT_FIELDS = {"timeout", "incompressible_mach_limit", "initial_guess"}
_RESULT_MODELS = {
    "StateResult",
    "NodeResult",
    "ElementResult",
    "ConvergencePoint",
    "WorstResidual",
    "SecondLawFlag",
    "IncompressibleFlag",
    "SolveResponse",
    "NodePosition",
}


def test_every_numeric_schema_field_has_an_input_unit() -> None:
    """Every numeric field a network can carry -- not only those the
    representative networks set -- has an input unit."""
    missing = []
    for name, cls in inspect.getmembers(schemas, inspect.isclass):
        if not issubclass(cls, schemas.BaseModel) or name in _RESULT_MODELS:
            continue
        for field, info in cls.model_fields.items():
            args = typing.get_args(info.annotation) or (info.annotation,)
            if field in _NOT_INPUT_FIELDS or not set(args) & {int, float, bool}:
                continue
            if input_unit(field) is None:
                missing.append(f"{name}.{field}")
    assert missing == []


def test_input_units_of_nested_keys() -> None:
    assert input_unit("composition.ambient_T") == "K"
    assert input_unit("composition.custom_fractions.CH4") == "-"
    assert input_unit("layers.0.thickness") == "m"
    assert input_unit("surface.alpha_deg") == "deg"
    assert input_unit("k") == "-"  # a discrete loss's slope, not a conductivity
    assert input_unit("composition.source") is None
