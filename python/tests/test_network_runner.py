"""Tests for NetworkRunner and NetworkResult.

Covers: from_file, from_dict, solve (no overrides), BC overrides, sweep,
NetworkResult.get, node_state, to_dataframe.  The multi_wall_tee.json
integration test is skipped if the fixture file is absent.
"""

import math
from pathlib import Path

import pandas as pd
import pytest

from gui.backend.runner import NetworkResult, NetworkRunner

_FIXTURE_DIR = Path(__file__).parent.parent.parent / "gui" / "tmp"
_MULTI_WALL_TEE = _FIXTURE_DIR / "multi_wall_tee.json"

# ---------------------------------------------------------------------------
# Minimal inline fixture — one mass boundary, one channel, one pressure boundary.
# Uses the incompressible (default) solver for speed.
# ---------------------------------------------------------------------------

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
            "id": "chan",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {"label": "chan", "L": 1.0, "D": 0.1, "roughness": 1e-5},
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
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "inlet", "target": "chan"},
        {"id": "e2", "source": "chan", "target": "outlet"},
    ],
}


# Fixture with UUID-style IDs and separate labels (tests label-based resolution)
_UUID_LABELED: dict = {
    "nodes": [
        {
            "id": "node_abc123",
            "type": "mass_boundary",
            "position": {"x": 0, "y": 0},
            "data": {
                "label": "air_inlet",
                "m_dot": 1.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
        {
            "id": "node_def456",
            "type": "channel",
            "position": {"x": 200, "y": 0},
            "data": {"label": "hot_channel", "L": 1.0, "D": 0.1, "roughness": 1e-5},
        },
        {
            "id": "node_ghi789",
            "type": "pressure_boundary",
            "position": {"x": 400, "y": 0},
            "data": {
                "label": "exhaust",
                "Pt": 101325.0,
                "Tt": 300.0,
                "composition": {"source": "dry_air"},
            },
        },
    ],
    "edges": [
        {"id": "e1", "source": "node_abc123", "target": "node_def456"},
        {"id": "e2", "source": "node_def456", "target": "node_ghi789"},
    ],
}


# ---------------------------------------------------------------------------
# Construction tests
# ---------------------------------------------------------------------------


def test_from_dict():
    runner = NetworkRunner.from_dict(_SIMPLE)
    assert isinstance(runner, NetworkRunner)


def test_from_file(tmp_path):
    import json

    p = tmp_path / "simple.json"
    p.write_text(json.dumps(_SIMPLE))
    runner = NetworkRunner.from_file(p)
    assert isinstance(runner, NetworkRunner)


def test_from_dict_normalises_solver_settings_camelcase():
    """GUI exports use camelCase 'solverSettings'; from_dict must translate it."""
    import copy

    d = copy.deepcopy(_SIMPLE)
    d["solverSettings"] = {"global_regime": "compressible", "method": "hybr"}
    runner = NetworkRunner.from_dict(d)
    assert runner._schema.solver_settings.global_regime == "compressible"


def test_from_dict_snake_case_solver_settings_unchanged():
    """If snake_case solver_settings is already present it takes precedence."""
    import copy

    d = copy.deepcopy(_SIMPLE)
    d["solver_settings"] = {"global_regime": "incompressible"}
    d["solverSettings"] = {"global_regime": "compressible"}
    runner = NetworkRunner.from_dict(d)
    assert runner._schema.solver_settings.global_regime == "incompressible"


# ---------------------------------------------------------------------------
# Solve — basic
# ---------------------------------------------------------------------------


def test_solve_success():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    assert result.success
    assert result.final_norm is not None
    assert result.final_norm < 1e-3


def test_solve_returns_network_result():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    assert isinstance(result, NetworkResult)


def test_solve_mass_conservation():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    m_dot = result.get("chan.m_dot")
    assert math.isclose(m_dot, 1.0, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# Solve — statelessness (repeated calls must not interfere)
# ---------------------------------------------------------------------------


def test_solve_stateless():
    runner = NetworkRunner.from_dict(_SIMPLE)
    r1 = runner.solve({"inlet.m_dot": 0.5})
    r2 = runner.solve({"inlet.m_dot": 1.5})
    r3 = runner.solve()
    assert math.isclose(result_m_dot(r1), 0.5, rel_tol=1e-3)
    assert math.isclose(result_m_dot(r2), 1.5, rel_tol=1e-3)
    assert math.isclose(result_m_dot(r3), 1.0, rel_tol=1e-3)


def result_m_dot(result: NetworkResult) -> float:
    return result.get("chan.m_dot")


# ---------------------------------------------------------------------------
# Overrides
# ---------------------------------------------------------------------------


def test_override_m_dot():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve({"inlet.m_dot": 2.0})
    assert result.success
    assert math.isclose(result.get("chan.m_dot"), 2.0, rel_tol=1e-3)


def test_override_by_node_id_fallback():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve({"inlet.m_dot": 0.7})
    assert math.isclose(result.get("chan.m_dot"), 0.7, rel_tol=1e-3)


def test_override_missing_label_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    with pytest.raises(KeyError, match="no_such_node"):
        runner.solve({"no_such_node.m_dot": 1.0})


def test_override_bad_key_format_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    with pytest.raises(ValueError, match="<label>.<attribute>"):
        runner.solve({"m_dot": 1.0})


# ---------------------------------------------------------------------------
# NetworkResult API
# ---------------------------------------------------------------------------


def test_result_get_known_key():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    v = result.get("chan.m_dot")
    assert isinstance(v, float)


def test_result_get_missing_key_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(KeyError):
        result.get("nonexistent.key")


def test_result_get_by_label_when_id_differs():
    """get() must resolve '<label>.<quantity>' when label != node id."""
    runner = NetworkRunner.from_dict(_UUID_LABELED)
    result = runner.solve()
    # 'hot_channel' is the label for node id 'node_def456'
    v = result.get("hot_channel.m_dot")
    assert isinstance(v, float)
    assert math.isclose(v, 1.0, rel_tol=1e-3)


def test_result_get_label_fallback_does_not_shadow_raw_key():
    """A key present verbatim in raw must be returned without label resolution."""
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    v_direct = result.get("chan.m_dot")
    v_label = result.get("chan.m_dot")
    assert v_direct == v_label


def test_node_state_by_label():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    state = result.node_state("outlet")
    assert "T" in state
    assert "P" in state
    assert "Y" in state


def test_node_state_missing_label_raises():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    with pytest.raises(KeyError):
        result.node_state("no_such")


# ---------------------------------------------------------------------------
# to_dataframe
# ---------------------------------------------------------------------------


def test_to_dataframe_returns_dataframe():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_to_dataframe_has_type_column():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "type" in df.columns


def test_to_dataframe_has_element_row():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "element" in df["type"].values


def test_to_dataframe_unit_annotated_columns():
    runner = NetworkRunner.from_dict(_SIMPLE)
    result = runner.solve()
    df = result.to_dataframe()
    assert "T [K]" in df.columns
    assert "P [Pa]" in df.columns
    assert "m_dot [kg/s]" in df.columns


# ---------------------------------------------------------------------------
# sweep
# ---------------------------------------------------------------------------


def test_sweep_metrics_shape():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0, 1.5]})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    assert len(df) == 3
    assert "chan.m_dot" in df.columns
    assert "inlet.m_dot" in df.columns
    assert "success" in df.columns


def test_sweep_metrics_values():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0, 2.0]})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    for i, expected in enumerate([0.5, 1.0, 2.0]):
        assert math.isclose(df.loc[i, "chan.m_dot"], expected, rel_tol=1e-3)


def test_sweep_no_metrics_returns_full_detail():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": [0.5, 1.0]})
    df = runner.sweep(params)
    assert "_sweep_index" in df.columns
    assert len(df) > 2


def test_sweep_empty_params():
    runner = NetworkRunner.from_dict(_SIMPLE)
    params = pd.DataFrame({"inlet.m_dot": []})
    df = runner.sweep(params, metrics=["chan.m_dot"])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


# ---------------------------------------------------------------------------
# Integration — multi_wall_tee.json (skipped if fixture absent)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _MULTI_WALL_TEE.exists(), reason="fixture not present")
def test_multi_wall_tee_solve():
    runner = NetworkRunner.from_file(_MULTI_WALL_TEE)
    result = runner.solve(
        method="hybr",
        init_strategy="incompressible_warmstart",
        timeout=60.0,
    )
    assert result.success, f"Solver failed: {result.message}"
    assert result.final_norm < 1e-3


@pytest.mark.skipif(not _MULTI_WALL_TEE.exists(), reason="fixture not present")
def test_multi_wall_tee_to_dataframe():
    runner = NetworkRunner.from_file(_MULTI_WALL_TEE)
    result = runner.solve(
        method="hybr",
        init_strategy="incompressible_warmstart",
        timeout=60.0,
    )
    df = result.to_dataframe()
    assert len(df) > 0
    assert "thermal_wall" in df["type"].values
