"""Pandas 3.0 compatibility tests for the GUI CSV export pipeline.

Exercises the exact DataFrame construction, column-rename, and to_csv code
paths in gui/backend/main.py::export_results without requiring an HTTP layer.
"""

import io

import pandas as pd

import combaero as cb
from gui.backend.schemas import ElementResult, NodeResult, StateResult
from gui.backend.units import label_with_unit

_N_SPECIES = len(cb.species.names)
_SPECIES_NAMES = list(cb.species.names)
_Y_AIR = [0.0] * _N_SPECIES
_Y_AIR[0] = 0.768  # N2
_Y_AIR[1] = 0.232  # O2


def _node(T: float = 300.0, P: float = 101325.0) -> NodeResult:
    return NodeResult(
        state=StateResult(T=T, P=P, Pt=P * 1.01, Tt=T, rho=1.18, mach=0.1, Y=_Y_AIR),
        success=True,
    )


def _build_data(
    node_results: dict,
    element_results: dict,
    edge_results: dict,
    id_to_label: dict,
    id_to_kind: dict,
) -> list[dict]:
    """Mirror the data-assembly logic in export_results."""
    _BOUNDARY_KINDS = {"mass_boundary", "pressure_boundary"}
    data: list[dict] = []

    for node_id, res in node_results.items():
        state_data = res.state.model_dump()
        Y = state_data.pop("Y", [])
        X = state_data.pop("X", None)
        kind = id_to_kind.get(node_id, "")
        row: dict = {
            "type": "node",
            "kind": kind,
            "is_boundary": kind in _BOUNDARY_KINDS,
            "combustion_method": "",
            "id": node_id,
            "label": id_to_label.get(node_id, ""),
            **state_data,
        }
        for i, y_val in enumerate(Y):
            name = _SPECIES_NAMES[i] if i < len(_SPECIES_NAMES) else str(i)
            row[f"Y[{name}]"] = y_val
        if X:
            for i, x_val in enumerate(X):
                name = _SPECIES_NAMES[i] if i < len(_SPECIES_NAMES) else str(i)
                row[f"X[{name}]"] = x_val
        data.append(row)

    for elem_id, res in element_results.items():
        data.append(
            {
                "type": "element",
                "kind": id_to_kind.get(elem_id, ""),
                "is_boundary": False,
                "combustion_method": "",
                "id": elem_id,
                "label": id_to_label.get(elem_id, ""),
                **res.model_dump(),
            }
        )

    for wall_id, res in edge_results.items():
        data.append(
            {
                "type": "thermal_wall",
                "kind": "thermal_wall",
                "is_boundary": False,
                "combustion_method": "",
                "id": wall_id,
                "label": id_to_label.get(wall_id, ""),
                **res,
            }
        )
    return data


def test_dataframe_from_mixed_rows_has_correct_shape() -> None:
    node_results = {
        "n_in": _node(T=300.0, P=200000.0),
        "n_out": _node(T=295.0, P=101325.0),
    }
    element_results = {"e1": ElementResult(m_dot=0.5, success=True)}
    data = _build_data(
        node_results,
        element_results,
        {},
        {"n_in": "Inlet", "n_out": "Outlet", "e1": ""},
        {"n_in": "pressure_boundary", "n_out": "pressure_boundary", "e1": "channel"},
    )
    df = pd.DataFrame(data)
    assert len(df) == 3
    assert "T" in df.columns
    assert "P" in df.columns
    assert "m_dot" in df.columns


def test_column_rename_applies_unit_labels() -> None:
    data = [
        {
            "T": 300.0,
            "P": 101325.0,
            "id": "n1",
            "type": "node",
            "kind": "plenum",
            "is_boundary": False,
            "combustion_method": "",
            "label": "",
        }
    ]
    df = pd.DataFrame(data)
    df = df.rename(columns={col: label_with_unit(col) for col in df.columns})
    assert "T [K]" in df.columns
    assert "P [Pa]" in df.columns
    assert "id" in df.columns


def test_to_csv_produces_valid_output() -> None:
    node_results = {"n_in": _node(T=300.0, P=200000.0)}
    element_results = {"e1": ElementResult(m_dot=1.23, success=True)}
    data = _build_data(
        node_results,
        element_results,
        {},
        {"n_in": "Inlet", "e1": ""},
        {"n_in": "pressure_boundary", "e1": "channel"},
    )
    df = pd.DataFrame(data)
    df = df.rename(columns={col: label_with_unit(col) for col in df.columns})
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    csv_text = stream.getvalue()

    lines = [ln for ln in csv_text.strip().split("\n") if ln]
    assert len(lines) == 3  # header + 2 data rows
    assert "P [Pa]" in lines[0]
    assert "T [K]" in lines[0]
    assert "m_dot [kg/s]" in lines[0]
    assert "200000" in csv_text
    assert "1.23" in csv_text


def test_species_composition_columns_use_name_index() -> None:
    """Y[N2], Y[O2] etc. (not Y[0]) appear in the exported DataFrame."""
    data = _build_data(
        {"n1": _node()},
        {},
        {},
        {"n1": ""},
        {"n1": "plenum"},
    )
    df = pd.DataFrame(data)
    assert f"Y[{_SPECIES_NAMES[0]}]" in df.columns
    assert "Y[0]" not in df.columns


def test_empty_edge_results_produces_no_thermal_wall_rows() -> None:
    data = _build_data(
        {"n1": _node()},
        {},
        {},
        {"n1": ""},
        {"n1": "plenum"},
    )
    df = pd.DataFrame(data)
    assert "type" in df.columns
    assert (df["type"] == "thermal_wall").sum() == 0


def test_meta_columns_pass_through_rename_unchanged() -> None:
    meta = ["id", "type", "kind", "is_boundary", "combustion_method", "label"]
    data = [dict.fromkeys(meta, "value")]
    df = pd.DataFrame(data)
    df = df.rename(columns={col: label_with_unit(col) for col in df.columns})
    for col in meta:
        assert col in df.columns, f"meta column '{col}' was unexpectedly renamed"
