"""#783: the CSV export writes the displayed result (no second solve), labels
the rows the builder made, and writes list cells as strict JSON."""

from __future__ import annotations

import copy
import io
import json
import math

import pandas as pd
import pytest
import test_gui_http_routes as http

import gui.backend.main as main
from gui.backend.results import load_results


@pytest.fixture(autouse=True)
def _fresh_cache():
    main._SOLVE_CACHE.clear()
    yield
    main._SOLVE_CACHE.clear()


def _export(graph: dict):
    r = http._post("/export", graph)
    assert r.status_code == 200, r.text[:300]
    return r, pd.read_csv(io.StringIO(r.text))


def _no_solve(monkeypatch: pytest.MonkeyPatch) -> None:
    async def refuse(schema):
        raise AssertionError("the export solved again")

    monkeypatch.setattr(main, "_run_solve", refuse)


def test_an_export_after_a_solve_writes_the_displayed_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    graph = http._palette_network("orifice")
    body = http._post("/solve", graph).json()
    assert body["success"], body["message"]
    _no_solve(monkeypatch)
    r, df = _export(graph)
    assert r.headers["x-combaero-export-source"] == "displayed"
    row = df[df["id"] == "x"].iloc[0]
    assert row["m_dot [kg/s]"] == body["element_results"]["x"]["m_dot"]


def test_layout_labels_and_stored_results_do_not_miss_the_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The GUI sends the graph back with each node's result in its data, at
    new positions and with new labels: still the displayed solve."""
    graph = http._palette_network("orifice")
    body = http._post("/solve", graph).json()
    shown = copy.deepcopy(graph)
    for n in shown["nodes"]:
        n["position"] = {"x": 123.0, "y": 456.0}
        n["data"]["label"] = f"renamed {n['id']}"
        n["data"]["result"] = body["element_results"].get(n["id"], {})
    _no_solve(monkeypatch)
    r, df = _export(shown)
    assert r.headers["x-combaero-export-source"] == "displayed"
    # the labels are the current ones
    assert df[df["id"] == "x"].iloc[0]["label"] == "renamed x"


def test_a_changed_input_misses_and_says_it_solved() -> None:
    graph = http._palette_network("orifice")
    http._post("/solve", graph)
    edited = copy.deepcopy(graph)
    for n in edited["nodes"]:
        if n["id"] == "x":
            n["data"]["diameter"] = 0.5 * float(n["data"].get("diameter", 0.08))
    r, _ = _export(edited)
    assert r.headers["x-combaero-export-source"] == "solved"
    # and the next export of that canvas is the one it just solved
    r2, _ = _export(edited)
    assert r2.headers["x-combaero-export-source"] == "displayed"


def test_the_cache_keeps_a_few_solves() -> None:
    for d in (0.02, 0.03, 0.04, 0.05, 0.06):
        g = http._palette_network("orifice")
        for n in g["nodes"]:
            if n["id"] == "x":
                n["data"]["diameter"] = d
        http._post("/solve", g)
    assert len(main._SOLVE_CACHE) == main._SOLVE_CACHE_SIZE


@pytest.mark.parametrize(
    "kind", ["impingement_array", "effusion_liner", "mpce_tee", "compressor", "fuel_boundary"]
)
def test_every_builder_row_names_its_drawn_parent(kind: str) -> None:
    graph = http._palette_network(kind)
    drawn = {n["id"] for n in graph["nodes"]} | {e["id"] for e in graph["edges"]}
    _, df = _export(graph)
    built = df[~df["id"].astype(str).isin(drawn) & (df["type"] != "network")]
    assert len(built) > 0, kind
    for _, row in built.iterrows():
        assert isinstance(row["kind"], str) and row["kind"], row["id"]
        assert isinstance(row["label"], str) and row["label"], row["id"]
        assert row["parent"] in drawn, (row["id"], row["parent"])
    drawn_rows = df[df["id"].astype(str).isin(drawn)]
    assert drawn_rows["parent"].isna().all() or (drawn_rows["parent"].fillna("") == "").all()


def test_list_cells_are_strict_json_and_round_trip() -> None:
    """Plain pandas reads them with json.loads; the results loader reads them
    back; NaN / inf are null."""
    graph = http._palette_network("impingement_array")
    r, df = _export(graph)
    list_cols = [c for c in df.columns if str(c).startswith("rows_")]
    assert list_cols
    df2 = pd.read_csv(
        io.StringIO(r.text),
        converters={c: lambda v: json.loads(v) if v else None for c in list_cols},
    )
    cell = df2[list_cols[0]].dropna().iloc[0]
    assert isinstance(cell, list) and all(x is None or isinstance(x, float) for x in cell)
    for c in list_cols:
        for v in df[c].dropna():
            assert "nan" not in v and "inf" not in v.lower(), (c, v)


def test_the_loader_reads_json_and_old_repr_list_cells(tmp_path) -> None:
    from gui.backend.results import _clean

    assert _clean("[1.0, null, 2.5]", True)[0] == 1.0
    assert math.isnan(_clean("[1.0, null, 2.5]", True)[1])
    assert math.isnan(_clean("[1.0, nan]", True)[1])  # a CSV exported before #783
    graph = http._palette_network("impingement_array")
    r, _ = _export(graph)
    path = tmp_path / "export.csv"
    path.write_text(r.text)
    table = load_results(path)
    assert len(table) == 1


def test_a_non_finite_list_value_is_written_as_null() -> None:
    """Python's repr wrote nan / inf, which JSON does not accept."""
    df = pd.DataFrame({"id": ["a"], "rows_T [K]": [[300.0, float("nan"), float("inf")]]})
    text = main._lists_as_json(df).to_csv(index=False)
    back = pd.read_csv(io.StringIO(text), converters={"rows_T [K]": json.loads})
    assert back["rows_T [K]"].iloc[0] == [300.0, None, None]
    # the frame itself keeps its list
    assert df["rows_T [K]"].iloc[0][0] == 300.0
