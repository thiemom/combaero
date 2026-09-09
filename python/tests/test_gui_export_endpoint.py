"""The /export route itself, not just the DataFrame it builds.

`test_gui_export_dataframe.py` covers the DataFrame construction "without
requiring an HTTP layer", which was exactly the gap: the route that calls it
went untested, and it had been broken since 2026-05-27.

`_solve_sync` gained a sixth return value (`diag`) in #168. The `/solve`
handler was updated; `/export` was not, and kept unpacking five. Every export
that got as far as solving died with

    HTTP 400  {"detail": "too many values to unpack (expected 5)"}

for three and a half months, reaching the user as the GUI's unhelpful
"Export failed. Check console for details."

These tests call the route coroutines directly rather than through
`fastapi.testclient`, which would pull in `httpx` -- not a dependency here, and
not one worth adding for this. The defect lives in the route body, so that is
what these exercise; the HTTP layer around it is FastAPI's to get right.
"""

from __future__ import annotations

import asyncio
import csv
import io

import pytest
from fastapi import HTTPException

from gui.backend.main import export_results, solve
from gui.backend.schemas import NetworkGraphSchema


def _schema() -> NetworkGraphSchema:
    """A minimal solvable network: two pressure boundaries and a channel."""
    return NetworkGraphSchema(
        nodes=[
            {
                "id": "inlet",
                "type": "pressure_boundary",
                "position": {"x": 0.0, "y": 0.0},
                "data": {"Pt": 2.0e5, "Tt": 300.0},
            },
            {
                "id": "ch",
                "type": "channel",
                "position": {"x": 200.0, "y": 0.0},
                "data": {"L": 0.5, "D": 0.05},
            },
            {
                "id": "outlet",
                "type": "pressure_boundary",
                "position": {"x": 400.0, "y": 0.0},
                "data": {"Pt": 1.9e5, "Tt": 300.0},
            },
        ],
        edges=[
            {"id": "e1", "source": "inlet", "target": "ch", "data": {}},
            {"id": "e2", "source": "ch", "target": "outlet", "data": {}},
        ],
    )


def _csv_text(response) -> str:
    async def drain() -> str:
        chunks = [c async for c in response.body_iterator]
        return "".join(c.decode() if isinstance(c, bytes) else c for c in chunks)

    return asyncio.run(drain())


@pytest.fixture(scope="module")
def exported():
    response = asyncio.run(export_results(_schema()))
    return response, _csv_text(response)


def test_export_returns_a_csv(exported):
    """The regression. Before the fix this raised HTTPException(400,
    'too many values to unpack (expected 5)')."""
    response, _ = exported

    assert response.status_code == 200
    assert "text/csv" in response.media_type


def test_the_csv_has_a_header_and_at_least_one_row(exported):
    """A 200 carrying an empty body would still pass the test above."""
    _, text = exported
    rows = list(csv.reader(io.StringIO(text)))

    assert len(rows) >= 2, "header only, no data rows"
    header = rows[0]
    assert "type" in header and "id" in header
    assert any(col.startswith("P [") for col in header), f"no pressure column in {header}"


def test_every_network_entity_reaches_the_csv(exported):
    """The export is per-entity, so a silently truncated frame would still
    look like a valid CSV."""
    _, text = exported
    ids = {r["id"] for r in csv.DictReader(io.StringIO(text))}

    assert {"inlet", "ch", "outlet"} <= ids, f"missing entities: {ids}"


def test_a_filename_is_offered(exported):
    """It is a download; without the disposition header a browser renders the
    CSV instead of saving it."""
    response, _ = exported
    disposition = response.headers.get("content-disposition", "")

    assert "attachment" in disposition and ".csv" in disposition


def test_an_invalid_network_is_a_clean_400_not_a_500():
    """The route's other job: a topology error should reach the user as a
    message rather than an unhandled exception."""
    broken = _schema()
    broken.edges.append(type(broken.edges[0])(id="e3", source="inlet", target="ch", data={}))

    with pytest.raises(HTTPException) as info:
        asyncio.run(export_results(broken))

    assert info.value.status_code == 400
    assert info.value.detail


def test_solve_and_export_agree_on_the_networks_they_accept():
    """The pair that would have caught the original bug. Both routes unpack
    `_solve_sync`; /solve was updated for the new arity and /export was not,
    so they diverged on a network both should handle."""
    schema = _schema()

    solved = asyncio.run(solve(schema))
    exported = asyncio.run(export_results(schema))

    assert solved.success, solved.message
    assert exported.status_code == 200
