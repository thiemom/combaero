"""The GUI's list of result keys is what the backend emits.

``gui/frontend/src/utils/resultKeys.generated.json`` (written by
``scripts/gui_result_keys.py``) is the frontend's view of the backend: the
Probe offers a kind's keys from it before a solve, a key with no catalogue
entry shows its SI unit from it, and ``quantities.test.ts`` checks the
hand-kept ``QUANTITY_CATALOGUE`` against it (no stale entry, no unit at odds
with the backend's). This file keeps the JSON itself honest.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))

import gui_result_keys as keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402

from gui.backend.graph_builder import ELEMENT_KINDS  # noqa: E402

_DOC = json.loads(keys.OUT.read_text(encoding="utf-8"))


def test_every_palette_kind_and_element_kind_is_listed() -> None:
    listed = set(_DOC["kinds"])
    palette = set(http._palette_kinds()) - {"probe"}
    assert palette <= listed, palette - listed
    assert listed >= ELEMENT_KINDS, ELEMENT_KINDS - listed
    assert keys.THERMAL_WALL in listed
    assert all(_DOC["kinds"].values()), [k for k, v in _DOC["kinds"].items() if not v]


def test_every_emitted_key_has_a_unit() -> None:
    """A key with no unit in gui/backend/units.py shows a bare number in the
    GUI and an unlabelled export column: add it to UNIT_MAP."""
    missing = sorted(k for k, u in _DOC["units"].items() if u is None)
    assert missing == []


def test_the_units_are_the_backends() -> None:
    from gui.backend.units import unit_of

    assert {k: unit_of(k) for k in _DOC["units"]} == _DOC["units"]


@pytest.mark.slow
def test_the_file_is_what_the_backend_emits_now() -> None:
    """Rerun ``uv run python scripts/gui_result_keys.py`` after adding,
    renaming or removing a result key (and update QUANTITY_CATALOGUE, which
    quantities.test.ts then checks against the new file)."""
    fresh = keys.collect()
    for kind, emitted in fresh["kinds"].items():
        committed = set(_DOC["kinds"].get(kind, []))
        assert set(emitted) == committed, (
            kind,
            {"new": sorted(set(emitted) - committed), "gone": sorted(committed - set(emitted))},
        )
    assert keys.render(fresh) == keys.OUT.read_text(encoding="utf-8")
