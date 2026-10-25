"""Every runtime dependency of combaero and combaero-gui declares a floor.

The nightly dependency-floors job installs each one at its floor and runs the
suite (scripts/dependency_floors.py). A requirement with no lower bound would
resolve to the first release ever published, so this check runs on every PR
rather than waiting for the night.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

from packaging.requirements import Requirement

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "dependency_floors", REPO / "scripts" / "dependency_floors.py"
)
floors = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(floors)


def test_every_runtime_dependency_has_a_floor() -> None:
    reqs = floors.requirements()
    assert reqs
    assert floors.unbounded(reqs) == []


def test_the_check_catches_an_unbounded_requirement() -> None:
    """Falsified: a bare name and an upper bound alone are both unbounded."""
    bad = [
        ("x", Requirement("fastapi")),
        ("x", Requirement("uvicorn<1")),
        ("x", Requirement("h11>=0.16")),
    ]
    assert floors.unbounded(bad) == ["x: fastapi", "x: uvicorn<1"]


def test_the_gui_web_stack_clears_the_known_advisories() -> None:
    """The floors that exist for a security reason do not drift down."""
    from packaging.version import Version

    want = {"starlette": "1.3.1", "h11": "0.16.0", "fastapi": "0.133.0", "pydantic": "2.7.0"}
    got = {r.name: floors.floor(r) for _, r in floors.requirements() if r.name in want}
    for name, v in want.items():
        assert Version(got[name]) >= Version(v), name
