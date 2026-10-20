"""One max-solve-time default for the GUI and NetworkRunner.

The frontend store has defaulted to 30 s while the backend schema and
NetworkRunner fell back to 180 s, so a request without solver settings (an
older saved file, a script, a direct API call) got six times the GUI's budget.
The default lives in ``gui.backend.schemas.DEFAULT_SOLVE_TIMEOUT``.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

from gui.backend.runner import NetworkRunner
from gui.backend.schemas import DEFAULT_SOLVE_TIMEOUT, SolverSettings

_STORE = Path(__file__).resolve().parents[2] / "gui" / "frontend" / "src" / "store" / "useStore.ts"


def test_the_default_is_thirty_seconds() -> None:
    assert DEFAULT_SOLVE_TIMEOUT == 30.0


def test_gui_backend_default_matches() -> None:
    assert SolverSettings().timeout == DEFAULT_SOLVE_TIMEOUT


def test_frontend_store_default_matches() -> None:
    text = _STORE.read_text(encoding="utf-8")
    # The value assignment (`timeout: 30,`), not the type (`timeout: number | null;`).
    found = re.findall(r"^\s*timeout:\s*([0-9.]+),", text, flags=re.MULTILINE)
    assert len(found) == 1, found
    assert float(found[0]) == DEFAULT_SOLVE_TIMEOUT


def test_network_runner_defaults_match() -> None:
    for method in (NetworkRunner.solve, NetworkRunner.sweep):
        default = inspect.signature(method).parameters["timeout"].default
        assert default == DEFAULT_SOLVE_TIMEOUT, f"{method.__qualname__}: {default!r}"
