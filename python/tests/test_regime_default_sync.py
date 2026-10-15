"""One regime default everywhere (2026-10-15 decision).

A script and a drawing of the same network must solve the same physics. The
default lives in ``combaero.network.components.DEFAULT_REGIME``; the GUI
backend's ``SolverSettings.global_regime`` and the frontend store's
``solverSettings.global_regime`` must equal it, as must every element
constructor that takes a ``regime``.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import combaero.network as network
from combaero.network.components import DEFAULT_REGIME
from gui.backend.schemas import SolverSettings

_STORE = Path(__file__).resolve().parents[2] / "gui" / "frontend" / "src" / "store" / "useStore.ts"


def test_the_default_is_incompressible_until_585() -> None:
    """Compressible is the goal, but with it as the default 60 tests failed
    to converge (2026-10-15). Flip this, and DEFAULT_REGIME, when #585's
    acceptance test passes."""
    assert DEFAULT_REGIME == "incompressible"


def test_gui_backend_default_matches() -> None:
    assert SolverSettings().global_regime == DEFAULT_REGIME


def test_frontend_store_default_matches() -> None:
    text = _STORE.read_text(encoding="utf-8")
    # The value assignment (`global_regime: "x",`), not the type declaration
    # (`global_regime: "incompressible" | "compressible";`).
    found = re.findall(r'global_regime:\s*"(\w+)",', text)
    assert len(found) == 1, found
    m = re.match(r"(\w+)", found[0])
    assert m is not None, "solverSettings.global_regime default not found in useStore.ts"
    assert m.group(1) == DEFAULT_REGIME


def test_every_element_regime_default_matches() -> None:
    checked = []
    for name in dir(network):
        obj = getattr(network, name)
        if not inspect.isclass(obj):
            continue
        param = inspect.signature(obj.__init__).parameters.get("regime")
        if param is None or param.default is inspect.Parameter.empty:
            continue
        checked.append(name)
        assert param.default == DEFAULT_REGIME, f"{name}.regime defaults to {param.default!r}"
    assert {"ChannelElement", "OrificeElement", "MomentumChamberNode"} <= set(checked)
