"""combaero.network.diagnostics imports without matplotlib.

matplotlib is the optional `examples` extra. A clean-venv install of the
wheel (the release check) failed to import this module because it imported
pyplot at module scope; only the plotting call needs it.
"""

from __future__ import annotations

import subprocess
import sys

_BLOCKED = """
import sys
sys.modules["matplotlib"] = None          # as if not installed
sys.modules["matplotlib.pyplot"] = None
import combaero.network.diagnostics as d
import tempfile, pathlib
with tempfile.TemporaryDirectory() as t:
    d.export_diagnostic_report({"ok": 1}, pathlib.Path(t) / "r.json")
try:
    d.create_heatmap(__import__("numpy").zeros((2, 2)), "t", __import__("pathlib").Path("x.png"), (2, 2))
except ModuleNotFoundError as e:
    assert "pip install combaero[examples]" in str(e), e
    print("refused-with-hint")
"""


def test_imports_and_exports_without_matplotlib() -> None:
    r = subprocess.run([sys.executable, "-c", _BLOCKED], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "refused-with-hint" in r.stdout
