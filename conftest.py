"""Repo-root conftest: make `validation/` importable from pytest.

The `validation/` tree is a private in-tree runner/dataset, not a published
package. This conftest adds the repo root to sys.path so tests under
python/tests/ can `from validation.junction... import ...`.
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
