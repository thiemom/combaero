#!/bin/bash
# Make the checkout's GUI launchable: rebuild the C++ extension when its
# sources changed, then check that the backend imports.
#
# The install is editable, so Python changes are live, but the compiled
# combaero._core is not: after a pull that touched C++, `uv run` started the
# backend against the old build and it died at import (2026-10-17: "cannot
# import name 'mixer_heat_per_unit_flow' from 'combaero._core'").
#
#   ./scripts/ensure-gui-build.sh          # called by run-gui.sh
#   ./scripts/ensure-gui-build.sh --check  # report only, never rebuild

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

VENV_PYTHON="./.venv/bin/python"
CHECK_ONLY=0
[ "$1" = "--check" ] && CHECK_ONLY=1

# What the compiled extension is built from. Python sources are not listed:
# the editable install reads them from the tree.
CPP_INPUTS=(src include python/combaero/_core.cpp CMakeLists.txt pyproject.toml)

STALE=""
CORE_SO="$(ls .venv/lib/python3*/site-packages/combaero/_core*.so 2>/dev/null | head -n 1 || true)"
if [ -z "$CORE_SO" ]; then
    STALE="no compiled extension found"
else
    NEWER="$(find "${CPP_INPUTS[@]}" -type f -newer "$CORE_SO" -print 2>/dev/null | head -n 1)"
    [ -n "$NEWER" ] && STALE="$NEWER is newer than the build"
fi

if [ -n "$STALE" ]; then
    if [ "$CHECK_ONLY" = 1 ]; then
        echo "C++ extension is stale: $STALE"
        exit 2
    fi
    echo "C++ extension is stale ($STALE): rebuilding, about a minute..."
    # --python pins the checkout's own venv: an inherited VIRTUAL_ENV (a
    # worktree's, say) would otherwise receive the build.
    uv pip install -e . --python "$VENV_PYTHON"
fi

# The real test: the backend imports against this build. combaero's
# __init__ catches a failed _core import, prints "Import failed" and carries
# on with a partial module, so that message counts as a failure too.
OUT="$(mktemp)"
trap 'rm -f "$OUT"' EXIT
if ! PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}" "$VENV_PYTHON" -c "import gui.backend.main" >"$OUT" 2>&1 \
    || grep -q "Import failed" "$OUT"; then
    echo "Error: the GUI backend does not import cleanly:"
    tail -n 5 "$OUT"
    echo "Try: uv pip install -e . (C++ extension) and uv sync (dependencies)."
    exit 1
fi
echo "C++ extension up to date; backend imports."
