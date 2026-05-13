#!/bin/bash
set -euo pipefail

# GUI Code Style Checker (Biome)
# ===============================
# Runs Biome linter and formatter on the frontend TypeScript/JavaScript code.
#
# Usage:
#   scripts/check-gui-style.sh [--fix]
#
# Options:
#   --fix   Auto-fix issues where possible (passes --write to biome)
#
# Default (no flags) is check-only — matches the CI lint-gui job.

cd "$(dirname "$0")/.."

FIX=false
for arg in "$@"; do
    case "$arg" in
        --fix) FIX=true ;;
        --help|-h)
            echo "Usage: $0 [--fix]"
            echo "  --fix   Auto-fix issues (biome check --write)"
            exit 0
            ;;
    esac
done

echo "======================================"
echo "===== GUI Code Style Checker (Biome)"
echo "======================================"
echo ""

FRONTEND_DIR="$(pwd)/gui/frontend"

if [ ! -d "$FRONTEND_DIR" ]; then
    echo "ERROR: gui/frontend not found. Run from the repo root."
    exit 1
fi

run_biome() {
    # pnpm exec must run from the package directory.
    if command -v pnpm > /dev/null 2>&1; then
        (cd "$FRONTEND_DIR" && pnpm exec biome "$@")
    elif [ -f "/Users/thiemo/.nvm/versions/node/v22.21.0/bin/npm" ]; then
        (cd "$FRONTEND_DIR" && "/Users/thiemo/.nvm/versions/node/v22.21.0/bin/npm" exec biome "$@")
    else
        echo "ERROR: neither pnpm nor npm found. Install Node/pnpm first."
        exit 1
    fi
}

if [ "$FIX" = true ]; then
    echo "Running Biome check with --write..."
    run_biome check ../frontend/src --write
    echo ""
    echo "✓ Biome check complete (fixes applied)"
else
    echo "Running Biome check (read-only)..."
    run_biome check ../frontend/src
    echo ""
    echo "✓ Biome check passed"
fi
