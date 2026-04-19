#!/bin/bash
set -euo pipefail

# GUI Code Style Checker (Biome)
# ===============================
# Runs Biome linter and formatter on the frontend TypeScript/JavaScript code.
# Automatically fixes issues where possible.

cd "$(dirname "$0")/.."

echo "======================================"
echo "===== GUI Code Style Checker (Biome)"
echo "======================================"
echo ""

# Use the npm from the user's NVM installation
NPM_PATH="/Users/thiemo/.nvm/versions/node/v22.21.0/bin/npm"

if [ ! -f "$NPM_PATH" ]; then
    echo "ERROR: npm not found at $NPM_PATH"
    echo "Please update NPM_PATH in this script."
    exit 1
fi

# Run Biome check with --fix to auto-fix issues
echo "Running Biome check with --fix..."
"$NPM_PATH" exec biome check gui/frontend/src --write

echo ""
echo "✓ Biome check complete"
