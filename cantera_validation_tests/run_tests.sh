#!/usr/bin/env bash
set -euo pipefail

# Helper to run validation tests in an isolated environment.
# This uses the root-managed 'uv' workspace environment.

echo "=== Running Validation Tests (via uv) ==="
# Use -n auto for parallel execution
# We use --project to ensure we're targeting the workspace member specifically
uv run --project cantera_validation_tests pytest -v --tb=short
