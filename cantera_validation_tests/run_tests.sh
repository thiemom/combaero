#!/bin/bash
set -e

cd "$(dirname "$0")"

echo "=== Installing dependencies with Poetry ==="
poetry install --quiet

echo ""
echo "=== Running Cantera Validation Tests ==="
poetry run pytest -v --tb=short

echo ""
echo "=== All tests passed! ==="
