#!/usr/bin/env bash
# Check that include/thermo_transport_data.h is in sync with the generator.
# Fails if the committed header differs from what the generator would produce.
# Run this manually or via pre-commit to catch accidental manual edits.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

TMPFILE=$(mktemp /tmp/thermo_transport_data_XXXXXX.h)
trap 'rm -f "$TMPFILE"' EXIT

.venv/bin/python thermo_data_generator/generate_thermo_data.py \
    --json thermo_data_generator/merged_species.json \
    --prefer-nasa9 \
    --species N2,O2,AR,CO2,H2O,CH4,C2H6,C3H8,IC4H10,NC5H12,NC6H14,NC7H16,CO,H2 \
    --output "$TMPFILE" \
    > /dev/null 2>&1

if ! diff -q "$TMPFILE" include/thermo_transport_data.h > /dev/null 2>&1; then
    echo "✗ include/thermo_transport_data.h is out of sync with the generator."
    echo "  Regenerate with:"
    echo "    python thermo_data_generator/generate_thermo_data.py \\"
    echo "        --json thermo_data_generator/merged_species.json --prefer-nasa9 \\"
    echo "        --species N2,O2,AR,CO2,H2O,CH4,C2H6,C3H8,IC4H10,NC5H12,NC6H14,NC7H16,CO,H2 \\"
    echo "        --output include/thermo_transport_data.h"
    echo ""
    diff "$TMPFILE" include/thermo_transport_data.h || true
    exit 1
fi

echo "✓ thermo_transport_data.h is up to date."
