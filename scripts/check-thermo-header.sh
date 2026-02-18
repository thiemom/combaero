#!/usr/bin/env bash
# Check that include/thermo_transport_data.h is in sync with the generator.
# Fails if the committed header differs from what the generator would produce.
# Run this manually or via pre-commit to catch accidental manual edits.
#
# Species list and source are read from thermo_data_generator/.generator_state.json,
# which is written automatically by generate_thermo_data.py on every run.
# No manual updates needed when species change — just regenerate and commit both files.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

SIDECAR="thermo_data_generator/.generator_state.json"
if [[ ! -f "$SIDECAR" ]]; then
    echo "✗ $SIDECAR not found."
    echo "  Run the generator first:"
    echo "    python thermo_data_generator/generate_thermo_data.py \\"
    echo "        --json thermo_data_generator/merged_species.json --prefer-nasa9 \\"
    echo "        --output include/thermo_transport_data.h"
    exit 1
fi

SPECIES=$(.venv/bin/python3 -c "import json; d=json.load(open('$SIDECAR')); print(','.join(d['species']))")
JSON_SOURCE=$(.venv/bin/python3 -c "import json; d=json.load(open('$SIDECAR')); print(d['json_source'])")

TMPFILE=$(mktemp /tmp/thermo_transport_data_XXXXXX.h)
trap 'rm -f "$TMPFILE"' EXIT

.venv/bin/python thermo_data_generator/generate_thermo_data.py \
    --json "$JSON_SOURCE" \
    --prefer-nasa9 \
    --species "$SPECIES" \
    --output "$TMPFILE" \
    > /dev/null 2>&1

if ! diff -q "$TMPFILE" include/thermo_transport_data.h > /dev/null 2>&1; then
    echo "✗ include/thermo_transport_data.h is out of sync with the generator."
    echo "  Regenerate with:"
    echo "    python thermo_data_generator/generate_thermo_data.py \\"
    echo "        --json $JSON_SOURCE --prefer-nasa9 \\"
    echo "        --species $SPECIES \\"
    echo "        --output include/thermo_transport_data.h"
    diff "$TMPFILE" include/thermo_transport_data.h || true
    exit 1
fi

echo "✓ thermo_transport_data.h is up to date."
