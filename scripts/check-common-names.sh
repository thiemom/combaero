#!/usr/bin/env bash
# check-common-names.sh
# Verify that every species in include/thermo_transport_data.h has an entry
# in include/common_names.h.
#
# Exits 0 if all active species have a common name defined.
# Exits 1 and prints the missing species if any are absent.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"

python3 - "${REPO_ROOT}" <<'EOF'
import re
import sys
from pathlib import Path

repo_root = Path(sys.argv[1])
thermo_h = repo_root / "include" / "thermo_transport_data.h"
names_h  = repo_root / "include" / "common_names.h"

if not thermo_h.exists():
    print(f"check-common-names: {thermo_h} not found, skipping.")
    sys.exit(0)

if not names_h.exists():
    print(f"check-common-names: {names_h} not found.")
    sys.exit(1)

# Extract active species from the species_names vector.
# The vector is a single line: const std::vector<std::string> species_names = {"N2", ...};
thermo_text = thermo_h.read_text()
m = re.search(r'species_names\s*=\s*\{([^;]+)\}', thermo_text)
if not m:
    print("check-common-names: could not find species_names in thermo_transport_data.h")
    sys.exit(1)
active = set(re.findall(r'"([^"]+)"', m.group(1)))

# Extract known formulas from formula_to_name map in common_names.h.
# Each entry is on its own line:  {"FORMULA", "Common name"},
# The formula key is always the first quoted string on lines inside the map.
names_text = names_h.read_text()
# Find the formula_to_name block: from the opening { to the matching };
m2 = re.search(r'formula_to_name\{(.*?)\};', names_text, re.DOTALL)
known = set(re.findall(r'\{\s*"([^"]+)"\s*,', m2.group(1))) if m2 else set()

missing = sorted(active - known)
if missing:
    print(f"check-common-names: ERROR - {len(missing)} active species have no entry in common_names.h:")
    for s in missing:
        print(f"  {s}")
    print()
    print("  Add them to include/common_names.h (both formula_to_name and name_to_formula maps).")
    sys.exit(1)

print(f"check-common-names: all {len(active)} active species have common names defined.")
sys.exit(0)
EOF
