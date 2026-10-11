"""Check that tool versions pinned in more than one place agree.

Dependabot bumps one file (e.g. gui/frontend/package.json) and leaves the
others, so a green PR can still drift the pins apart. Run in CI, in
pre-commit and by scripts/bulk_merge_safe.sh before it merges.

    uv run scripts/check_version_sync.py            # working tree
    uv run scripts/check_version_sync.py --ref REF  # any git ref, e.g. a PR head
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def reader(ref: str | None) -> Callable[[str], str]:
    if ref is None:
        return lambda path: (ROOT / path).read_text()
    return lambda path: (
        subprocess.run(
            ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout
    )


def biome_pins(read: Callable[[str], str]) -> dict[str, str | None]:
    """Every place the Biome version is pinned, by file."""

    def find(pattern: str, text: str) -> str | None:
        m = re.search(pattern, text)
        return m.group(1) if m else None

    pkg = json.loads(read("gui/frontend/package.json"))
    deps = {**pkg.get("dependencies", {}), **pkg.get("devDependencies", {})}
    return {
        "gui/frontend/package.json": deps.get("@biomejs/biome", "").lstrip("^~") or None,
        "biome.json ($schema)": find(
            r"biomejs\.dev/schemas/([0-9.]+)/schema\.json", read("biome.json")
        ),
        ".pre-commit-config.yaml": find(
            r"@biomejs/biome@([0-9.]+)", read(".pre-commit-config.yaml")
        ),
        ".github/workflows/ci.yml (setup-biome)": find(
            r"setup-biome@[^\n]*\n\s*with:\s*\n\s*version:\s*['\"]?([0-9.]+)",
            read(".github/workflows/ci.yml"),
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", help="check this git ref instead of the working tree")
    args = ap.parse_args()
    pins = biome_pins(reader(args.ref))
    versions = set(pins.values())
    if len(versions) == 1 and None not in versions:
        print(f"version sync OK: biome {versions.pop()}")
        return 0
    print("version sync FAILED: biome is pinned differently (Version Sync Checklist in CLAUDE.md):")
    for where, v in pins.items():
        print(f"  {where:42s} {v or 'NOT FOUND'}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
