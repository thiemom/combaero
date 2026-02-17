#!/usr/bin/env python3
"""Fail fast if the active interpreter is not the repo-local .venv Python."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only print output when validation fails.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    expected_python = (repo_root / ".venv" / "bin" / "python").resolve()
    actual_python = Path(sys.executable).resolve()

    if actual_python != expected_python:
        if not args.quiet:
            print("ERROR: wrong Python interpreter for this repository.")
            print(f"  expected: {expected_python}")
            print(f"  actual:   {actual_python}")
            print("Fix:")
            print("  1) ./scripts/bootstrap.sh")
            print("  2) source .venv/bin/activate")
        return 1

    if not args.quiet:
        print(f"Using project venv interpreter: {actual_python}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
