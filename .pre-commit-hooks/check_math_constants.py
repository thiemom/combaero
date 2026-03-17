#!/usr/bin/env python3
"""
Pre-commit hook to ensure M_PI and related constants are used via math_constants.h

MSVC does not define M_PI by default (POSIX extension), so we must use our
portable math_constants.h header for cross-platform compatibility.
"""

import re
import sys
from pathlib import Path

# Constants that require math_constants.h
MATH_CONSTANTS = [
    "M_PI",
    "M_PI_2",
    "M_PI_4",
    "M_SQRT2",
    "M_LN10",
]

def check_file(filepath: Path) -> list[str]:
    """Check if file uses math constants without including math_constants.h"""
    errors = []

    try:
        content = filepath.read_text(encoding='utf-8')
    except (UnicodeDecodeError, FileNotFoundError):
        return errors

    # Check if any math constants are used
    constants_used = set()
    for const in MATH_CONSTANTS:
        # Match constant but not in comments or strings
        if re.search(rf'\b{const}\b', content):
            constants_used.add(const)

    if not constants_used:
        return errors

    # Check if math_constants.h is included (matches both "math_constants.h" and "../include/math_constants.h")
    has_include = bool(re.search(r'#include\s+[<"].*math_constants\.h[>"]', content))

    if not has_include:
        errors.append(
            f"{filepath}: Uses {', '.join(sorted(constants_used))} "
            f"but missing #include \"math_constants.h\" (required for MSVC compatibility)"
        )

    return errors

def main():
    """Check all C++ files passed as arguments"""
    errors = []

    for arg in sys.argv[1:]:
        filepath = Path(arg)

        # Only check C++ source/header files
        if filepath.suffix not in {'.cpp', '.h', '.hpp', '.cc', '.cxx'}:
            continue

        # Skip the math_constants.h file itself
        if filepath.name == 'math_constants.h':
            continue

        errors.extend(check_file(filepath))

    if errors:
        print("\n❌ Math constants portability check failed:\n")
        for error in errors:
            print(f"  {error}")
        print("\nMSVC does not define M_PI by default. Use #include \"math_constants.h\" for portability.")
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())
