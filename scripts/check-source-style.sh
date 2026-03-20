#!/usr/bin/env bash
set -euo pipefail

# Check for source code style violations.
# Returns non-zero if violations are found.
#
# Usage:
#   scripts/check-source-style.sh [--verbose] [--strict]
#
# Options:
#   --verbose   Show all matches with line content
#   --strict    Also flag non-ASCII in string literals (default: allow in strings)
#
# Checks:
#   1. Non-ASCII characters (0x80-0xFF) outside comments and strings
#   2. Block comments (/* ... */) - only line comments (//) are allowed
#   3. M_PI usage without math_constants.h (MSVC compatibility)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Directories to scan (space-separated for grep -R)
SCAN_DIRS=("src" "include" "examples" "tests")

# File extensions to check
INCLUDE_PATTERNS=("*.cpp" "*.h" "*.hpp" "*.c")

print_usage() {
    echo "Usage: $0 [--verbose] [--strict]"
    echo "  --verbose   Show all matches with line content"
    echo "  --strict    Also flag non-ASCII in string literals"
}

VERBOSE=false
STRICT=false
for arg in "$@"; do
    case "$arg" in
        --verbose) VERBOSE=true ;;
        --strict) STRICT=true ;;
        --help|-h) print_usage; exit 0 ;;
    esac
done

cd "${ROOT_DIR}"

# Build include arguments for grep
INCLUDE_ARGS=()
for pat in "${INCLUDE_PATTERNS[@]}"; do
    INCLUDE_ARGS+=("--include=$pat")
done

# Temporary file for results
TMPFILE=$(mktemp)
trap 'rm -f "$TMPFILE"' EXIT

# Use grep to find non-ASCII characters (bytes with high bit set: 0x80-0xFF)
# LC_ALL=C ensures byte-level matching
# The pattern [^[:print:][:space:]] catches non-printable, but we want non-ASCII specifically
# Using octal range: [\200-\377] matches bytes 0x80-0xFF
# Exclude CMakeFiles directories (CMake-generated code)
if LC_ALL=C grep -r -n "${INCLUDE_ARGS[@]}" --exclude-dir=CMakeFiles $'[\x80-\xff]' "${SCAN_DIRS[@]}" > "$TMPFILE" 2>/dev/null; then
    HAS_MATCHES=true
else
    HAS_MATCHES=false
fi

# Filter out matches that are inside comments (and optionally strings)
# This is a simplified check - it won't catch all edge cases but handles common patterns
VIOLATIONS=()
while IFS= read -r line; do
    # Extract filename and line number
    file_line="${line%%:*}"
    rest="${line#*:}"
    lineno="${rest%%:*}"
    content="${rest#*:}"

    # Skip lines that are clearly inside block comments (start with * or /**)
    trimmed="${content#"${content%%[![:space:]]*}"}"  # trim leading whitespace
    if [[ "$trimmed" == \** ]] || [[ "$trimmed" == /\*\** ]]; then
        continue  # Likely a block comment line
    fi

    # Skip if the non-ASCII char appears only in a // comment
    # Extract code before any // comment
    code_part="${content%%//*}"

    # Check if non-ASCII exists in the code part (before //)
    if ! echo "$code_part" | LC_ALL=C grep -q $'[\x80-\xff]' 2>/dev/null; then
        continue  # Only in comment, skip
    fi

    # Unless --strict, allow non-ASCII inside string literals
    if [[ "$STRICT" != "true" ]]; then
        # Remove content inside double quotes (simple heuristic)
        # This won't handle escaped quotes perfectly but covers most cases
        stripped=$(echo "$code_part" | sed 's/"[^"]*"//g')
        if ! echo "$stripped" | LC_ALL=C grep -q $'[\x80-\xff]' 2>/dev/null; then
            continue  # Only in strings, skip
        fi
    fi

    VIOLATIONS+=("$line")
done < "$TMPFILE"

# Check for block comments (/* ... */)
BLOCK_COMMENT_FILES=()
BLOCK_TMPFILE=$(mktemp)
if grep -r -l "${INCLUDE_ARGS[@]}" --exclude-dir=CMakeFiles '/\*' "${SCAN_DIRS[@]}" > "$BLOCK_TMPFILE" 2>/dev/null; then
    while IFS= read -r file; do
        BLOCK_COMMENT_FILES+=("$file")
    done < "$BLOCK_TMPFILE"
fi
rm -f "$BLOCK_TMPFILE"

# Summary
HAS_ERRORS=false

if [[ ${#VIOLATIONS[@]} -gt 0 ]]; then
    HAS_ERRORS=true
    echo ""
    echo "=========================================="
    echo "ERROR: Non-ASCII characters found in code"
    echo "=========================================="
    echo ""

    if [[ "$VERBOSE" == "true" ]]; then
        echo "Violations:"
        for v in "${VIOLATIONS[@]}"; do
            echo "  $v"
        done
        echo ""
    fi

    # Extract unique files
    declare -A seen_files
    for v in "${VIOLATIONS[@]}"; do
        f="${v%%:*}"
        seen_files["$f"]=1
    done

    echo "Files with violations (${#VIOLATIONS[@]} total):"
    for f in "${!seen_files[@]}"; do
        echo "  - $f"
    done
    echo ""

    if [[ "$VERBOSE" != "true" ]]; then
        echo "Run with --verbose to see details."
        echo ""
    fi

    echo "Manual check command:"
    echo "  LC_ALL=C grep -r -n --include='*.cpp' --include='*.h' \$'[\\x80-\\xff]' src/ include/"
    echo ""
fi

if [[ ${#BLOCK_COMMENT_FILES[@]} -gt 0 ]]; then
    HAS_ERRORS=true
    echo ""
    echo "=========================================="
    echo "ERROR: Block comments (/* */) found"
    echo "=========================================="
    echo "Only line comments (//) are allowed."
    echo ""
    echo "Files with block comments:"
    for f in "${BLOCK_COMMENT_FILES[@]}"; do
        echo "  - $f"
    done
    echo ""
    echo "Manual check command:"
    echo "  grep -r -n --include='*.cpp' --include='*.h' '/\\*' src/ include/"
    echo ""
fi

# Check for POSIX math constants used without math_constants.h include.
# These are undefined on MSVC by default: M_PI, M_PI_2, M_PI_4, M_SQRT2, M_LN10.
# The list is derived directly from include/math_constants.h — keep in sync.
MATH_CONST_PATTERN='M_PI\b|M_PI_2\b|M_PI_4\b|M_SQRT2\b|M_LN10\b'
MATH_CONST_VIOLATIONS=()
for dir in "${SCAN_DIRS[@]}"; do
    if [[ -d "$dir" ]]; then
        while IFS= read -r file; do
            # Skip math_constants.h itself
            [[ "$file" == *"math_constants.h" ]] && continue
            # Check if file uses any of the POSIX math constants
            if grep -qE "$MATH_CONST_PATTERN" "$file" 2>/dev/null; then
                # Check if it includes math_constants.h
                if ! grep -q 'math_constants\.h' "$file" 2>/dev/null; then
                    MATH_CONST_VIOLATIONS+=("$file")
                fi
            fi
        done < <(find "$dir" -type f \( -name "*.cpp" -o -name "*.h" -o -name "*.hpp" \) 2>/dev/null)
    fi
done

if [[ ${#MATH_CONST_VIOLATIONS[@]} -gt 0 ]]; then
    HAS_ERRORS=true
    echo ""
    echo "=========================================="
    echo "ERROR: POSIX math constant used without math_constants.h"
    echo "=========================================="
    echo "M_PI, M_PI_2, M_PI_4, M_SQRT2, M_LN10 are POSIX extensions"
    echo "not available on MSVC. Include math_constants.h for cross-platform compatibility."
    echo ""
    echo "Files missing include:"
    for f in "${MATH_CONST_VIOLATIONS[@]}"; do
        # Show which constant(s) triggered the violation
        consts=$(grep -oE "$MATH_CONST_PATTERN" "$f" | sort -u | tr '\n' ' ')
        echo "  - $f  (uses: ${consts% })"
    done
    echo ""
    echo "Fix: Add this include to each file:"
    echo '  #include "math_constants.h"  // MSVC compatibility for M_PI, M_LN10, etc.'
    echo ""
fi

if [[ "$HAS_ERRORS" == "true" ]]; then
    exit 1
else
    if [[ "$STRICT" == "true" ]]; then
        echo "✓ All checks passed (strict mode)."
    else
        echo "✓ All checks passed (non-ASCII allowed in comments and strings)."
    fi
    exit 0
fi
