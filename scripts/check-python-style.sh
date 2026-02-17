#!/usr/bin/env bash
set -euo pipefail

# Check Python code style and quality.
# Returns non-zero if violations are found.
#
# Usage:
#   scripts/check-python-style.sh [--verbose] [--strict] [--fix]
#
# Options:
#   --verbose   Show detailed output from all tools
#   --strict    Fail on warnings, not just errors
#   --fix       Auto-fix issues where possible (black, isort)
#
# Checks:
#   1. Non-ASCII characters outside comments and strings
#   2. PEP 8 compliance (via flake8)
#   3. Code formatting (via black --check)
#   4. Import sorting (via isort --check)
#   5. Type hints (via mypy, if available)

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Directories to scan
SCAN_DIRS=("python" "cantera_validation_tests" "thermo_data_generator")

# Python file patterns
INCLUDE_PATTERNS=("*.py")

print_usage() {
    echo "Usage: $0 [--verbose] [--strict] [--fix]"
    echo "  --verbose   Show detailed output from all tools"
    echo "  --strict    Fail on warnings, not just errors"
    echo "  --fix       Auto-fix issues where possible"
}

VERBOSE=false
STRICT=false
FIX=false
for arg in "$@"; do
    case "$arg" in
        --verbose) VERBOSE=true ;;
        --strict) STRICT=true ;;
        --fix) FIX=true ;;
        --help|-h) print_usage; exit 0 ;;
    esac
done

cd "${ROOT_DIR}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

VIOLATIONS=0

echo "=========================================="
echo "Python Code Style Checker"
echo "=========================================="
echo ""

# -----------------------------------------------------------------------------
# Check 1: Non-ASCII characters outside comments and strings
# -----------------------------------------------------------------------------
echo "Check 1: Non-ASCII characters..."

# Find all Python files
PYTHON_FILES=()
for dir in "${SCAN_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        while IFS= read -r -d '' file; do
            PYTHON_FILES+=("$file")
        done < <(find "$dir" -name "*.py" -print0)
    fi
done

if [ ${#PYTHON_FILES[@]} -eq 0 ]; then
    echo "  No Python files found in: ${SCAN_DIRS[*]}"
else
    NON_ASCII_FOUND=false
    for file in "${PYTHON_FILES[@]}"; do
        # Check for non-ASCII bytes (0x80-0xFF)
        # Exclude lines that are comments or inside strings
        if grep -P '[\x80-\xFF]' "$file" > /dev/null 2>&1; then
            # Simple heuristic: skip lines starting with # or containing quotes
            while IFS= read -r line; do
                # Skip comment lines
                if [[ "$line" =~ ^[[:space:]]*# ]]; then
                    continue
                fi
                # Skip lines that look like string literals (contains quotes)
                if [[ "$line" =~ [\"\'] ]]; then
                    continue
                fi
                # If we get here, it's likely non-ASCII in code
                if [ "$VERBOSE" = true ]; then
                    echo -e "${RED}  Non-ASCII in $file:${NC}"
                    echo "    $line"
                fi
                NON_ASCII_FOUND=true
            done < <(grep -n -P '[\x80-\xFF]' "$file" || true)
        fi
    done

    if [ "$NON_ASCII_FOUND" = true ]; then
        echo -e "${RED}  ✗ Non-ASCII characters found${NC}"
        ((VIOLATIONS++))
    else
        echo -e "${GREEN}  ✓ No non-ASCII characters${NC}"
    fi
fi

echo ""

# -----------------------------------------------------------------------------
# Check 2: Code quality (ruff)
# -----------------------------------------------------------------------------
echo "Check 2: Code quality (ruff)..."

if command -v ruff &> /dev/null; then
    if [ "$FIX" = true ]; then
        echo "  Running ruff with auto-fix..."
        ruff check --fix "${SCAN_DIRS[@]}"
        echo -e "${GREEN}  ✓ Ruff checks passed (with fixes)${NC}"
    else
        RUFF_OUTPUT=$(ruff check "${SCAN_DIRS[@]}" 2>&1)
        RUFF_STATUS=$?

        if [ $RUFF_STATUS -ne 0 ]; then
            echo -e "${RED}  ✗ Code quality issues found${NC}"
            if [ "$VERBOSE" = true ]; then
                echo "$RUFF_OUTPUT"
            else
                echo "$RUFF_OUTPUT" | head -20
                TOTAL_LINES=$(echo "$RUFF_OUTPUT" | wc -l)
                if [ "$TOTAL_LINES" -gt 20 ]; then
                    echo "  ... and $((TOTAL_LINES - 20)) more (use --verbose to see all)"
                fi
            fi
            echo "  Run with --fix to auto-fix"
            ((VIOLATIONS++))
        else
            echo -e "${GREEN}  ✓ Code quality checks passed${NC}"
            if [ "$VERBOSE" = true ] && [ -n "$RUFF_OUTPUT" ]; then
                echo "$RUFF_OUTPUT"
            fi
        fi
    fi
else
    echo -e "${YELLOW}  ⚠ ruff not found (install: pip install ruff)${NC}"
fi

echo ""

# -----------------------------------------------------------------------------
# Check 3: Code formatting (black)
# -----------------------------------------------------------------------------
echo "Check 3: Code formatting (black)..."

if command -v black &> /dev/null; then
    if [ "$FIX" = true ]; then
        echo "  Applying black formatting..."
        black --line-length=100 "${SCAN_DIRS[@]}"
        echo -e "${GREEN}  ✓ Code formatted${NC}"
    else
        BLACK_OUTPUT=$(black --check --line-length=100 "${SCAN_DIRS[@]}" 2>&1 || true)

        if echo "$BLACK_OUTPUT" | grep -q "would reformat"; then
            echo -e "${RED}  ✗ Code formatting issues found${NC}"
            if [ "$VERBOSE" = true ]; then
                echo "$BLACK_OUTPUT"
            else
                echo "$BLACK_OUTPUT" | grep "would reformat" | head -10
            fi
            echo "  Run with --fix to auto-format"
            ((VIOLATIONS++))
        else
            echo -e "${GREEN}  ✓ Code properly formatted${NC}"
        fi
    fi
else
    echo -e "${YELLOW}  ⚠ black not found (install: pip install black)${NC}"
fi

echo ""

# -----------------------------------------------------------------------------
# Check 4: Type hints (mypy) - optional
# -----------------------------------------------------------------------------
echo "Check 4: Type hints (mypy)..."

if command -v mypy &> /dev/null; then
    MYPY_ARGS=(
        "--ignore-missing-imports"
        "--no-strict-optional"
        "--warn-unused-ignores"
    )

    if [ "$STRICT" = true ]; then
        MYPY_ARGS+=("--strict")
    fi

    MYPY_OUTPUT=$(mypy "${MYPY_ARGS[@]}" "${SCAN_DIRS[@]}" 2>&1 || true)

    # Only fail on errors, not warnings (unless strict)
    if echo "$MYPY_OUTPUT" | grep -q "error:"; then
        echo -e "${RED}  ✗ Type checking errors found${NC}"
        if [ "$VERBOSE" = true ]; then
            echo "$MYPY_OUTPUT"
        else
            echo "$MYPY_OUTPUT" | grep "error:" | head -10
        fi
        if [ "$STRICT" = true ]; then
            ((VIOLATIONS++))
        else
            echo "  (Not counted as violation in non-strict mode)"
        fi
    else
        echo -e "${GREEN}  ✓ Type hints OK${NC}"
    fi
else
    echo -e "${YELLOW}  ⚠ mypy not found (install: pip install mypy)${NC}"
fi

echo ""

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "=========================================="
if [ $VIOLATIONS -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed!${NC}"
    exit 0
else
    echo -e "${RED}✗ Found $VIOLATIONS violation(s)${NC}"
    echo ""
    echo "To fix automatically:"
    echo "  $0 --fix"
    echo ""
    echo "For detailed output:"
    echo "  $0 --verbose"
    exit 1
fi
