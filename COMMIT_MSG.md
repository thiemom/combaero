fix: Replace non-ASCII characters and fix macOS pre-commit portability

## Non-ASCII Character Fixes
- Replace kg/m³ → kg/m^3
- Replace J/(kg·K) → J/(kg*K)
- Replace Σ → sum()
- Fixes CI lint failure on GitHub Actions

## Pre-commit Script Portability Fix
- Replace `grep -P` (GNU grep only) with `LC_ALL=C grep` (works on both macOS BSD grep and Linux GNU grep)
- Root cause: macOS ships with BSD grep which doesn't support `-P` flag
- Previous behavior: Non-ASCII check silently failed on macOS but passed on Linux CI
- New behavior: Works correctly on both platforms

## Additional Improvements
- Add math constants portability pre-commit hook (checks for M_PI usage without math_constants.h)
- Update commit-and-push workflow with manual fix step for issues ruff can't auto-fix
