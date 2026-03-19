#!/usr/bin/env bash
# Check for existing COMMIT_MSG.md file before starting commit workflow.
# This prevents accidentally committing stale message files.
#
# Usage:
#   scripts/check-commit-msg.sh
#
# Exit codes:
#   0 - No COMMIT_MSG.md found (safe to proceed)
#   1 - COMMIT_MSG.md exists (user should review/delete)

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

if [ -f COMMIT_MSG.md ]; then
  echo "⚠️  Found existing COMMIT_MSG.md"
  echo "This may be from a previous commit attempt."
  echo "Please review and delete if stale, or keep if you want to reuse it."
  exit 1
fi

exit 0
