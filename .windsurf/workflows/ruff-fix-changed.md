---
description: Auto-fix linting issues in changed Python files with ruff
---

# Ruff Auto-Fix Changed Files

This workflow automatically fixes linting issues in all modified Python files using ruff.

## Steps

// turbo
1. Run ruff check with --fix on all changed Python files
```bash
git status --porcelain | grep '^[AMR].*\.py$' | cut -c4- | xargs -r uv run ruff check --fix
```

// turbo
2. Run ruff format on all changed Python files
```bash
git status --porcelain | grep '^[AMR].*\.py$' | cut -c4- | xargs -r uv run ruff format
```

// turbo
3. Verify all issues are fixed
```bash
git status --porcelain | grep '^[AMR].*\.py$' | cut -c4- | xargs -r uv run ruff check
```

## Usage

Run this workflow before committing to ensure all Python files pass linting checks.

## Notes

- Only fixes files that are Added, Modified, or Renamed (AMR status)
- Uses `--fix` flag to automatically apply safe fixes
- Also runs `ruff format` to ensure consistent formatting
- The `-r` flag prevents xargs from failing when no files match
- All steps are auto-turbo enabled for seamless execution
- Always review the changes after auto-fix
