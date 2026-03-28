---
description: Auto-fix issues in changed Frontend files with Biome
---

# Biome Auto-Fix Changed Files

This workflow automatically fixes linting and formatting issues in all modified Frontend files using Biome.

## Steps

// turbo
1. Run Biome check with --write on all changed JS/TS/CSS/JSON files
```bash
git status --porcelain | grep '^[AMR].*gui/frontend/.*\.\(js\|ts\|tsx\|css\|json\)$' | cut -c4- | xargs -r npx @biomejs/biome check --write --no-errors-on-unmatched --files-max-size=1000000
```

// turbo
2. Verify all issues are fixed
```bash
git status --porcelain | grep '^[AMR].*gui/frontend/.*\.\(js\|ts\|tsx\|css\|json\)$' | cut -c4- | xargs -r npx @biomejs/biome check --no-errors-on-unmatched
```

## Usage

Run this workflow before committing to ensure all Frontend files pass Biome checks.

## Notes

- Only fixes files that are Added, Modified, or Renamed (AMR status) within the `gui/frontend` directory.
- Uses `check --write` to automatically apply formatting and safe lint fixes.
- The `-r` flag prevents xargs from failing when no files match.
- All steps are auto-turbo enabled for seamless execution.
- Always review the changes after auto-fix.
