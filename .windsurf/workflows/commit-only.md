---
description: Run style checks, fix issues, commit with message file, but NO push
---

This workflow automates the commit process with proper style checking and message handling.

## Steps

### 1. Check for old commit message files
```bash
./scripts/check-commit-msg.sh
```
Note: If a COMMIT_MSG.md exists, review it and delete if stale, or keep to reuse it.
// turbo
### 2. Run Python style checks with auto-fix (ruff format + ruff check)
```bash
./scripts/check-python-style.sh --fix
```

// turbo
### 3. Run GUI style checks with auto-fix (Biome)
```bash
./scripts/check-gui-style.sh
```

// turbo
### 4. Check for non-ASCII characters in C++ source
```bash
./scripts/check-source-style.sh
```
Note: Ignore errors from `tests/CMakeFiles/` (CMake-generated files). Only fix issues in `src/`, `include/`, `examples/`.

// turbo
### 5. Stage all changes
```bash
git add -A
```

// turbo
### 6. Run pre-commit hooks (should be fast now, most issues already fixed)
```bash
pre-commit run --all-files
```

// turbo
### 7. Re-stage any files modified by pre-commit hooks
```bash
git add -A
```

### 8. Fix any remaining issues that ruff/biome couldn't auto-fix
If pre-commit reported errors that couldn't be auto-fixed, review the output and fix them manually. Common issues:
- Unused imports that need removal
- Type errors requiring annotation fixes
- Logic issues flagged by linters

Once fixed, re-run pre-commit and stage:
```bash
pre-commit run --all-files
git add -A
```

### 9. Create commit message file
Create a file named `COMMIT_MSG.md` in the project root with your commit message.

The file should follow conventional commit format:
```
type: brief description

Detailed explanation if needed.
Can be multiple paragraphs.

- Bullet points work
- List changes clearly
```

Common types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `style`

### 10. Commit using message file
// turbo
```bash
git commit -F COMMIT_MSG.md
```

### 11. Clean up message file
// turbo
```bash
rm COMMIT_MSG.md
```

## Usage

When you run this workflow, Cascade will:
1. Check for old commit message files
2. Run Python style checks with auto-fix (ruff format + ruff check via check-python-style.sh)
3. Run GUI style checks with auto-fix (Biome via check-gui-style.sh)
4. Check for non-ASCII characters in C++ source
5. Stage all changes
6. Run pre-commit hooks (should be fast now, most issues already fixed)
7. Re-stage any files modified by pre-commit hooks
8. Fix any remaining issues that scripts couldn't auto-fix
9. Create commit message file
10. Commit using that file (avoiding shell escaping issues)
11. Clean up the temporary file

**Why run style scripts before pre-commit?**
Pre-commit runs the full test suite which is slow (~60s). Running the style scripts first catches 95% of issues in ~1s, making the final pre-commit run much faster since tests only run if hooks pass.
