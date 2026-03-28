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
### 2. Run ruff format on Python files (faster than pre-commit)
```bash
.venv/bin/ruff format python/
```

// turbo
### 3. Run ruff check with auto-fix on Python files
```bash
.venv/bin/ruff check --fix python/
```

// turbo
### 4. Run Frontend formatting and linting (Biome)
```bash
cd gui/frontend && npm run format && cd ../..
```

// turbo
### 5. Check for non-ASCII characters in C++ source
```bash
./scripts/check-source-style.sh
```
Note: Ignore errors from `tests/CMakeFiles/` (CMake-generated files). Only fix issues in `src/`, `include/`, `examples/`.

// turbo
### 6. Stage all changes
```bash
git add -A
```

// turbo
### 7. Run pre-commit hooks (should be fast now, most issues already fixed)
```bash
pre-commit run --all-files
```

// turbo
### 8. Re-stage any files modified by pre-commit hooks
```bash
git add -A
```

### 9. Fix any remaining issues that ruff/biome couldn't auto-fix
If pre-commit reported errors that couldn't be auto-fixed, review the output and fix them manually. Common issues:
- Unused imports that need removal
- Type errors requiring annotation fixes
- Logic issues flagged by linters

Once fixed, re-run pre-commit and stage:
```bash
pre-commit run --all-files
git add -A
```

### 10. Create commit message file
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

### 11. Commit using message file
// turbo
```bash
git commit -F COMMIT_MSG.md
```

### 12. Clean up message file
// turbo
```bash
rm COMMIT_MSG.md
```

## Usage

When you run this workflow, Cascade will:
1. Check for old commit message files
2. Run ruff format on Python files (faster than pre-commit)
3. Run ruff check --fix on Python files
4. Run Biome format/check on Frontend (JS/TS)
5. Check for non-ASCII characters in C++ source
6. Stage all changes
7. Run pre-commit hooks (should be fast now, most issues already fixed)
8. Re-stage any files modified by pre-commit hooks
9. Fix any remaining issues that ruff/biome couldn't auto-fix
10. Create commit message file
11. Commit using that file (avoiding shell escaping issues)
12. Clean up the temporary file

**Why run ruff/biome manually before pre-commit?**
Pre-commit runs the full test suite which is slow (~60s). Running ruff format/check, Biome format/check, and ASCII validation first catches 95% of issues in ~1s, making the final pre-commit run much faster since tests only run if hooks pass.
