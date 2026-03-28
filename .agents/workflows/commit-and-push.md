---
description: Run style checks, units sync validation, commit with message file, and push
---

// turbo-all
This workflow automates the commit and push process with comprehensive checks and message handling.

## Steps

### 1. Check for old commit message files
```bash
./scripts/check-commit-msg.sh
```
Note: If a `COMMIT_MSG.md` exists, review it and delete if stale, or keep to reuse it.

### 2. Run Python formatting and linting
```bash
.venv/bin/ruff format python/
.venv/bin/ruff check --fix python/
```

### 3. Run Frontend formatting and linting (Biome)
```bash
cd gui/frontend && npm run format && cd ../..
```

### 4. Check for units synchronization (Mandatory)
```bash
.venv/bin/python -m pytest python/tests/test_units_sync.py
```

### 5. Check for non-ASCII characters in C++ source
```bash
./scripts/check-source-style.sh
```

### 6. Stage and run pre-commit hooks
```bash
git add -A
pre-commit run --all-files
git add -A
```

### 7. Create commit message file
Create a file named `COMMIT_MSG.md` in the project root with your commit message.
Follow conventional commit format: `type: brief description`.

### 8. Commit and cleanup
```bash
git commit -F COMMIT_MSG.md
rm COMMIT_MSG.md
```

### 9. Push to current branch
```bash
git push
```
