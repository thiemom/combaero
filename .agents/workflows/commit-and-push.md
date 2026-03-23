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

### 3. Check for units synchronization (Mandatory)
```bash
.venv/bin/python -m pytest python/tests/test_units_sync.py
```

### 4. Check for non-ASCII characters in C++ source
```bash
./scripts/check-source-style.sh
```

### 5. Stage and run pre-commit hooks
```bash
git add -A
pre-commit run --all-files
git add -A
```

### 6. Create commit message file
Create a file named `COMMIT_MSG.md` in the project root with your commit message.
Follow conventional commit format: `type: brief description`.

### 7. Commit and cleanup
```bash
git commit -F COMMIT_MSG.md
rm COMMIT_MSG.md
```

### 8. Push to current branch
```bash
git push
```
