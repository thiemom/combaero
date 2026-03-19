---
description: Run style checks, fix issues, commit with message file, and push
---

This workflow automates the commit and push process with proper style checking and message handling.

## Steps

// turbo
### 1. Check for old commit message files
```bash
if [ -f COMMIT_MSG.md ]; then
  echo "⚠️  Found existing COMMIT_MSG.md"
  echo "This may be from a previous commit attempt."
  echo "Please review and delete if stale, or keep if you want to reuse it."
  exit 1
fi
```
Note: keeping this file will commit it. Do not leave it.
If a previous COMMIT_MSG.md exists in the repository: review it and delete if stale.

// turbo
### 2. Stage all changes
```bash
git add -A
```

// turbo
### 3. Run pre-commit style checks
```bash
pre-commit run --all-files
```

// turbo
### 4. Re-stage any files modified by pre-commit hooks
```bash
git add -A
```

### 5. Fix any remaining issues that ruff couldn't auto-fix
If pre-commit reported errors that couldn't be auto-fixed, review the output and fix them manually. Common issues:
- Unused imports that need removal
- Type errors requiring annotation fixes
- Logic issues flagged by linters

Once fixed, re-run pre-commit and stage:
```bash
pre-commit run --all-files
git add -A
```

### 6. Create commit message file
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

### 7. Commit using message file
// turbo
```bash
git commit -F COMMIT_MSG.md
```

// turbo
### 8. Push to remote
```bash
git push
```

### 9. Clean up message file
// turbo
```bash
rm COMMIT_MSG.md
```

## Usage

When you run this workflow, Cascade will:
When you run this workflow, Cascade will:
1. Check for old commit message files
2. Stage all changes
3. Run pre-commit hooks (ruff, trailing whitespace, etc.)
4. Re-stage any files modified by pre-commit hooks
5. Fix any remaining issues that ruff couldn't auto-fix
6. Create commit message file
7. Commit using that file (avoiding shell escaping issues)
8. Push to remote
9. Clean up the temporary file

This avoids issues with special characters, backticks, and multi-line messages in shell commands.
