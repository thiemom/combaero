---
description: Run style checks, fix issues, commit with message file, and push
---

This workflow automates the commit and push process with proper style checking and message handling.

## Steps

### 1. Stage all changes
```bash
git add -A
```

### 2. Run pre-commit style checks
```bash
pre-commit run --all-files
```

### 3. Re-stage any files modified by pre-commit hooks
```bash
git add -A
```

### 4. Fix any remaining issues that ruff couldn't auto-fix
If pre-commit reported errors that couldn't be auto-fixed, review the output and fix them manually. Common issues:
- Unused imports that need removal
- Type errors requiring annotation fixes
- Logic issues flagged by linters

Once fixed, re-run pre-commit and stage:
```bash
pre-commit run --all-files
git add -A
```

### 5. Create commit message file
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

### 6. Commit using message file
// turbo
```bash
git commit -F COMMIT_MSG.md
```

### 7. Push to remote
```bash
git push
```

### 8. Clean up message file
// turbo
```bash
rm COMMIT_MSG.md
```

## Usage

When you run this workflow, Cascade will:
1. Stage all changes
2. Run pre-commit hooks (ruff, trailing whitespace, etc.)
3. Re-stage any auto-fixed files
4. **Fix any remaining issues** that ruff couldn't auto-fix (agent steps in here)
5. Prompt you to create `COMMIT_MSG.md`
6. Commit using that file (avoiding shell escaping issues)
7. Push to remote
8. Clean up the temporary file

This avoids issues with special characters, backticks, and multi-line messages in shell commands.
