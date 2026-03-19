---
description: Run style checks, fix issues, commit with message file, and push
---

This workflow automates the commit and push process with proper style checking and message handling.

## Steps

### 1-8. Complete commit workflow
Follow all steps from @[/commit-only] workflow (check for old commit messages, stage, pre-commit, fix issues, create message file, commit, cleanup).

// turbo
### 9. Push to remote
```bash
git push
```

## Usage

This workflow runs the complete @[/commit-only] workflow and then pushes to the remote repository.

When you run this workflow, Cascade will:
1. Run all commit-only steps (check, stage, pre-commit, fix, create message, commit, cleanup)
2. Push to remote

This avoids issues with special characters, backticks, and multi-line messages in shell commands.
