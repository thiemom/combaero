---
description: Squash merge PR and cleanup branch (default for small tasks)
---

// turbo-all
This workflow handles the final merge and cleanup. It uses squash merge by default as preferred for small task branches.

## Merge Logic Selection
- **Squash Merge (`--squash`)**: Default for small task branches, bugfixes, and feature refinements to keep history clean.
- **Normal Merge (`--merge`)**: Used for major architectural changes or multi-contributor feature branches where preserving intermediate commit history is vital.

## Steps

### 1. Verification Checklist
- [ ] All CI checks are green (confirm with `gh pr view`).
- [ ] Any manual verification/benchmarks are completed.

### 2. Execute Squash Merge (Preferred)
```bash
gh pr merge --squash --delete-branch
```

### 3. Alternative: Normal Merge (Major Work)
If the work is major according to project conventions, use:
```bash
gh pr merge --merge --delete-branch
```

### 4. Cleanup Local Branch
```bash
git checkout main
git pull
# Local branch deletion is handled by --delete-branch if remote, but check local:
# git branch -d <branch_name>
```
