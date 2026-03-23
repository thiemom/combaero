---
description: Create a Pull Request and monitor CI status
---

// turbo-all
This workflow automates the creation of a Pull Request using the GitHub CLI and tracks CI progress.

## Steps

### 1. Create Pull Request
```bash
gh pr create --title "type: brief description" --body-file COMMIT_MSG.md
```
Note: Re-use the `COMMIT_MSG.md` if available, or specify a summary.

### 2. Monitor CI Status
```bash
gh pr view --json statusCheckRollup --jq '.statusCheckRollup[]'
```

### 3. Check specific run logs if needed
```bash
gh run list --limit 3
# Then view specific run
# gh run view <ID> --json jobs --jq '.jobs[] | {name, status, conclusion}'
```
