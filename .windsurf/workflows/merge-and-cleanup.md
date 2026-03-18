---
description: Squash merge PR and cleanup branches
---

# Merge PR and Cleanup Branches

This workflow squash merges a pull request and cleans up both local and remote branches.

## Prerequisites
- PR exists and all checks have passed
- You have the PR number ready

## Steps

// turbo
### 1. Squash merge PR and delete remote branch
```bash
gh pr merge <PR_NUMBER> --squash --delete-branch
```

Replace `<PR_NUMBER>` with the actual PR number (e.g., `51`).

This will:
- Squash all commits into one
- Merge into the base branch (usually `main`)
- Delete the remote feature branch
- Attempt to delete the local branch and switch to main

// turbo
### 2. Pull latest changes from main
```bash
git pull --rebase origin main
```

This ensures your local main branch is up-to-date with the squashed commit.

// turbo
### 3. Verify branch cleanup
```bash
git branch -a | grep -v "main\|master"
```

This shows any remaining branches. If the feature branch still exists locally, delete it manually:
```bash
git branch -D <branch-name>
```

## Usage

When you run this workflow:
1. Provide the PR number when prompted
2. The PR will be squash merged
3. Remote branch will be deleted
4. Local main will be updated
5. Local feature branch will be cleaned up

## Notes

- The `--squash` flag combines all commits into a single commit
- The `--delete-branch` flag removes the remote branch after merge
- If local branch deletion fails, step 3 helps identify and clean up manually
- All steps are auto-turbo enabled for seamless execution
