---
description: Create PR and monitor GitHub Actions checks until complete
---

# Create PR and Monitor Checks

This workflow creates a pull request and monitors GitHub Actions checks until they complete.
Does NOT include merging or branch cleanup - those are separate steps.

## Prerequisites
- Changes are committed on a feature branch
- `gh` CLI is installed and authenticated
- Branch is pushed to origin

## Steps

// turbo
### 1. Ensure branch is pushed
```bash
git push -u origin $(git branch --show-current)
```

// turbo
### 2. Create pull request

The AI will generate an appropriate PR title and body based on recent commits and changes.

**For simple PRs (single-line description):**
```bash
gh pr create --title "feat: brief description" --body "Simple one-line summary."
```

**For complex PRs (multi-line description - RECOMMENDED):**

Create a PR body file programmatically:
```bash
PR_BODY_FILE="/tmp/pr_body_$(date +%s).md"
printf "## Summary\n\nBrief description here\n\n## Changes\n- Change 1\n- Change 2\n\n## Testing\n- Test results\n" > "$PR_BODY_FILE"
gh pr create --title "TITLE_HERE" --body-file "$PR_BODY_FILE"
rm "$PR_BODY_FILE"
```

**Notes:**
- AI-generated titles and descriptions are consistently good
- Use unique filenames (`$(date +%s)`) to avoid conflicts
- The AI will analyze commits to generate appropriate PR content

// turbo
### 3. Monitor checks until complete
```bash
gh pr checks --watch
```

This will:
- Display all GitHub Actions checks
- Auto-refresh every few seconds
- Exit when all checks complete (pass or fail)
- Show final status

### 4. Review check results
If any checks failed, view detailed logs:
```bash
gh pr checks
```

Or view specific run logs:
```bash
gh run view <run-id> --log-failed
```

## Example Usage

```bash
# Push current branch
git push -u origin feat/my-feature

# Create PR with inline body
gh pr create --title "feat: add new feature" --body "Implements X, Y, Z. All tests pass."

# Monitor until complete
gh pr checks --watch

# If checks pass, proceed to merge (separate workflow)
# If checks fail, review logs and fix issues
```

## Notes
- The `--watch` flag blocks until checks complete
- You can Ctrl+C to stop watching without affecting the checks
- Checks continue running on GitHub even if you stop watching
- Use `gh pr view` to see PR details and status
- Merging and branch cleanup are intentionally separate workflows
