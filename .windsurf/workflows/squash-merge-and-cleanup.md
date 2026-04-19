# Squash Merge and Branch Cleanup

This workflow performs squash merge with branch cleanup after all checks have passed.

## Prerequisites
- Pull request is created and all checks are passing
- No merge conflicts detected
- Sufficient permissions to merge the PR

## Steps

### 1. Check PR status and conflicts
```bash
# Get PR info in single call for efficiency
PR_INFO=$(gh pr view --json number,mergeStateStatus,state,mergeable)
PR_NUMBER=$(echo "$PR_INFO" | jq -r '.number')
CONFLICT_STATUS=$(echo "$PR_INFO" | jq -r '.mergeStateStatus')
CHECK_STATUS=$(echo "$PR_INFO" | jq -r '.state + " " + .mergeable')

echo "Checking PR #$PR_NUMBER..."
echo "PR status: $CHECK_STATUS"

# Early exit for conflicts
if [ "$CONFLICT_STATUS" = "DIRTY" ]; then
    echo "❌ Merge conflicts detected in PR #$PR_NUMBER"
    echo "Please resolve conflicts manually before merging."
    exit 1
fi

# Early exit if not mergeable
if echo "$CHECK_STATUS" | grep -q "NOT_MERGEABLE"; then
    echo "❌ PR is not mergeable"
    exit 1
fi
```

### 2. Verify checks are passing (with repair attempts)
```bash
# Get check results - fixed jq query for proper array handling
FAILED_CHECKS=$(gh pr checks --json state | jq '[.[] | select(.state != "SUCCESS")] | length')

if [ "$FAILED_CHECKS" != "0" ]; then
    echo "❌ $FAILED_CHECKS checks are not passing:"
    gh pr checks --json state | jq '.[] | select(.state != "SUCCESS") | {name, state}'
    echo "🔧 Attempting to repair common issues..."

    # Get failed check names for targeted repairs
    FAILED_CHECK_NAMES=$(gh pr checks --json state | jq -r '.[] | select(.state != "SUCCESS") | .name' | tr '\n' ' ')

    # Attempt common repairs
    REPAIR_ATTEMPTED=false

    # Check for style/lint issues
    if echo "$FAILED_CHECK_NAMES" | grep -q -E "Lint|Style|format"; then
        echo "📝 Attempting style fixes..."
        uv run pre-commit run --all-files || true
        REPAIR_ATTEMPTED=true
    fi

    # Check for build issues that might need rebuild
    if echo "$FAILED_CHECK_NAMES" | grep -q -E "Build|Compile"; then
        echo "🔨 Attempting clean rebuild..."
        if [ -f "CMakeLists.txt" ]; then
            rm -rf build && mkdir build && cd build && cmake .. && make -j$(nproc) || true
            cd ..
        fi
        REPAIR_ATTEMPTED=true
    fi

    # Check for Python import issues (might need reinstall)
    if echo "$FAILED_CHECK_NAMES" | grep -q -E "Python|Import"; then
        echo "🐍 Attempting Python package reinstall..."
        uv pip install -e . || true
        REPAIR_ATTEMPTED=true
    fi

    # Check for dependency issues
    if echo "$FAILED_CHECK_NAMES" | grep -q -E "dependency|Dependency"; then
        echo "📦 Attempting dependency refresh..."
        uv sync || true
        REPAIR_ATTEMPTED=true
    fi

    if [ "$REPAIR_ATTEMPTED" = true ]; then
        echo "🔄 Repair attempts made, committing any fixes..."
        git add . || true
        if git diff --staged --quiet; then
            echo "No changes to commit from repairs."
        else
            git commit -m "fix: automated repair of CI issues" || true
            git push || true
        fi

        echo "⏳ Waiting for checks to re-run after repairs..."
        sleep 30

        # Re-check status with fixed query
        echo "🔍 Re-checking CI status after repairs..."
        FAILED_CHECKS_AFTER=$(gh pr checks --json state | jq '[.[] | select(.state != "SUCCESS")] | length')

        if [ "$FAILED_CHECKS_AFTER" != "0" ]; then
            echo "❌ $FAILED_CHECKS_AFTER checks still failing after repair attempts:"
            gh pr checks --json state | jq '.[] | select(.state != "SUCCESS") | {name, state}'
            echo "🚫 Cannot merge until all checks pass. Manual intervention required."
            exit 1
        else
            echo "✅ All checks now passing after repairs!"
        fi
    else
        echo "🚫 No applicable repair strategies found. Manual intervention required."
        exit 1
    fi
else
    echo "✅ All checks are passing"
fi
```

### 3. Perform squash merge with automatic cleanup
```bash
echo "🔄 Performing squash merge with automatic branch cleanup..."

# Get current branch name for confirmation
CURRENT_BRANCH=$(git branch --show-current)

# Squash merge with automatic branch deletion
gh pr merge --squash --delete-branch

if [ $? -eq 0 ]; then
    echo "✅ Squash merge completed successfully"
    echo "🧹 Branch $CURRENT_BRANCH automatically cleaned up"
else
    echo "❌ Merge failed"
    exit 1
fi
```

### 4. Update local repository state
```bash
# Get main branch name
MAIN_BRANCH=$(git remote show origin | awk '/HEAD branch/ {print $NF}')

echo "📥 Updating local repository state..."

# Switch to main branch and pull latest changes
git checkout "$MAIN_BRANCH" || git checkout main || git checkout master
git pull origin "$MAIN_BRANCH"

# Clean up any remaining local branch if it still exists
if git show-ref --verify --quiet "refs/heads/$CURRENT_BRANCH"; then
    git branch -D "$CURRENT_BRANCH"
    echo "🧹 Cleaned up remaining local branch: $CURRENT_BRANCH"
fi

echo "✅ Local repository updated"
```

### 5. Confirm merge
```bash
echo "🎉 Merge completed successfully!"
echo "PR #$PR_NUMBER has been squash merged into $MAIN_BRANCH"
echo "Branch $CURRENT_BRANCH has been cleaned up (both remote and local)"
```

## Usage

```bash
# After PR checks pass, run this workflow
# It will automatically detect the current PR and merge if ready

# Manual execution:
./path/to/this/workflow.sh

# Or run the individual steps:
gh pr view --json number,mergeStateStatus
gh pr checks
gh pr merge --squash --delete-branch  # Handles both merge and remote cleanup
git checkout main && git pull && git branch -D feature-branch  # Local cleanup only
```

## Exit Codes

- **0**: Success - PR merged and branch cleaned up
- **1**: Failure - Conflicts detected, checks failing, or merge failed

## Safety Features

- **Conflict Detection**: Checks for merge conflicts before attempting merge
- **Check Verification**: Ensures all GitHub Actions checks are passing
- **Automated Repair**: Attempts to fix common CI issues before failing
- **Branch Protection**: Respects branch protection rules
- **Clean Exit**: Leaves repository in consistent state if any step fails

## Performance Optimizations

- **Single PR info call**: Combines multiple `gh pr view` calls into one for efficiency
- **Early exits**: Checks conflicts and mergeability first to fail fast
- **Fixed jq queries**: Uses correct array syntax to avoid parsing errors
- **Graceful cleanup**: Handles cases where branches are already cleaned up
- **Targeted repairs**: Only attempts repairs relevant to failed checks

## Notes

- Requires `gh` CLI with appropriate permissions
- Works with any branch name (auto-detects current branch)
- Preserves main branch name (auto-detects from remote)
- Uses squash merge to keep commit history clean
- Uses `gh pr merge --delete-branch` for automatic remote cleanup
- Handles local repository state cleanup separately

## Error Handling

If merge conflicts are detected:
```
❌ Merge conflicts detected in PR #58
Please resolve conflicts manually before merging.
```

If checks are not passing (with repair attempts):
```
❌ Some checks are not passing:
- CI/Lint & Style Checks
🔧 Attempting to repair common issues...
📝 Attempting style fixes...
🔄 Repair attempts made, committing any fixes...
⏳ Waiting for checks to re-run after repairs...
✅ All checks now passing after repairs!
```

If repairs fail:
```
❌ Checks still failing after repair attempts:
- CI/Build & Test C++ (ubuntu-latest)
🚫 Cannot merge until all checks pass. Manual intervention required.
```

If no repair strategies apply:
```
❌ Some checks are not passing:
- Custom Validation Test
🚫 No applicable repair strategies found. Manual intervention required.
```

## Troubleshooting

### Common Issues and Solutions

**jq query errors**:
```bash
# If you get "Cannot index array with string" errors:
# The workflow now uses the correct array syntax: '[.[] | select(.state != "SUCCESS")] | length'

# Alternative manual check:
gh pr checks --json state | jq '.[0].state'  # Check first check structure
```

**gh CLI authentication issues**:
```bash
# Re-authenticate if needed:
gh auth login

# Check current auth status:
gh auth status
```

**Branch already cleaned up**:
```bash
# The workflow handles this gracefully - it will confirm branch is already cleaned up
```

**Merge conflicts**:
```bash
# Manual resolution required:
git checkout main
git pull
git checkout feature-branch
git merge main
# Resolve conflicts, then:
git add . && git commit && git push
```
