#!/bin/bash
# scripts/bulk_merge_safe.sh
# Automates the discovery and merging of "safe" Dependabot updates (minor/patch/grouped).

GH_BIN="/opt/homebrew/bin/gh"
DRY_RUN=false
VERBOSE=false

usage() {
    echo "Usage: $0 [--dry-run] [--verbose]"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN=true ;;
        --verbose) VERBOSE=true ;;
        *) usage ;;
    esac
    shift
done

echo "🔍 Fetching open Dependabot PRs..."
PRS=$($GH_BIN pr list --author "app/dependabot" --json number,title,statusCheckRollup,mergeable --limit 50)

# Function to extract SemVer change type from title
# Example: "bump foo from 1.0.0 to 1.1.0"
get_bump_type() {
    local title="$1"
    # match patterns:
    # "from 1.0.0 to 1.1.0"
    # "requirement from >=1.0 to >=1.1"
    if [[ "$title" =~ from[[:space:]]([>=]*[0-9.]+)[[:space:]]to[[:space:]]([>=]*[0-9.]+) ]]; then
        local old_ver=$(echo "${BASH_REMATCH[1]}" | tr -d '>=')
        local new_ver=$(echo "${BASH_REMATCH[2]}" | tr -d '>=')

        # Simple extraction of first two parts for comparison
        # (This avoids complex IFS logic failures on short versions like '1.0')
        local o1=$(echo "$old_ver" | cut -d. -f1)
        local o2=$(echo "$old_ver" | cut -d. -f2 || echo "0")
        local o3=$(echo "$old_ver" | cut -d. -f3 || echo "0")

        local n1=$(echo "$new_ver" | cut -d. -f1)
        local n2=$(echo "$new_ver" | cut -d. -f2 || echo "0")
        local n3=$(echo "$new_ver" | cut -d. -f3 || echo "0")

        if [[ "$n1" -gt "$o1" ]]; then echo "major";
        elif [[ "$n2" -gt "$o2" ]]; then echo "minor";
        elif [[ "$n3" -gt "$o3" ]]; then echo "patch";
        else echo "patch"; fi # Assume patch if syntax matches but versions same or just ambiguous
    elif [[ "$title" == *"safe-updates"* || "$title" == *"actions group"* ]]; then
        echo "grouped-safe"
    else
        echo "unknown"
    fi
}

echo "--------------------------------------------------"
printf "%-7s | %-10s | %-10s | %s\n" "PR #" "TYPE" "STATUS" "TITLE"
echo "--------------------------------------------------"

candidates=()

while read -r pr; do
    num=$(echo "$pr" | jq -r '.number')
    title=$(echo "$pr" | jq -r '.title')
    status=$(echo "$pr" | jq -r '.statusCheckRollup // [] | if length == 0 then "PENDING" else (map(.conclusion) | if all(. == "SUCCESS") then "SUCCESS" else "FAIL" end) end')
    mergeable=$(echo "$pr" | jq -r '.mergeable')

    bump_type=$(get_bump_type "$title")

    # Logic for "Safe":
    # 1. Grouped "safe-updates"
    # 2. Individual "minor" or "patch"
    is_safe=false
    if [[ "$title" == *"safe-updates"* ]]; then is_safe=true; fi
    if [[ "$bump_type" == "minor" || "$bump_type" == "patch" ]]; then is_safe=true; fi

    printf "%-7s | %-10s | %-10s | %s\n" "#$num" "$bump_type" "$status" "$title"

    if [[ "$is_safe" == "true" && "$status" == "SUCCESS" && "$mergeable" == "MERGEABLE" ]]; then
        candidates+=("$num")
    fi
done < <(echo "$PRS" | jq -c '.[]')

echo "--------------------------------------------------"

if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "✅ No safe, verified candidates found for merging."
    exit 0
fi

if [[ "$DRY_RUN" == "true" ]]; then
    echo "📝 Dry-run complete. Found ${#candidates[@]} candidates for bulk merge: ${candidates[*]}"
    exit 0
fi

echo "🚀 Starting bulk squash merge of ${#candidates[@]} PRs..."
for pr_num in "${candidates[@]}"; do
    echo "➡️ Merging #$pr_num..."
    if $GH_BIN pr merge "$pr_num" --squash --delete-branch; then
        echo "✅ Merged #$pr_num"
    else
        echo "⚠️ Failed to merge #$pr_num (likely conflict). Skipping..."
    fi
done

echo "--------------------------------------------------"
echo "🔄 Synchronizing local main branch..."
git checkout main && git pull
echo "✨ Bulk merge complete."
