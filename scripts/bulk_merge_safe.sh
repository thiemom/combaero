#!/bin/bash
# scripts/bulk_merge_safe.sh
# Merge "safe" Dependabot PRs (grouped safe-updates, or minor/patch bumps),
# ONE AT A TIME: main is branch-protected and requires up-to-date branches,
# so every merge leaves the next PR behind. For each PR, in order:
#   1. author must be app/dependabot and the bump safe;
#   2. if behind main: ask Dependabot to rebase ("@dependabot rebase"), wait
#      for the new head, then for CI. NOT gh pr update-branch: the Iron
#      ruleset requires signed commits and a GitHub-side rebase is unsigned;
#      Dependabot re-signs its own;
#   3. all checks finished, none failed (SKIPPED/NEUTRAL are fine);
#   4. mergeStateStatus CLEAN;
#   5. pinned tool versions agree at the PR head (scripts/check_version_sync.py);
#   6. squash-merge exactly the verified head commit (--match-head-commit).
# Anything else is reported and skipped; nothing is forced.
#
# Usage: scripts/bulk_merge_safe.sh [--dry-run] [PR ...]

set -euo pipefail

GH_BIN="$(command -v gh)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DRY_RUN=false
ONLY=()

for arg in "$@"; do
    case $arg in
        --dry-run) DRY_RUN=true ;;
        [0-9]*) ONLY+=("$arg") ;;
        *) echo "Usage: $0 [--dry-run] [PR ...]"; exit 1 ;;
    esac
done

is_safe_title() {
    local title="$1"
    [[ "$title" == *"safe-updates"* || "$title" == *"actions group"* ]] && return 0
    if [[ "$title" =~ from[[:space:]]v?([0-9]+)\.([0-9]+)[^[:space:]]*[[:space:]]to[[:space:]]v?([0-9]+)\. ]]; then
        [[ "${BASH_REMATCH[1]}" == "${BASH_REMATCH[3]}" ]] && return 0  # same major
    fi
    return 1
}

pr_field() { $GH_BIN pr view "$1" --json "$2" --jq ".$2"; }

checks_state() {
    # PASS when every check finished and none failed; PENDING; or FAIL.
    $GH_BIN pr view "$1" --json statusCheckRollup --jq '
      [.statusCheckRollup[] | (.conclusion // .state // "")] as $c
      | if ($c | length) == 0 then "PENDING"
        elif any($c[]; . == "" or . == "PENDING" or . == "QUEUED" or . == "IN_PROGRESS" or . == "EXPECTED") then "PENDING"
        elif all($c[]; . == "SUCCESS" or . == "SKIPPED" or . == "NEUTRAL") then "PASS"
        else "FAIL" end'
}

wait_for_new_head() {
    # Echo the new head oid once it differs from $2, or TIMEOUT (~20 min).
    local n="$1" old="$2" i=0 h
    while :; do
        sleep 30
        h=$(pr_field "$n" headRefOid)
        [[ "$h" != "$old" ]] && { echo "$h"; return; }
        (( i++ > 40 )) && { echo "TIMEOUT"; return; }
    done
}

wait_for_checks() {
    local n="$1" i=0 s
    sleep 20  # let GitHub register the checks of a new head
    while :; do
        s=$(checks_state "$n")
        [[ "$s" != "PENDING" ]] && { echo "$s"; return; }
        (( i++ > 120 )) && { echo "TIMEOUT"; return; }  # ~60 min
        sleep 30
    done
}

merge_one() {
    local n="$1" title state head
    title=$(pr_field "$n" title)
    echo "--- #$n  $title"
    if [[ "$($GH_BIN pr view "$n" --json author --jq .author.login)" != "app/dependabot" ]]; then
        echo "    skip: not a Dependabot PR"; return 1
    fi
    is_safe_title "$title" || { echo "    skip: not a safe (grouped/minor/patch) bump"; return 1; }

    state=$(pr_field "$n" mergeStateStatus)
    if [[ "$state" == "BEHIND" ]]; then
        if $DRY_RUN; then echo "    would ask Dependabot to rebase, wait for CI, then re-check"; return 0; fi
        echo "    behind main: asking Dependabot to rebase"
        local old; old=$(pr_field "$n" headRefOid)
        $GH_BIN pr comment "$n" --body "@dependabot rebase" >/dev/null
        if [[ "$(wait_for_new_head "$n" "$old")" == "TIMEOUT" ]]; then
            echo "    skip: Dependabot did not rebase (if it says the branch was edited, comment '@dependabot recreate')"
            return 1
        fi
        echo "    rebased; waiting for CI..."
        [[ "$(wait_for_checks "$n")" == "PASS" ]] || { echo "    skip: CI did not pass after rebase"; return 1; }
        state=$(pr_field "$n" mergeStateStatus)
    fi
    if [[ "$(checks_state "$n")" == "PENDING" ]] && ! $DRY_RUN; then
        echo "    CI still running; waiting..."
        wait_for_checks "$n" >/dev/null
        state=$(pr_field "$n" mergeStateStatus)
    fi
    [[ "$(checks_state "$n")" == "PASS" ]] || { echo "    skip: checks not all passed"; return 1; }
    local tries=0
    while [[ "$state" == "UNKNOWN" ]] && (( tries++ < 8 )); do
        sleep 15  # GitHub is still computing mergeability for the new head
        state=$(pr_field "$n" mergeStateStatus)
    done
    [[ "$state" == "CLEAN" ]] || { echo "    skip: mergeStateStatus $state"; return 1; }

    head=$(pr_field "$n" headRefOid)
    git -C "$ROOT" fetch -q origin "$head"
    if ! uv run --project "$ROOT" "$ROOT/scripts/check_version_sync.py" --ref "$head"; then
        echo "    skip: pinned versions disagree at the PR head"; return 1
    fi
    if $DRY_RUN; then echo "    would merge $head"; return 0; fi
    $GH_BIN pr merge "$n" --squash --delete-branch --match-head-commit "$head"
    echo "    merged"
}

if [[ ${#ONLY[@]} -eq 0 ]]; then
    # bash 3.2 (macOS): no mapfile
    while read -r n; do ONLY+=("$n"); done < <($GH_BIN pr list --author "app/dependabot" --json number --jq '.[].number' | sort -n)
fi
[[ ${#ONLY[@]} -eq 0 ]] && { echo "No open Dependabot PRs."; exit 0; }

merged=0; skipped=0
for n in "${ONLY[@]}"; do
    if merge_one "$n"; then merged=$((merged + 1)); else skipped=$((skipped + 1)); fi
done
if $DRY_RUN; then echo "=== dry run: $merged eligible, $skipped skipped"; else echo "=== done: $merged merged, $skipped skipped"; fi
