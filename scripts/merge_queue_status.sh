#!/bin/bash
# scripts/merge_queue_status.sh
# Read-only view of the merge queue on main: one fixed GraphQL query, no
# mutation, so it can be allow-listed on its own (a general `gh api` rule
# would also admit writes).
#
#   scripts/merge_queue_status.sh        # every entry: position, PR, state
#   scripts/merge_queue_status.sh 978    # one PR: queued (position, state) or not
#
# Entry states: QUEUED, AWAITING_CHECKS, MERGEABLE, UNMERGEABLE, LOCKED.

set -euo pipefail

OWNER="thirox-lab"
NAME="combaero"

if [[ $# -gt 1 || ( $# -eq 1 && ! "$1" =~ ^[0-9]+$ ) ]]; then
    echo "Usage: $0 [PR]"; exit 1
fi

entries=$(gh api graphql -f owner="$OWNER" -f name="$NAME" -f query='
  query($owner: String!, $name: String!) {
    repository(owner: $owner, name: $name) {
      mergeQueue(branch: "main") {
        entries(first: 50) {
          nodes { position state enqueuedAt pullRequest { number title } }
        }
      }
    }
  }' --jq '.data.repository.mergeQueue.entries.nodes[]
           | "\(.position)\t#\(.pullRequest.number)\t\(.state)\t\(.enqueuedAt)\t\(.pullRequest.title)"')

if [[ $# -eq 1 ]]; then
    line=$(printf '%s\n' "$entries" | awk -F'\t' -v pr="#$1" '$2 == pr')
    if [[ -n "$line" ]]; then
        printf '%s\n' "$line" | awk -F'\t' '{print $2 " queued: position " $1 ", " $3}'
    else
        echo "#$1 not in the queue (merged, dequeued or never queued)"
    fi
else
    [[ -z "$entries" ]] && { echo "merge queue on main: empty"; exit 0; }
    printf 'pos\tPR\tstate\tenqueued\ttitle\n%s\n' "$entries"
fi
