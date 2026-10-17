"""Summarize a stress run (scripts/stress/run.py), optionally against another.

    uv run python scripts/stress/summarize.py tmp/stress/after.jsonl
    uv run python scripts/stress/summarize.py tmp/stress/after.jsonl --compare tmp/stress/before.jsonl

Prints per family: solves, converged, failed, raised out of solve(), caught
exceptions, slowest solve, outcomes and row-kind verdicts; then every
DISAGREE row (a misfiled residual row, RESIDUAL-ROW-SCALING), every
exception out of solve() and the caught exceptions grouped by site. With
``--compare``: per family the converged count before -> after and the seeds
gained and lost. A "converged" result is not yet a physical one: look at
implausible internal states too (#619: zero flow at P_py ~ 1e13 Pa).
"""

from __future__ import annotations

import argparse
import collections
import json
import re
from pathlib import Path
from typing import Any


def load(path: Path) -> list[dict[str, Any]]:
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def _norm(msg: str, n: int) -> str:
    return re.sub(r"[-+]?\d[\d.e+-]*", "#", msg)[:n]


def summarize(recs: list[dict[str, Any]]) -> None:
    by: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for r in recs:
        by[r["family"]].append(r)
    print(
        f"{'family':14s} {'n':>4s} {'conv':>5s} {'fail':>5s} {'raise':>5s} {'caught':>6s} {'t_max':>6s}"
        "  outcomes / row kinds"
    )
    for fam, rs in sorted(by.items()):
        conv = sum(bool(r.get("converged")) for r in rs)
        raised = sum("exc" in r for r in rs)
        caught = sum(r.get("caught_count", 0) for r in rs)
        oc = collections.Counter(
            r.get("outcome") for r in rs if not r.get("converged") and "exc" not in r
        )
        rk: collections.Counter = collections.Counter()
        for r in rs:
            rk.update(r.get("rowkind") or {})
        print(
            f"{fam:14s} {len(rs):4d} {conv:5d} {len(rs) - conv - raised:5d} {raised:5d} {caught:6d} "
            f"{max(r['dt'] for r in rs):6.1f}  {dict(oc)} {dict(rk)}"
        )
    print("\nDISAGREE rows (misfiled residual units):")
    for r in recs:
        for b in r.get("rowkind_bad") or []:
            print(f"   {r['family']} seed {r['seed']}: {b}")
        if r.get("rowkind_exc"):
            print(f"   check raised: {r['family']} seed {r['seed']}: {r['rowkind_exc']}")
    print("\nExceptions out of solve():")
    ex: dict[str, list] = collections.defaultdict(list)
    for r in recs:
        if "exc" in r:
            ex[_norm(r["exc"], 140)].append((r["family"], r["seed"]))
    for k, v in sorted(ex.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(v):4d} {k}\n       e.g. {v[:4]}")
    print("\nCaught exceptions, grouped (phase, site, type, message, innermost frame):")
    g: dict[tuple, list] = collections.defaultdict(list)
    for r in recs:
        for c in r.get("caught", []):
            key = (
                c["phase"],
                (c["where"] or "-").split(" ")[0],
                c["type"],
                _norm(c["message"], 90),
                (c["frames"] or [""])[-1],
            )
            g[key].append((r["family"], r["seed"]))
    for k, v in sorted(g.items(), key=lambda kv: -len(kv[1])):
        print(
            f"  {len(v):4d} {k[0]:14s} {k[1]:8s} {k[2]:12s} {k[3]}\n       at {k[4]}  e.g. {sorted(set(v))[:4]}"
        )


def compare(after: list[dict[str, Any]], before: list[dict[str, Any]]) -> None:
    def conv(recs: list[dict[str, Any]]) -> dict[tuple[str, int], bool]:
        return {(r["family"], r["seed"]): bool(r.get("converged")) for r in recs}

    a, b = conv(after), conv(before)
    print("\nAgainst the baseline (same family and seed in both):")
    for fam in sorted({k[0] for k in a} & {k[0] for k in b}):
        keys = sorted(k for k in a if k[0] == fam and k in b)
        gained = [k[1] for k in keys if a[k] and not b[k]]
        lost = [k[1] for k in keys if b[k] and not a[k]]
        print(
            f"  {fam:14s} converged {sum(b[k] for k in keys):3d} -> {sum(a[k] for k in keys):3d}"
            f"   +{len(gained)} {gained[:10]}   -{len(lost)} {lost[:10]}"
        )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("run", type=Path)
    p.add_argument("--compare", type=Path, help="a baseline run to compare against")
    a = p.parse_args(argv)
    after = load(a.run)
    summarize(after)
    if a.compare:
        compare(after, load(a.compare))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
