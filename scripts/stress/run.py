"""Solve stress networks cold and record what happened, one JSON line each.

    uv run python scripts/stress/run.py --out tmp/stress/after.jsonl
    uv run python scripts/stress/run.py --families ejector,tee3 --seeds 0:100 --jobs 8 --out ...
    uv run python scripts/stress/summarize.py tmp/stress/after.jsonl --compare tmp/stress/before.jsonl

Per solve: converged/success/outcome, the caught-exception log
(``__caught_exceptions__``, #615), an exception that escaped ``solve()``,
inconsistent elements, the second-law verdicts and ledger, and the
RESIDUAL-ROW-SCALING check (``check_row_kinds``) at the solution. Lenient by
default, so a run collects every finding; ``--strict`` raises at the first
(the suite's mode, COMBAERO_SOLVER_STRICT).

The solver's stall detector and timeouts are wall-clock: parallel jobs
share the CPU and can flip a slow case. For an A/B comparison near the
timeout run both sides with the same ``--jobs`` (``--jobs 1`` to be sure).
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import multiprocessing as mp
import sys
import time
import traceback
import warnings
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from networks import FAMILIES  # noqa: E402


def _clean(v: Any) -> Any:
    return repr(v) if isinstance(v, float) and not math.isfinite(v) else v


def run_one(task: tuple[str, int, float, bool]) -> dict[str, Any]:
    family, seed, timeout, strict = task
    from combaero.network import NetworkSolver

    rec: dict[str, Any] = {"family": family, "seed": seed}
    t0 = time.time()
    try:
        net, tags = FAMILIES[family](seed)
        rec["tags"] = sorted(set(tags))
        rec["n_el"] = len(net.elements)
        rec["classes"] = sorted({type(e).__name__ for e in net.elements.values()})
        solver = NetworkSolver(net)
        solver.raise_internal_errors = strict
        # As the GUI solves: with the Newton-homotopy and boundary-condition
        # ramp retries (#642).
        solver.homotopy_retry = True
        solver.bc_ramp_retry = True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = solver.solve(timeout=timeout)
        rec["converged"] = bool(r.get("__converged__"))
        rec["success"] = bool(r.get("__success__"))
        rec["outcome"] = str(r.get("__outcome__"))
        rec["final_norm"] = _clean(r.get("__final_norm__"))
        rec["caught"] = r.get("__caught_exceptions__", [])
        rec["caught_count"] = r.get("__caught_exception_count__", 0)
        rec["inconsistent"] = list(r.get("__inconsistent_elements__") or [])
        rec["second_law"] = {
            eid: {k: _clean(x) for k, x in v.items() if k != "sigma"}
            for eid, v in (r.get("__second_law__") or {}).items()
        }
        rec["ledger"] = {
            k: _clean(v) for k, v in (r.get("__ledger__") or {}).items() if k != "streams"
        }
        if rec["converged"] and "__x_solution__" in r:
            try:
                rep = solver.check_row_kinds(r["__x_solution__"])
                rec["rowkind"] = dict(collections.Counter(x["verdict"] for x in rep))
                rec["rowkind_bad"] = [x for x in rep if x["verdict"] == "DISAGREE"]
            except Exception as exc:  # noqa: BLE001 -- recorded, never fatal
                rec["rowkind_exc"] = f"{type(exc).__name__}: {exc}"[:300]
    except Exception as exc:  # noqa: BLE001 -- an exception out of solve() is a finding
        rec["exc"] = f"{type(exc).__name__}: {exc}"[:300]
        rec["tb_tail"] = traceback.format_exc()[-1500:]
    rec["dt"] = round(time.time() - t0, 2)
    return rec


def _seeds(spec: str) -> range:
    lo, _, hi = spec.partition(":")
    return range(int(lo), int(hi)) if hi else range(int(lo))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--families", default="all", help=f"comma list or 'all' ({', '.join(FAMILIES)})")
    p.add_argument("--seeds", default="0:40", help="'lo:hi' or a count (default 0:40)")
    p.add_argument("--jobs", type=int, default=4, help="worker processes (default 4)")
    p.add_argument("--timeout", type=float, default=20.0, help="per-solve timeout [s] (default 20)")
    p.add_argument("--strict", action="store_true", help="raise at the first finding (suite mode)")
    p.add_argument("--out", type=Path, required=True, help="JSONL output (appended)")
    a = p.parse_args(argv)
    fams = list(FAMILIES) if a.families == "all" else a.families.split(",")
    unknown = [f for f in fams if f not in FAMILIES]
    if unknown:
        p.error(f"unknown families {unknown}; known: {list(FAMILIES)}")
    tasks = [(f, s, a.timeout, a.strict) for s in _seeds(a.seeds) for f in fams]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(a.out, "a") as fh, mp.get_context("spawn").Pool(a.jobs) as pool:
        for i, rec in enumerate(pool.imap_unordered(run_one, tasks), 1):
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            if i % 50 == 0 or i == len(tasks):
                print(f"{i}/{len(tasks)} solves, {time.time() - t0:.0f} s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
