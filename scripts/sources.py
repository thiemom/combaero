"""The source register: every local reference document, what it holds and
what has been done with it.

The documents themselves live under docs/ and are gitignored (copyright);
validation/sources.yaml records, per file, its citation, what it contains,
which questions it has been read for and with what verdict, and where the
repo uses it. One place to answer "do we already have this?", "what did we
conclude?" and "what else covers this?".

    uv run scripts/sources.py find impingement effusion   # search (all terms)
    uv run scripts/sources.py list --status unscreened    # filter
    uv run scripts/sources.py list --topic film --status implemented
    uv run scripts/sources.py check                       # schema + local docs
    uv run scripts/sources.py new docs/x/new_paper.pdf    # print a stub entry

`check` validates the register (schema, statuses, that every record and
code path it names exists) and, where docs/ is present, lists local
documents with no entry and entries whose file is missing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
REGISTER = ROOT / "validation" / "sources.yaml"
DOCS = ROOT / "docs"
DOC_SUFFIXES = {".pdf", ".djvu", ".epub"}

STATUSES = {
    "implemented": "a shipped closure or kernel is built from it",
    "fidelity": "scored as fidelity on its own data",
    "accuracy": "used or scored as a cross-source accuracy set",
    "extracted": "data transcribed or digitised into the repo, not yet scored",
    "screened": "read, verdict recorded, not used yet",
    "unscreened": "nobody has read it for this project yet",
}
TOPICS = {
    "film",
    "effusion",
    "impingement",
    "pin_fin",
    "rib",
    "dimple",
    "channel",
    "junction",
    "tee",
    "ejector",
    "orifice",
    "combustion",
    "area_change",
    "vortex",
    "acoustics",
    "transport",
    "review",
    "other",
}
REQUIRED = ("file", "citation", "topic", "content", "status", "verdict", "read_for", "issues", "records")
OPTIONAL = ("id", "code", "notes")


def load(path: Path = REGISTER) -> list[dict[str, Any]]:
    entries = yaml.safe_load(path.read_text())
    if not isinstance(entries, list):
        raise SystemExit(f"{path}: expected a YAML list of entries")
    return entries


def problems(entries: list[dict[str, Any]]) -> list[str]:
    """Schema errors in the register; empty when it is valid."""
    out: list[str] = []
    seen: set[str] = set()
    for i, e in enumerate(entries):
        where = f"entry {i} ({e.get('file', '?')})"
        if not isinstance(e, dict):
            out.append(f"{where}: not a mapping")
            continue
        for k in REQUIRED:
            if k not in e:
                out.append(f"{where}: missing '{k}'")
        for k in e:
            if k not in REQUIRED and k not in OPTIONAL:
                out.append(f"{where}: unknown key '{k}'")
        f = str(e.get("file", ""))
        if not f.startswith("docs/") or Path(f).suffix.lower() not in DOC_SUFFIXES:
            out.append(f"{where}: file must be a document under docs/")
        if f in seen:
            out.append(f"{where}: duplicate file")
        seen.add(f)
        if e.get("status") not in STATUSES:
            out.append(f"{where}: status '{e.get('status')}' not one of {sorted(STATUSES)}")
        topics = e.get("topic") or []
        if not isinstance(topics, list) or not topics or any(t not in TOPICS for t in topics):
            out.append(f"{where}: topic must be a non-empty list from {sorted(TOPICS)}")
        for k in ("read_for", "issues", "records"):
            if k in e and not isinstance(e[k], list):
                out.append(f"{where}: '{k}' must be a list")
        if any(not isinstance(n, int) for n in e.get("issues") or []):
            out.append(f"{where}: issues are integers")
        for r in e.get("records") or []:
            if not (ROOT / str(r)).exists():
                out.append(f"{where}: record '{r}' does not exist")
        for c in e.get("code") or []:
            path = str(c).split(":", 1)[0].strip()
            if not (ROOT / path).exists():
                out.append(f"{where}: code path '{path}' does not exist")
        if e.get("status") in ("implemented", "fidelity", "extracted") and not (
            e.get("records") or e.get("code")
        ):
            out.append(f"{where}: status '{e.get('status')}' needs a record or code reference")
    return out


def local_documents(docs: Path = DOCS) -> list[str]:
    """Documents under `docs` (the repo's docs/ by default; a worktree has
    none, so point it at the main checkout's), as paths starting docs/."""
    if not docs.is_dir():
        return []
    return sorted(
        "docs/" + str(p.relative_to(docs)) for p in docs.rglob("*") if p.suffix.lower() in DOC_SUFFIXES and p.is_file()
    )


def show(e: dict[str, Any], long: bool = False) -> str:
    head = f"[{e['status']}] {e['file']}\n    {e['citation']}"
    if not long:
        return head + f"\n    -> {e['verdict']}"
    lines = [head, f"    topic: {', '.join(e['topic'])}", f"    content: {e['content']}", f"    verdict: {e['verdict']}"]
    for k in ("read_for", "issues", "records", "code"):
        if e.get(k):
            lines.append(f"    {k}: {', '.join(str(v) for v in e[k])}")
    if e.get("notes"):
        lines.append(f"    notes: {e['notes']}")
    return "\n".join(lines)


def cmd_find(args: argparse.Namespace) -> int:
    terms = [t.lower() for t in args.terms]
    hits = [
        e
        for e in load()
        if all(t in yaml.safe_dump(e, allow_unicode=True).lower() for t in terms)
    ]
    for e in hits:
        print(show(e, args.long))
    print(f"{len(hits)} match(es)")
    return 0


def cmd_list(args: argparse.Namespace) -> int:
    rows = [
        e
        for e in load()
        if (not args.status or e["status"] == args.status) and (not args.topic or args.topic in e["topic"])
    ]
    for e in rows:
        print(show(e, args.long))
    print(f"{len(rows)} entr{'y' if len(rows) == 1 else 'ies'}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    entries = load()
    errs = problems(entries)
    for p in errs:
        print(f"error: {p}")
    local = local_documents(Path(args.docs) if args.docs else DOCS)
    if local:
        listed = {e["file"] for e in entries}
        for f in local:
            if f not in listed:
                print(f"unregistered: {f}")
        for f in sorted(listed - set(local)):
            print(f"not found locally: {f}")
    else:
        print("docs/ not present: local coverage not checked")
    counts: dict[str, int] = {}
    for e in entries:
        counts[e["status"]] = counts.get(e["status"], 0) + 1
    print(f"{len(entries)} entries: " + ", ".join(f"{k} {counts.get(k, 0)}" for k in STATUSES))
    return 1 if errs else 0


def cmd_new(args: argparse.Namespace) -> int:
    stub = {
        "file": args.file,
        "citation": "?",
        "topic": ["other"],
        "content": "?",
        "status": "unscreened",
        "verdict": "unscreened",
        "read_for": [],
        "issues": [],
        "records": [],
    }
    print(yaml.safe_dump([stub], sort_keys=False, allow_unicode=True), end="")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("find", help="entries matching ALL terms (any field)")
    f.add_argument("terms", nargs="+")
    f.add_argument("--long", action="store_true")
    f.set_defaults(fn=cmd_find)
    ls_ = sub.add_parser("list", help="entries by status and/or topic")
    ls_.add_argument("--status", choices=sorted(STATUSES))
    ls_.add_argument("--topic", choices=sorted(TOPICS))
    ls_.add_argument("--long", action="store_true")
    ls_.set_defaults(fn=cmd_list)
    c = sub.add_parser("check", help="validate the register; list unregistered local documents")
    c.add_argument("--docs", help="the docs/ directory to audit (default: this checkout's)")
    c.set_defaults(fn=cmd_check)
    n = sub.add_parser("new", help="print a stub entry for a new document")
    n.add_argument("file")
    n.set_defaults(fn=cmd_new)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
