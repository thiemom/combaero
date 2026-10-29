"""One-shot migration (#846): move untracked files under docs/ into sources/.

Copyrighted third-party material (papers, handbooks, scans, figure crops,
digitising folders, supplementary archives) used to sit in gitignored
folders under docs/. It now lives in sources/, which git ignores whatever
the file type, while docs/ became an allowlist for our own .md files.

Every file under docs/ that git does not track moves to the same subpath
under sources/ (docs/materials/x.pdf -> sources/materials/x.pdf), with one
rename: docs/sources/ (page scans) becomes sources/scans/. Tracked files
are never touched, an existing destination is never overwritten, and
nothing moves without --apply.

Untracked .md files move too: under the new allowlist an untracked .md in
docs/ is no longer ignored, and the ones there today are transcriptions and
working notes about the sources. Run this right after pulling the change.

    uv run --no-project scripts/migrate_docs_to_sources.py               # dry run, this checkout
    uv run --no-project scripts/migrate_docs_to_sources.py --root <main checkout>
    uv run --no-project scripts/migrate_docs_to_sources.py --apply
    uv run --no-project scripts/migrate_docs_to_sources.py --list        # every move

Afterwards: `uv run scripts/sources.py check` should be clean.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections import Counter
from pathlib import Path

# Top-level docs/ folders that change name on the way (sources/sources/
# would read badly).
RENAME = {"sources": "scans"}
# Finder litter: deleted, not moved.
JUNK = {".DS_Store"}


def tracked(root: Path) -> set[str]:
    r = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z", "--", "docs"],
        capture_output=True,
        text=True,
        check=True,
    )
    return {f for f in r.stdout.split("\0") if f}


def destination(rel: Path) -> Path:
    """sources/... for a path relative to docs/."""
    parts = list(rel.parts)
    if len(parts) > 1:
        parts[0] = RENAME.get(parts[0], parts[0])
    return Path("sources", *parts)


def plan(root: Path) -> tuple[list[tuple[str, str]], list[str], list[str]]:
    """(moves, junk, conflicts): moves as (src, dst) repo-relative POSIX paths."""
    docs = root / "docs"
    keep = tracked(root)
    moves: list[tuple[str, str]] = []
    junk: list[str] = []
    conflicts: list[str] = []
    for p in sorted(docs.rglob("*")):
        if not p.is_file() and not p.is_symlink():
            continue
        src = p.relative_to(root).as_posix()
        if src in keep:
            continue
        if p.name in JUNK:
            junk.append(src)
            continue
        dst = destination(p.relative_to(docs))
        if (root / dst).exists() or (root / dst).is_symlink():
            conflicts.append(f"{src} -> {dst.as_posix()} (destination exists)")
            continue
        moves.append((src, dst.as_posix()))
    return moves, junk, conflicts


def group(path: str) -> str:
    """docs/<topic> for a file in a folder, else docs/ itself."""
    parts = path.split("/")
    return "/".join(parts[:2]) + "/" if len(parts) > 2 else "docs/"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--root",
        default=str(Path(__file__).resolve().parent.parent),
        help="repository checkout to migrate (default: this one)",
    )
    ap.add_argument("--apply", action="store_true", help="move the files (default: dry run)")
    ap.add_argument("--list", action="store_true", help="print every move, not only counts")
    args = ap.parse_args(argv)
    root = Path(args.root).resolve()
    if not (root / "docs").is_dir():
        print(f"{root}/docs: not found, nothing to do")
        return 0
    moves, junk, conflicts = plan(root)

    by_group = Counter(group(src) for src, _ in moves)
    by_suffix = Counter(Path(src).suffix.lower() or "(none)" for src, _ in moves)
    for g in sorted(by_group):
        dst = destination(Path(g).relative_to("docs") / "x").parent.as_posix() + "/"
        print(f"{by_group[g]:5d}  {g:28s} -> {dst if g != 'docs/' else 'sources/'}")
    print(
        f"{len(moves):5d}  files to move; by type: "
        + ", ".join(f"{k} {v}" for k, v in by_suffix.most_common())
    )
    print(f"{len(junk):5d}  {', '.join(sorted(JUNK))} to delete")
    if args.list:
        for src, dst in moves:
            print(f"  {src} -> {dst}")
    for c in conflicts:
        print(f"conflict: {c}")
    if conflicts:
        print(f"{len(conflicts)} conflict(s): resolve them by hand; nothing was moved")
        return 1
    if not args.apply:
        print("dry run: nothing moved (pass --apply)")
        return 0

    for src, dst in moves:
        target = root / dst
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():  # appeared since planning
            print(f"refusing to overwrite {dst}; stopped")
            return 1
        (root / src).rename(target)
    for j in junk:
        (root / j).unlink()
    # Prune the now-empty folders (tracked ones such as docs/archive/ keep
    # their files and stay).
    for d in sorted((root / "docs").rglob("*"), key=lambda p: len(p.parts), reverse=True):
        if d.is_dir() and not any(d.iterdir()):
            d.rmdir()
    print(f"moved {len(moves)} file(s), deleted {len(junk)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
