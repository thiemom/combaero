"""Reject copyrighted third-party material from commits (#846).

Papers, handbooks, scans and figure crops live in sources/, which .gitignore
ignores whatever the file type. A `git add -f`, a `--no-verify` commit or a
web edit gets past .gitignore, so this check runs twice: as a pre-commit
hook on the staged paths and in CI on every tracked file.

A path fails when it is
  * under sources/ (except sources/.gitignore, the second lock), or
  * a document, archive, office or image file outside ALLOWLIST, the places
    where we commit images we made ourselves.

    uv run --no-project scripts/check_no_sources.py PATH...   # pre-commit
    uv run --no-project scripts/check_no_sources.py --all     # every tracked file

Standard library only, so CI runs it without building the package.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
SOURCES_DIR = "sources/"
SOURCES_LOCK = "sources/.gitignore"

# Document, archive, office and raster image types: what a copyrighted
# source arrives as. Lower-case suffixes without the dot; .tar.gz is caught
# by "gz".
BLOCKED_SUFFIXES: frozenset[str] = frozenset(
    {
        "pdf",
        "djvu",
        "epub",
        "zip",
        "tar",
        "tgz",
        "gz",
        "7z",
        "rar",
        "doc",
        "docx",
        "xls",
        "xlsx",
        "ppt",
        "pptx",
        "png",
        "jpg",
        "jpeg",
        "tif",
        "tiff",
        "gif",
        "bmp",
        "webp",
    }
)

IMAGE_SUFFIXES: frozenset[str] = frozenset(
    {"png", "jpg", "jpeg", "gif", "webp", "tif", "tiff", "bmp"}
)

# The only places a blocked type may be committed: images we made. Directory
# prefixes (POSIX, relative to the repository root) mapped to the suffixes
# allowed there. Add a directory here, with a reason, before committing an
# image of our own anywhere else.
ALLOWLIST: dict[str, frozenset[str]] = {
    # plots the validation harness renders from digitised coordinates
    "validation/cooling/plots/": IMAGE_SUFFIXES,
    # README banner
    ".github/assets/": IMAGE_SUFFIXES,
    # GUI artwork
    "gui/frontend/public/": IMAGE_SUFFIXES,
    "gui/frontend/src/assets/": IMAGE_SUFFIXES,
}


def problem(path: str) -> str | None:
    """Why `path` must not be committed, or None when it may be."""
    p = path.replace("\\", "/")
    if p.startswith("./"):
        p = p[2:]
    if p == SOURCES_LOCK:
        return None
    if p.startswith(SOURCES_DIR):
        return "under sources/ (copyrighted third-party material, never committed)"
    suffix = PurePosixPath(p).suffix.lower().lstrip(".")
    if suffix not in BLOCKED_SUFFIXES:
        return None
    for prefix, allowed in ALLOWLIST.items():
        if p.startswith(prefix) and suffix in allowed:
            return None
    return (
        f".{suffix} outside the allowlist: a copyrighted source belongs in sources/; "
        "an image we made needs its directory added to ALLOWLIST in scripts/check_no_sources.py"
    )


def problems(paths: list[str]) -> list[str]:
    out = []
    for path in paths:
        why = problem(path)
        if why:
            out.append(f"{path}: {why}")
    return out


def tracked_files() -> list[str]:
    """Every tracked path, relative to the repository root whatever the cwd
    (CI runs pytest from python/tests)."""
    r = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True
    )
    return [f for f in r.stdout.split("\0") if f]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("paths", nargs="*", help="paths to check (pre-commit passes the staged ones)")
    ap.add_argument("--all", action="store_true", help="check every tracked file")
    args = ap.parse_args(argv)
    paths = list(args.paths) + (tracked_files() if args.all else [])
    errs = problems(paths)
    for e in errs:
        print(f"error: {e}")
    if errs:
        print(
            f"{len(errs)} path(s) look like third-party material (#846). Unstage them "
            "(git rm --cached <path>) and keep the files in sources/."
        )
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
