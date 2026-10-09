"""Changelog fragments: one file per change, folded into CHANGELOG.md at release.

Every PR used to insert its entry at the top of CHANGELOG.md's [Unreleased]
block, so any two open PRs conflicted there. Instead each change adds its own
file under changelog.d/ (see changelog.d/README.md):

    changelog.d/<id>.<category>.md

<id> is the issue or PR number (or a short slug); <category> is one of
CATEGORIES. The file holds Keep a Changelog bullets ("- ...", nested as
usual).

    uv run scripts/changelog.py check            # validate the fragments
    uv run scripts/changelog.py preview          # the would-be release section
    uv run scripts/changelog.py release 0.7.0    # fold into CHANGELOG.md

`release` merges the fragments with whatever is still written under
[Unreleased] in CHANGELOG.md (entries from before the fragments, or added by
hand), one section per category in Keep a Changelog order, dates the block,
moves the compare links on, and deletes the fragments.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHANGELOG = ROOT / "CHANGELOG.md"
FRAGMENTS = ROOT / "changelog.d"
REPO_URL = "https://github.com/thiemom/combaero"

# Keep a Changelog order, plus the "Documentation" section this changelog uses.
CATEGORIES = ("added", "changed", "deprecated", "removed", "fixed", "security", "documentation")
NAME_RE = re.compile(r"^(?P<id>[A-Za-z0-9][A-Za-z0-9_-]*)\.(?P<cat>[a-z]+)\.md$")
VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")


def fragments(directory: Path = FRAGMENTS) -> list[Path]:
    return sorted(p for p in directory.glob("*.md") if p.name != "README.md")


def check(directory: Path = FRAGMENTS) -> list[str]:
    """Problems with the fragments, one message each (empty: all good)."""
    problems = []
    for p in fragments(directory):
        m = NAME_RE.match(p.name)
        if not m:
            problems.append(f"{p.name}: name must be <id>.<category>.md")
            continue
        if m["cat"] not in CATEGORIES:
            problems.append(f"{p.name}: category '{m['cat']}' is not one of {', '.join(CATEGORIES)}")
        text = p.read_text(encoding="utf-8").strip()
        if not text:
            problems.append(f"{p.name}: empty")
        elif not text.startswith("- "):
            problems.append(f"{p.name}: must start with a '- ' bullet")
        elif re.search(r"^#", text, re.M):
            problems.append(f"{p.name}: no headings -- the category comes from the file name")
    return problems


def _sort_key(p: Path) -> tuple[int, int, str]:
    """Newest first: higher issue/PR numbers above lower; slugs after numbers."""
    ident = NAME_RE.match(p.name)["id"]
    return (0, -int(ident), "") if ident.isdigit() else (1, 0, ident)


def _split_unreleased(text: str) -> tuple[str, str, str]:
    """(before, [Unreleased] body, from the next release heading on)."""
    m = re.search(r"^## \[Unreleased\][^\n]*\n", text, re.M)
    if not m:
        raise SystemExit("CHANGELOG.md has no '## [Unreleased]' heading")
    nxt = re.search(r"^## \[", text[m.end() :], re.M)
    end = m.end() + nxt.start() if nxt else len(text)
    return text[: m.end()], text[m.end() : end], text[end:]


def _sections(body: str) -> dict[str, list[str]]:
    """[Unreleased] body -> {category: [block, ...]}, duplicate headings merged."""
    out: dict[str, list[str]] = {}
    current = None
    chunk: list[str] = []

    def flush() -> None:
        block = "\n".join(chunk).strip("\n")
        if current is not None and block.strip():
            out.setdefault(current, []).append(block)

    for line in body.split("\n"):
        h = re.match(r"^### (.+?)\s*$", line)
        if h:
            flush()
            current, chunk = h.group(1).strip().lower(), []
        else:
            chunk.append(line)
    flush()
    return out


def render(body: str, frags: list[Path]) -> str:
    """The merged section: CHANGELOG's own [Unreleased] entries per category,
    then the fragments (newest first), categories in Keep a Changelog order."""
    sections = _sections(body)
    by_cat: dict[str, list[Path]] = {}
    for p in frags:
        by_cat.setdefault(NAME_RE.match(p.name)["cat"], []).append(p)
    order = list(CATEGORIES) + [c for c in sections if c not in CATEGORIES]
    parts = []
    for cat in order:
        blocks = [p.read_text(encoding="utf-8").strip() for p in sorted(by_cat.get(cat, []), key=_sort_key)]
        blocks += sections.get(cat, [])
        if blocks:
            parts.append(f"### {cat.capitalize()}\n\n" + "\n\n".join(blocks))
    return "\n\n".join(parts)


def release(version: str, date: str, changelog: Path = CHANGELOG, directory: Path = FRAGMENTS) -> None:
    if not VERSION_RE.match(version):
        raise SystemExit(f"version must be X.Y.Z, got {version!r}")
    problems = check(directory)
    if problems:
        raise SystemExit("fragment problems:\n  " + "\n  ".join(problems))
    text = changelog.read_text(encoding="utf-8")
    if re.search(rf"^## \[{re.escape(version)}\]", text, re.M):
        raise SystemExit(f"CHANGELOG.md already has a [{version}] section")
    head, body, rest = _split_unreleased(text)
    frags = fragments(directory)
    merged = render(body, frags)
    if not merged.strip():
        raise SystemExit("nothing to release: no fragments and an empty [Unreleased]")
    new = head + "\n" + f"## [{version}] - {date}\n\n" + merged + "\n\n" + rest
    link = re.search(r"^\[Unreleased\]: \S+/compare/v(\d+\.\d+\.\d+)\.\.\.HEAD$", new, re.M)
    if link:
        prev = link.group(1)
        new = new.replace(
            link.group(0),
            f"[Unreleased]: {REPO_URL}/compare/v{version}...HEAD\n"
            f"[{version}]: {REPO_URL}/compare/v{prev}...v{version}",
        )
    changelog.write_text(new, encoding="utf-8")
    for p in frags:
        p.unlink()
    print(f"CHANGELOG.md: [{version}] - {date}, {len(frags)} fragment(s) folded in and removed.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="validate changelog.d/ fragments")
    sub.add_parser("preview", help="print the would-be release section")
    rel = sub.add_parser("release", help="fold the fragments into CHANGELOG.md")
    rel.add_argument("version")
    rel.add_argument("--date", default=dt.date.today().isoformat())
    args = ap.parse_args(argv)
    if args.cmd == "check":
        problems = check()
        for msg in problems:
            print(f"changelog.d/{msg}", file=sys.stderr)
        return 1 if problems else 0
    if args.cmd == "preview":
        _, body, _ = _split_unreleased(CHANGELOG.read_text(encoding="utf-8"))
        print(render(body, fragments()))
        return 0
    release(args.version, args.date)
    return 0


if __name__ == "__main__":
    sys.exit(main())
