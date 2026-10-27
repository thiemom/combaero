"""Reject non-ASCII characters in C++ and Python sources (#791).

Non-ASCII is rejected everywhere: comments, docstrings, identifiers and
code. The one exception is a string literal in a file on ALLOWLIST, the
files that carry display units ("m^2/s" written with a real superscript
reads better in review than a \\u escape). Comments and docstrings in an
allowlisted file must still be ASCII.

    uv run scripts/check_ascii.py FILE...      # what pre-commit runs
    uv run scripts/check_ascii.py --all        # every tracked source
    uv run scripts/check_ascii.py --all --stats   # counts, never fails

Python is classified with `tokenize`: STRING and FSTRING_MIDDLE tokens are
literals, except a bare string statement (a STRING in statement position,
i.e. a module, class or function docstring), which is a docstring. C++ is
classified by a small lexer that knows line and block comments, ordinary
and raw string literals (with the u8/u/U/L prefixes), character literals
and digit separators.
"""

from __future__ import annotations

import argparse
import io
import subprocess
import sys
import tokenize
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

# Files whose string literals may hold non-ASCII (display units). Paths are
# relative to the repository root, POSIX style.
ALLOWLIST: frozenset[str] = frozenset(
    {
        "include/units_data.h",
        "python/combaero/_core.cpp",
        "gui/backend/units.py",
    }
)

PY_SUFFIXES: frozenset[str] = frozenset({".py"})
CPP_SUFFIXES: frozenset[str] = frozenset({".h", ".hpp", ".cpp", ".cc"})
GIT_PATTERNS: tuple[str, ...] = ("*.py", "*.h", "*.hpp", "*.cpp", "*.cc")

COMMENT = "comment"
DOCSTRING = "docstring"
CODE = "identifier/code"
LITERAL = "string literal"
CATEGORIES: tuple[str, ...] = (COMMENT, DOCSTRING, CODE, LITERAL)


@dataclass(frozen=True)
class Offence:
    path: str
    line: int  # 1-based
    col: int  # 1-based, in characters
    char: str
    category: str

    def format(self) -> str:
        shown = self.char if self.char.isprintable() else "?"
        return (
            f"{self.path}:{self.line}:{self.col}: non-ASCII "
            f"U+{ord(self.char):04X} '{shown}' in {self.category}"
        )


# ---------------------------------------------------------------- Python


def _lines(text: str) -> list[str]:
    """Physical lines as tokenize sees them: split on \\n only, not on U+2028."""
    return io.StringIO(text).readlines()


def _paint(
    spans: dict[int, list[tuple[int, int, str]]],
    start: tuple[int, int],
    end: tuple[int, int],
    lines: list[str],
    category: str,
) -> None:
    """Record that [start, end) (tokenize coordinates) has `category`."""
    (r0, c0), (r1, c1) = start, end
    for row in range(r0, r1 + 1):
        lo = c0 if row == r0 else 0
        hi = c1 if row == r1 else len(lines[row - 1]) if row - 1 < len(lines) else 0
        spans[row].append((lo, hi, category))


_SKIP = frozenset(
    {
        tokenize.NL,
        tokenize.COMMENT,
        tokenize.INDENT,
        tokenize.DEDENT,
        tokenize.ENCODING,
    }
)
_STATEMENT_START = frozenset({tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT})


def _python_spans(text: str) -> dict[int, list[tuple[int, int, str]]] | None:
    """Category spans per line, later entries winning; None if untokenizable."""
    lines = _lines(text)
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (tokenize.TokenError, SyntaxError):
        return None

    fstring_start = getattr(tokenize, "FSTRING_START", -1)
    fstring_middle = getattr(tokenize, "FSTRING_MIDDLE", -2)
    fstring_end = getattr(tokenize, "FSTRING_END", -3)

    # Indices of significant tokens, to look one back and ahead of a STRING.
    sig = [i for i, t in enumerate(toks) if t.type not in _SKIP]
    pos_in_sig = {tok_index: k for k, tok_index in enumerate(sig)}

    def prev_sig_type(i: int) -> int | None:
        # INDENT/DEDENT are skipped as "insignificant" for adjacency but do
        # mark a statement start, so look at the raw stream for them.
        j = i - 1
        while j >= 0 and toks[j].type in (tokenize.NL, tokenize.COMMENT, tokenize.ENCODING):
            j -= 1
        return toks[j].type if j >= 0 else None

    def is_docstring(i: int) -> bool:
        before = prev_sig_type(i)
        if before is not None and before not in _STATEMENT_START:
            return False
        # The statement must be nothing but (implicitly concatenated) strings.
        k = pos_in_sig[i] + 1
        while k < len(sig) and toks[sig[k]].type == tokenize.STRING:
            k += 1
        if k >= len(sig):
            return True
        nxt = toks[sig[k]]
        return nxt.type in (tokenize.NEWLINE, tokenize.ENDMARKER) or (
            nxt.type == tokenize.OP and nxt.string == ";"
        )

    spans: dict[int, list[tuple[int, int, str]]] = defaultdict(list)
    in_docstring_run = False
    for i, tok in enumerate(toks):
        if tok.type == tokenize.COMMENT:
            _paint(spans, tok.start, tok.end, lines, COMMENT)
        elif tok.type == tokenize.STRING:
            category = DOCSTRING if in_docstring_run or is_docstring(i) else LITERAL
            _paint(spans, tok.start, tok.end, lines, category)
            in_docstring_run = category == DOCSTRING
            continue
        elif tok.type in (fstring_start, fstring_middle, fstring_end):
            # Token positions are exact in 3.12: only the doubled brace of an
            # escaped {{ or }} falls between two FSTRING_MIDDLE tokens, and it
            # is ASCII. Expressions inside the braces are ordinary tokens.
            _paint(spans, tok.start, tok.end, lines, LITERAL)
        elif (
            tok.type not in (tokenize.NL, tokenize.NEWLINE, tokenize.ENDMARKER)
            and tok.start != tok.end
        ):
            _paint(spans, tok.start, tok.end, lines, CODE)
        if tok.type not in _SKIP:
            in_docstring_run = False
    return spans


def classify_python(text: str) -> list[tuple[int, int, str, str]]:
    """(line, col0, char, category) for every non-ASCII char in `text`."""
    spans = _python_spans(text)
    out: list[tuple[int, int, str, str]] = []
    for lineno, line in enumerate(_lines(text), start=1):
        if line.isascii():
            continue
        for col, ch in enumerate(line):
            if ord(ch) < 128:
                continue
            category = CODE
            if spans is not None:
                for lo, hi, cat in reversed(spans.get(lineno, [])):
                    if lo <= col < hi:
                        category = cat
                        break
            out.append((lineno, col, ch, category))
    return out


# ---------------------------------------------------------------- C++

_RAW_PREFIXES = frozenset({"R", "u8R", "uR", "UR", "LR"})
_PLAIN_PREFIXES = frozenset({"u8", "u", "U", "L"})


def _is_ident_char(ch: str) -> bool:
    return ch.isalnum() or ch == "_" or ord(ch) > 127


def classify_cpp(text: str) -> list[tuple[int, int, str, str]]:
    """(line, col0, char, category) for every non-ASCII char in `text`."""
    n = len(text)
    cats = [CODE] * n  # category per character index; only non-ASCII read

    def mark(lo: int, hi: int, category: str) -> None:
        for k in range(lo, min(hi, n)):
            if ord(text[k]) > 127:
                cats[k] = category

    def quoted_end(i: int, quote: str) -> int:
        """Index past the closing quote of a literal opening at text[i]."""
        j = i + 1
        while j < n:
            ch = text[j]
            if ch == "\\":
                j += 2
                continue
            if ch == quote:
                return j + 1
            if ch == "\n":  # unterminated; stop at the line end
                return j
            j += 1
        return n

    i = 0
    while i < n:
        ch = text[i]
        if ch == "/" and text.startswith("//", i):
            j = i
            while j < n and text[j] != "\n":
                # A backslash-newline continues a line comment.
                if text[j] == "\\" and text.startswith("\n", j + 1):
                    j += 2
                    continue
                j += 1
            mark(i, j, COMMENT)
            i = j
        elif ch == "/" and text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            mark(i, j, COMMENT)
            i = j
        elif ch == '"':
            j = quoted_end(i, '"')
            mark(i, j, LITERAL)
            i = j
        elif ch == "'":
            j = quoted_end(i, "'")
            mark(i, j, LITERAL)
            i = j
        elif ch.isdigit() or (ch == "." and i + 1 < n and text[i + 1].isdigit()):
            # pp-number: digits, letters, '.', exponent signs and the C++14
            # digit separator ' (only between two alphanumerics).
            j = i + 1
            while j < n:
                c = text[j]
                exponent_sign = c in "+-" and text[j - 1] in "eEpP"
                separator = c == "'" and j + 1 < n and text[j + 1].isalnum()
                if not (c.isalnum() or c in "._" or exponent_sign or separator):
                    break
                j += 1
            i = j
        elif _is_ident_char(ch):
            j = i + 1
            while j < n and _is_ident_char(text[j]):
                j += 1
            word = text[i:j]
            if word in _RAW_PREFIXES and j < n and text[j] == '"':
                open_paren = text.find("(", j + 1)
                if open_paren < 0:
                    i = j
                    continue
                delim = text[j + 1 : open_paren]
                close = text.find(")" + delim + '"', open_paren + 1)
                end = n if close < 0 else close + len(delim) + 2
                mark(j, end, LITERAL)
                i = end
            elif word in _PLAIN_PREFIXES and j < n and text[j] in "\"'":
                end = quoted_end(j, text[j])
                mark(j, end, LITERAL)
                i = end
            else:
                i = j  # identifier; its non-ASCII stays CODE
        else:
            i += 1

    out: list[tuple[int, int, str, str]] = []
    line, col = 1, 0
    for k, c in enumerate(text):
        if c == "\n":
            line, col = line + 1, 0
            continue
        if ord(c) > 127:
            out.append((line, col, c, cats[k]))
        col += 1
    return out


# ---------------------------------------------------------------- driver


def repo_root() -> Path:
    out = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=True,
    )
    return Path(out.stdout.strip())


def rel_path(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def tracked_sources(root: Path) -> list[Path]:
    out = subprocess.run(
        ["git", "ls-files", "-z", "--", *GIT_PATTERNS],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    )
    return sorted(root / p for p in out.stdout.split("\0") if p)


def scan_text(rel: str, text: str) -> list[Offence]:
    """Every non-ASCII char in `text`, classified; the allowlist is not applied."""
    if text.isascii():
        return []
    suffix = Path(rel).suffix
    if suffix in PY_SUFFIXES:
        found = classify_python(text)
    elif suffix in CPP_SUFFIXES:
        found = classify_cpp(text)
    else:
        found = [
            (ln, col, ch, CODE)
            for ln, line in enumerate(_lines(text), start=1)
            for col, ch in enumerate(line)
            if ord(ch) > 127
        ]
    return [Offence(rel, ln, col + 1, ch, cat) for ln, col, ch, cat in found]


def is_allowed(offence: Offence) -> bool:
    return offence.category == LITERAL and offence.path in ALLOWLIST


def check_text(rel: str, text: str) -> list[Offence]:
    """The offences in `text`, with the allowlist applied."""
    return [o for o in scan_text(rel, text) if not is_allowed(o)]


def read_source(path: Path) -> str:
    return path.read_bytes().decode("utf-8", errors="replace")


def print_stats(found: list[Offence]) -> None:
    def label(o: Offence) -> str:
        return f"{o.category} (allowlisted)" if is_allowed(o) else o.category

    def kind(o: Offence) -> str:
        return "py" if Path(o.path).suffix in PY_SUFFIXES else "C++"

    chars: Counter[tuple[str, str]] = Counter()
    files: dict[tuple[str, str], set[str]] = defaultdict(set)
    for o in found:
        key = (kind(o), label(o))
        chars[key] += 1
        files[key].add(o.path)
    print("Non-ASCII characters per language and category (chars / files):")
    for key in sorted(chars):
        print(f"  {key[0]:<4} {key[1]:<30} {chars[key]:>6} / {len(files[key]):>3}")

    top_chars: Counter[str] = Counter()
    top_files: dict[str, set[str]] = defaultdict(set)
    failing: set[str] = set()
    for o in found:
        top = o.path.split("/", 1)[0] if "/" in o.path else "."
        top_chars[top] += 1
        top_files[top].add(o.path)
        if not is_allowed(o):
            failing.add(o.path)
    print("Per top-level directory (chars / files):")
    for top in sorted(top_chars):
        print(f"  {top:<28} {top_chars[top]:>6} / {len(top_files[top]):>3}")
    names = Counter(unicodedata.name(o.char, f"U+{ord(o.char):04X}") for o in found)
    print("Most common characters:")
    for name, count in names.most_common(15):
        print(f"  {count:>6}  {name}")
    print(f"Files that fail the check: {len(failing)}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="*", type=Path, help="files to check")
    parser.add_argument(
        "--all", action="store_true", help="check every git-tracked C++/Python source"
    )
    parser.add_argument(
        "--stats", action="store_true", help="print counts per category and directory; never fail"
    )
    args = parser.parse_args(argv)

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(errors="backslashreplace")

    root = repo_root()
    paths: list[Path] = tracked_sources(root) if args.all else list(args.files)
    suffixes = PY_SUFFIXES | CPP_SUFFIXES

    found: list[Offence] = []
    for path in paths:
        if path.suffix not in suffixes or not path.is_file():
            continue
        found.extend(scan_text(rel_path(path, root), read_source(path)))

    if args.stats:
        print_stats(found)
        return 0

    offences = [o for o in found if not is_allowed(o)]
    for o in offences:
        print(o.format())
    if offences:
        n_files = len({o.path for o in offences})
        print(f"{len(offences)} non-ASCII character(s) in {n_files} file(s); see #791.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
