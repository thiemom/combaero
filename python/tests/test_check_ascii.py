"""Tests for scripts/check_ascii.py, the ascii-only pre-commit hook (#791).

Sources are synthetic strings written with \\u escapes, so this file stays
ASCII itself.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("check_ascii", ROOT / "scripts" / "check_ascii.py")
assert _spec is not None and _spec.loader is not None
check_ascii = importlib.util.module_from_spec(_spec)
sys.modules["check_ascii"] = check_ascii  # dataclasses resolve the module by name
_spec.loader.exec_module(check_ascii)

E_ACUTE = "\u00e9"
SUP2 = "\u00b2"
EM_DASH = "\u2014"
MU = "\u03bc"


def found(rel: str, text: str) -> list[tuple[int, int, str, str]]:
    """(line, col, char, category) of every offence, allowlist applied."""
    return [(o.line, o.col, o.char, o.category) for o in check_ascii.check_text(rel, text)]


def categories(rel: str, text: str) -> list[str]:
    return [o.category for o in check_ascii.check_text(rel, text)]


# ---------------------------------------------------------------- Python


def test_python_ascii_source_is_clean() -> None:
    assert found("a.py", 'x = "abc"  # fine\n') == []


def test_python_comment() -> None:
    assert found("a.py", f"x = 1  # 5 m{SUP2}\n") == [(1, 13, SUP2, "comment")]


def test_python_docstrings() -> None:
    src = (
        f'"""Module {EM_DASH} doc."""\n'
        "\n"
        "\n"
        "class A:\n"
        f'    """Class {EM_DASH} doc."""\n'
        "\n"
        "    def f(self) -> None:\n"
        "        # a comment before the docstring\n"
        f'        """Function {EM_DASH} doc."""\n'
        "        return None\n"
    )
    assert categories("a.py", src) == ["docstring"] * 3


def test_python_plain_literal() -> None:
    assert found("a.py", f'unit = "m{SUP2}/s"\n') == [(1, 10, SUP2, "string literal")]


def test_python_string_not_in_statement_position_is_literal() -> None:
    src = f'def f() -> str:\n    return "{MU}m"\n'
    assert categories("a.py", src) == ["string literal"]
    # A string that starts a statement but is not the whole statement.
    assert categories("a.py", f'"{MU}".join([])\n') == ["string literal"]


def test_python_fstring_middle() -> None:
    src = f'x = f"{{{{a}}}} {{y}} {E_ACUTE} {{z:>4}} {SUP2}"\n'
    assert categories("a.py", src) == ["string literal", "string literal"]
    # A non-ASCII name inside the braces is code, not literal.
    assert categories("a.py", f'x = f"a {{{E_ACUTE}}} b"\n') == ["identifier/code"]


def test_python_multiline_string_positions() -> None:
    src = f'x = """line one\nline {MU} two"""\n'
    assert found("a.py", src) == [(2, 6, MU, "string literal")]


def test_python_identifier() -> None:
    assert found("a.py", f"{E_ACUTE}t{E_ACUTE} = 1\n") == [
        (1, 1, E_ACUTE, "identifier/code"),
        (1, 3, E_ACUTE, "identifier/code"),
    ]


def test_python_untokenizable_falls_back_to_code() -> None:
    # A no-break space is not a valid token; the whole file degrades to code.
    assert categories("a.py", 'x = (\n"\u00a0"\n') == ["identifier/code"]


# ---------------------------------------------------------------- C++


def test_cpp_line_comment() -> None:
    assert found("a.cpp", f"int x = 1; // {MU}m\n") == [(1, 15, MU, "comment")]


def test_cpp_block_comment() -> None:
    src = f"/* first\n * second {EM_DASH} line */ int x;\n"
    assert found("a.h", src) == [(2, 11, EM_DASH, "comment")]


def test_cpp_normal_string_with_escapes() -> None:
    src = f'const char* s = "a\\"{SUP2}\\\\"; // {MU}\n'
    assert categories("a.cpp", src) == ["string literal", "comment"]


def test_cpp_raw_string_with_delimiter() -> None:
    # The raw string holds a quote, a ')"' and a '//' that must not end it.
    src = f'auto d = R"doc(a ")" // not a comment {EM_DASH}\nmore {MU})doc"; int y; // {SUP2}\n'
    assert categories("a.cpp", src) == ["string literal", "string literal", "comment"]


@pytest.mark.parametrize("prefix", ["u8R", "LR", "uR", "UR"])
def test_cpp_prefixed_raw_string(prefix: str) -> None:
    src = f'auto d = {prefix}"x({MU})x"; // {SUP2}\n'
    assert categories("a.cpp", src) == ["string literal", "comment"]


def test_cpp_char_literal() -> None:
    src = f"char q = '\"'; char c = '{E_ACUTE}'; // {MU}\n"
    assert categories("a.cpp", src) == ["string literal", "comment"]


def test_cpp_digit_separator_does_not_open_char_literal() -> None:
    # If 1'000 opened a char literal, the comment would be misread as a literal.
    src = f"long n = 1'000'000; double h = 0x1'F; // it's {MU}\nint {E_ACUTE} = 2;\n"
    assert found("a.cpp", src) == [(1, 47, MU, "comment"), (2, 5, E_ACUTE, "identifier/code")]


def test_cpp_apostrophe_in_comment_does_not_open_literal() -> None:
    src = f"// don't\nint x = 1; // {MU}\n"
    assert categories("a.cpp", src) == ["comment"]


def test_cpp_identifier_and_code() -> None:
    assert categories("a.hpp", f"int {MU}_value = 1;\n") == ["identifier/code"]


# ---------------------------------------------------------------- allowlist


def test_allowlisted_cpp_literal_passes_and_comment_fails() -> None:
    src = f'{{"P", "m{SUP2}"}}, // {EM_DASH}\nR"doc(x {MU})doc";\n'
    assert categories("include/units_data.h", src) == ["comment"]
    assert categories("python/combaero/_core.cpp", src) == ["comment"]


def test_allowlisted_python_literal_passes_but_docstring_and_comment_fail() -> None:
    src = f'"""Units {EM_DASH} table."""\nUNIT = "m{SUP2}"  # {MU}\n'
    assert categories("gui/backend/units.py", src) == ["docstring", "comment"]


def test_literal_fails_outside_allowlist() -> None:
    assert categories("include/other.h", f'const char* u = "m{SUP2}";\n') == ["string literal"]
    assert categories("gui/backend/other.py", f'UNIT = "m{SUP2}"\n') == ["string literal"]


def test_offence_format() -> None:
    (offence,) = check_ascii.check_text("a.py", f"# {MU}\n")
    assert offence.format() == f"a.py:1:3: non-ASCII U+03BC '{MU}' in comment"


def test_main_exit_codes(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    good = tmp_path / "good.py"
    good.write_text("x = 1\n", encoding="utf-8")
    bad = tmp_path / "bad.cpp"
    bad.write_text(f"// {MU}\n", encoding="utf-8")
    assert check_ascii.main([str(good)]) == 0
    assert check_ascii.main([str(good), str(bad)]) == 1
    out = capsys.readouterr().out
    assert "bad.cpp:1:4: non-ASCII U+03BC" in out
    assert "1 non-ASCII character(s) in 1 file(s)" in out
    assert check_ascii.main(["--stats", str(bad)]) == 0
