"""scripts/changelog.py: changelog fragments folded into CHANGELOG.md."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("changelog", ROOT / "scripts" / "changelog.py")
changelog = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(changelog)

BASE = """# Changelog

Intro.

## [Unreleased]

### Fixed

- **Old fix.** Written before the fragments.

### Changed

- **Old change.**

### Changed

- **A second Changed section.**

## [0.6.0] - 2026-09-09

### Changed

- **Released.**

[Unreleased]: https://github.com/thiemom/combaero/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/thiemom/combaero/compare/v0.5.0...v0.6.0
"""


@pytest.fixture
def repo(tmp_path: Path) -> tuple[Path, Path]:
    cl = tmp_path / "CHANGELOG.md"
    cl.write_text(BASE, encoding="utf-8")
    d = tmp_path / "changelog.d"
    d.mkdir()
    (d / "README.md").write_text("# not a fragment\n", encoding="utf-8")
    (d / "496.fixed.md").write_text("- **Newer fix (#496).**\n  - **Now:** detail.\n")
    (d / "463.fixed.md").write_text("- **Older fix (#463).**\n")
    (d / "tooling.added.md").write_text("- **A new script.**\n")
    return cl, d


def test_release_folds_fragments_into_a_dated_block(repo: tuple[Path, Path]) -> None:
    cl, d = repo
    changelog.release("0.7.0", "2026-10-09", changelog=cl, directory=d)
    text = cl.read_text(encoding="utf-8")
    unreleased, rest = text.split("## [0.7.0] - 2026-10-09\n", 1)
    assert unreleased.rstrip().endswith("## [Unreleased]")  # left empty
    block = rest.split("## [0.6.0]")[0]
    # Keep a Changelog order; one section per category; duplicates merged.
    heads = [line for line in block.splitlines() if line.startswith("### ")]
    assert heads == ["### Added", "### Changed", "### Fixed"]
    # Newest fragment first, then what was already under [Unreleased].
    fixed = block.split("### Fixed")[1]
    assert fixed.index("#496") < fixed.index("#463") < fixed.index("Old fix")
    assert "Old change" in block and "A second Changed section" in block
    # Compare links moved on.
    assert "[Unreleased]: https://github.com/thiemom/combaero/compare/v0.7.0...HEAD" in text
    assert "[0.7.0]: https://github.com/thiemom/combaero/compare/v0.6.0...v0.7.0" in text
    # Fragments consumed, the README kept.
    assert [p.name for p in d.iterdir()] == ["README.md"]


def test_check_reports_malformed_fragments(tmp_path: Path) -> None:
    (tmp_path / "1.fixed.md").write_text("- fine\n")
    (tmp_path / "2.bugfix.md").write_text("- wrong category\n")
    (tmp_path / "notes.md").write_text("- no category\n")
    (tmp_path / "3.added.md").write_text("### Added\n- heading inside\n")
    (tmp_path / "4.added.md").write_text("   \n")
    problems = changelog.check(tmp_path)
    assert len(problems) == 4
    assert any("2.bugfix.md" in p for p in problems)
    assert not any(p.startswith("1.fixed.md") for p in problems)


def test_release_refuses_an_existing_version_and_a_bad_number(repo: tuple[Path, Path]) -> None:
    cl, d = repo
    with pytest.raises(SystemExit):
        changelog.release("0.6.0", "2026-10-09", changelog=cl, directory=d)
    with pytest.raises(SystemExit):
        changelog.release("v0.7", "2026-10-09", changelog=cl, directory=d)
    assert cl.read_text(encoding="utf-8") == BASE  # untouched


def test_the_repositorys_fragments_are_valid() -> None:
    assert changelog.check() == []
