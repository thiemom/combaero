"""The source register (validation/sources.yaml) stays valid.

Every local reference document has one entry: citation, content, what it was
read for, the verdict and where the repo uses it. The documents are not in
git, so CI checks the register's own consistency -- schema, statuses, and
that every record and code path it names exists. `scripts/sources.py check`
also lists local documents with no entry.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sources", REPO / "scripts" / "sources.py")
sources = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(sources)


def test_the_register_is_valid() -> None:
    entries = sources.load()
    assert entries, "the register is empty"
    assert sources.problems(entries) == []


def test_the_schema_check_catches_what_it_should() -> None:
    """Falsified: a bad status, an unknown topic, a missing record path and
    an 'implemented' entry naming no code are each reported."""
    good = sources.load()[0]
    bad = [
        {**good, "file": "sources/a.pdf", "status": "read"},
        {**good, "file": "sources/b.pdf", "topic": ["plasma"]},
        {**good, "file": "sources/c.pdf", "records": ["validation/nowhere.md"]},
        {**good, "file": "sources/d.pdf", "status": "implemented", "records": [], "code": []},
    ]
    found = sources.problems(bad)
    assert any("status 'read'" in p for p in found)
    assert any("topic must be" in p for p in found)
    assert any("does not exist" in p for p in found)
    assert any("needs a record or code reference" in p for p in found)


def test_a_file_outside_sources_is_reported() -> None:
    """Falsified (#846): a register entry still pointing into docs/ is an
    error -- copyrighted documents live in the ignored sources/ only."""
    good = sources.load()[0]
    assert good["file"].startswith("sources/")
    found = sources.problems([{**good, "file": "docs/junction/x.pdf"}])
    assert any("must be under sources/" in p for p in found)


def test_local_files_are_listed_under_sources(tmp_path: Path) -> None:
    """Every file below the sources dir is seen (any type, nested), as a
    sources/ path; only registrable types count as documents."""
    root = tmp_path / "sources"
    (root / "junction" / "crops").mkdir(parents=True)
    (root / "junction" / "paper.pdf").write_bytes(b"%PDF")
    (root / "junction" / "crops" / "fig.png").write_bytes(b"png")
    (root / ".gitignore").write_text("*\n")
    assert sources.local_files(root) == [
        "sources/junction/crops/fig.png",
        "sources/junction/paper.pdf",
    ]
    assert sources.local_documents(root) == ["sources/junction/paper.pdf"]


def test_the_wanted_list_is_valid() -> None:
    assert sources.wanted_problems(sources.load_wanted()) == []


def test_the_wanted_check_catches_what_it_should() -> None:
    """Falsified: an unknown priority, a held file, a non-boolean 'checked'
    and string issue numbers are each reported."""
    good = sources.load_wanted()[0]
    bad = [
        {**good, "priority": "urgent"},
        {**good, "file": "sources/x.pdf"},
        {**good, "checked": "no"},
        {**good, "issues": ["468"]},
    ]
    found = sources.wanted_problems(bad)
    assert any("priority 'urgent'" in p for p in found)
    assert any("belongs in 'documents'" in p for p in found)
    assert any("'checked' is true or false" in p for p in found)
    assert any("list of integers" in p for p in found)
