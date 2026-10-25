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
        {**good, "file": "docs/a.pdf", "status": "read"},
        {**good, "file": "docs/b.pdf", "topic": ["plasma"]},
        {**good, "file": "docs/c.pdf", "records": ["validation/nowhere.md"]},
        {**good, "file": "docs/d.pdf", "status": "implemented", "records": [], "code": []},
    ]
    found = sources.problems(bad)
    assert any("status 'read'" in p for p in found)
    assert any("topic must be" in p for p in found)
    assert any("does not exist" in p for p in found)
    assert any("needs a record or code reference" in p for p in found)
