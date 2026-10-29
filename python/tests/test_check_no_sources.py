"""scripts/check_no_sources.py keeps copyrighted material out of commits (#846).

Falsified both ways: each forbidden path is reported, and each path we do
commit (our plots, our docs, the sources/ lock file) passes.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "check_no_sources", REPO / "scripts" / "check_no_sources.py"
)
check = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(check)


@pytest.mark.parametrize(
    "path",
    [
        "sources/new/paper.pdf",
        "sources/junction/notes.md",
        "sources/materials/table.csv",
        "docs/paper.pdf",
        "docs/heat_transfer/film/scan.PDF",
        "docs/somenew/fig.png",
        "docs/somenew/data.xlsx",
        "validation/cooling/data/x/supp.zip",
        "validation/cooling/plots/x/handbook.pdf",
        "python/tests/data/page.tif",
        "release.tar.gz",
    ],
)
def test_forbidden_paths_fail(path: str) -> None:
    assert check.problem(path) is not None
    assert check.main([path]) == 1


@pytest.mark.parametrize(
    "path",
    [
        "validation/cooling/plots/han2012/han2012_fig4_47.png",
        ".github/assets/combaero_banner.png",
        "gui/frontend/src/assets/hero.png",
        "docs/VALIDATION_POLICY.md",
        "docs/new_topic/notes.md",
        "sources/.gitignore",
        "validation/cooling/data/han2012/fig.csv",
        "gui/frontend/public/icons.svg",
    ],
)
def test_allowed_paths_pass(path: str) -> None:
    assert check.problem(path) is None
    assert check.main([path]) == 0


def test_every_tracked_file_passes() -> None:
    assert check.main(["--all"]) == 0
