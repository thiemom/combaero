"""The fidelity count: #389's question 1, which had no column.

What this guards is that the registry describes reality -- every declared
check names a test that exists, every set that scores is accounted for,
and the distinction between "no check" and "not looked at" stays visible.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from validation.cooling import fidelity  # noqa: E402
from validation.cooling.schema import load_dataset  # noqa: E402


@pytest.fixture(scope="module")
def dataset():
    return load_dataset()


def test_every_declared_check_exists() -> None:
    """A registry that names a test which does not exist is worse than
    none: it reports evidence nobody can follow.

    Python nodes are resolved by collection. ctest names are checked
    against the C++ sources, because building to collect them would make
    this test depend on a build.
    """
    python_nodes = [c.pinned_by for c in fidelity.CHECKS if not c.pinned_by.startswith("ctest:")]
    assert python_nodes, "no python-pinned checks; the resolution below proves nothing"

    collected = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", *python_nodes],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert collected.returncode == 0, (
        "a declared fidelity check names a test that does not collect:\n"
        + collected.stdout[-2000:]
        + collected.stderr[-2000:]
    )

    cpp = "\n".join(p.read_text() for p in (REPO / "tests").glob("*.cpp"))
    for c in fidelity.CHECKS:
        if not c.pinned_by.startswith("ctest:"):
            continue
        suite, name = c.pinned_by[len("ctest:") :].split(".", 1)
        assert f"TEST({suite}, {name})" in cpp, (
            f"{c.pinned_by} is declared but no such gtest exists"
        )


def test_every_declared_check_describes_a_printed_quantity() -> None:
    """The definition this registry rests on.

    A fidelity check compares against something the SOURCE PRINTS -- a
    constant, a table entry, an equation the author evaluates. Scoring
    against a digitised measurement is the accuracy basis the scorecard
    already reports and cannot separate transcription from model error.
    """
    for c in fidelity.CHECKS:
        assert c.printed and c.location and c.agreement, c
        # The location must name where on the page, not just the paper.
        assert any(
            word in c.location for word in ("Eq.", "Eqs.", "Table", "Fig.", "figure", "branch")
        ), f"{c.pinned_by}: location '{c.location}' does not say where"


def test_the_count_accounts_for_every_scoring_set(dataset) -> None:
    """A set with no check must be NAMED, not silently absent.

    "No printed-quantity check" and "nobody looked" read identically in a
    count of zero, and the first is a fact while the second is a gap.
    """
    scoring = {s.scores for s in dataset if s.scores}
    counts = fidelity.count_by_set(dataset)
    empty = set(fidelity.sets_without_checks(dataset))

    # Every scoring set is either counted or named as having none.
    unaccounted = scoring - set(counts) - empty
    assert not unaccounted, sorted(unaccounted)

    # And never both -- the contradiction a hard-coded list produced:
    # han_park_1988_angled reported one machine-counted check while being
    # listed as having none, in adjacent lines of the same report.
    assert not (set(counts) & empty), sorted(set(counts) & empty)
    assert "han_park_1988_angled" in counts
    assert "han_park_1988_angled" not in empty


def test_the_counts_are_what_the_registry_says(dataset) -> None:
    """The headline numbers, pinned so a silent drop is visible."""
    counts = fidelity.count_by_set(dataset)
    assert counts["andrews_1986_effusion_internal"] == 7
    assert counts["baldauf_2002_sellers"] == 2
    assert counts["han_1988_orthogonal"] >= 1
    assert sum(1 for c in fidelity.CHECKS) == 10


def test_a_check_whose_result_is_a_contradiction_still_counts() -> None:
    """The one that would be perverse to treat as a defect.

    Baldauf's Eq. (31) disagrees with the paper's own Table 4 by 36%, and
    the equation is implemented as printed. That is fidelity evidence of
    the most useful kind -- it says exactly where the source is
    self-inconsistent -- so it belongs in the count, not outside it.
    """
    contradiction = next(
        c
        for c in fidelity.CHECKS
        if "contradicts" in c.pinned_by.lower() or "DISAGREES" in c.agreement
    )
    assert contradiction.correlation_set == "baldauf_2002_sellers"
    assert "AS PRINTED" in contradiction.agreement
    assert contradiction in fidelity.CHECKS


def test_the_machine_counted_part_is_real_but_small(dataset) -> None:
    """#389 expected this to be mostly aggregation. Measured, it is not.

    `verify.py` finds 9 `printed-curve`/`printed-exponent` findings across
    the whole dataset and only one lands on a scored set -- the rest sit
    on `scores: null` correlation curves belonging to no set. Recorded so
    the registry's existence is justified by a number rather than an
    assertion.
    """
    machine = fidelity.machine_counted(dataset)
    assert sum(machine.values()) <= 3, (
        f"the machine-countable part grew to {machine}; if it is now "
        "substantial, reconsider how much of this registry needs declaring"
    )
    assert len(fidelity.CHECKS) > 3 * max(sum(machine.values()), 1)


def test_render_states_both_what_is_checked_and_what_is_not(dataset) -> None:
    text = "\n".join(fidelity.render(dataset))
    assert "EVIDENCE, not a score" in text
    assert "No printed-quantity check at all" in text
    for name in fidelity.sets_without_checks(dataset):
        assert name in text
    # And a set WITH a check must not appear in that list.
    assert "han_park_1988_angled" not in text.split("No printed-quantity")[1]
    assert "0.26983 against 0.27" in text
