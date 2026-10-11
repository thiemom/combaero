"""Tests for extract_species_data.py's source merging (#485).

Flagged transport (a mechanism's own 'WARNING ... theoret trans') must yield
to a later source's unflagged parameters; thermo keeps first-wins.
"""

from __future__ import annotations

from pathlib import Path

from extract_species_data import merge_yaml_sources, transport_is_flagged


def _entry(eps: float, source: str, flagged: bool, thermo: str) -> dict:
    return {
        "name": "O2",
        "thermo_nasa7": {"tag": thermo},
        "transport": {"well_depth": eps, "diameter": 3.0},
        "transport_source": source,
        "transport_flagged": flagged,
    }


def test_the_nuig_note_is_recognised_as_flagged() -> None:
    assert transport_is_flagged(
        {"note": "\\AUTHOR: WARNING !\\REF: WARNING !\\COMMENT: theoret trans"}
    )
    assert not transport_is_flagged({"note": ""})
    assert not transport_is_flagged({})


def test_flagged_transport_yields_to_a_later_unflagged_source() -> None:
    merged = merge_yaml_sources(
        [
            (Path("NUIG.yaml"), {"O2": _entry(676.424, "NUIG.yaml", True, "first")}),
            (Path("JetSurf2.yaml"), {"O2": _entry(107.4, "JetSurf2.yaml", False, "second")}),
        ]
    )
    o2 = merged["O2"]
    assert o2["transport"]["well_depth"] == 107.4
    assert o2["transport_source"] == "JetSurf2.yaml"
    assert not o2["transport_flagged"]
    # Thermo is untouched: first source still wins there.
    assert o2["thermo_nasa7"]["tag"] == "first"


def test_unflagged_first_source_keeps_its_transport() -> None:
    merged = merge_yaml_sources(
        [
            (Path("A.yaml"), {"O2": _entry(107.4, "A.yaml", False, "first")}),
            (Path("B.yaml"), {"O2": _entry(110.0, "B.yaml", False, "second")}),
        ]
    )
    assert merged["O2"]["transport"]["well_depth"] == 107.4


def test_a_flagged_entry_with_no_alternative_is_kept() -> None:
    """NH3 exists only in NUIGMech; flagged or not, it is what there is."""
    merged = merge_yaml_sources(
        [(Path("NUIG.yaml"), {"NH3": _entry(481.0, "NUIG.yaml", True, "x")})]
    )
    assert merged["NH3"]["transport"]["well_depth"] == 481.0
