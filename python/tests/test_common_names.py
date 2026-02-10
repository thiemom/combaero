import pytest

import combaero as cb


def test_formula_to_name_basic() -> None:
    expected_keys = {
        "N2",
        "O2",
        "AR",
        "CO2",
        "H2O",
        "CH4",
        "C2H6",
        "C3H8",
        "IC4H10",
        "NC5H12",
        "NC6H14",
        "NC7H16",
        "CO",
        "H2",
    }

    mapping = cb.formula_to_name()
    keys = set(mapping.keys())
    assert expected_keys.issubset(keys)

    assert mapping["CH4"].lower().startswith("methane")
    assert "heptane" in mapping["NC7H16"].lower()
    assert "nitrogen" in mapping["N2"].lower()


def test_name_to_formula_basic() -> None:
    mapping = cb.name_to_formula()
    assert mapping["Methane"] == "CH4"
    assert mapping["Nitrogen"] == "N2"
    assert mapping["n-Heptane"] == "NC7H16"


def test_common_name_function() -> None:
    assert cb.common_name("CH4") == "Methane"
    assert cb.common_name("N2") == "Nitrogen"
    assert cb.common_name("H2O") == "Water"


def test_formula_function() -> None:
    assert cb.formula("Methane") == "CH4"
    assert cb.formula("Nitrogen") == "N2"
    assert cb.formula("Water") == "H2O"


def test_common_name_unknown_raises() -> None:
    with pytest.raises(IndexError):
        cb.common_name("UNKNOWN")


def test_formula_unknown_raises() -> None:
    with pytest.raises(IndexError):
        cb.formula("Unknown Gas")
