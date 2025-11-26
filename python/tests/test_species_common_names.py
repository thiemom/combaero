from combaero import species_common_names


def test_species_common_names_basic() -> None:
    expected_keys = {
        "N2", "O2", "AR", "CO2", "H2O", "CH4", "C2H6",
        "C3H8", "IC4H10", "NC5H12", "NC6H14", "NC7H16",
        "CO", "H2",
    }

    mapping = species_common_names()
    keys = set(mapping.keys())
    assert expected_keys.issubset(keys)

    assert mapping["CH4"].lower().startswith("methane")
    assert "heptane" in mapping["NC7H16"].lower()
    assert "nitrogen" in mapping["N2"].lower()
