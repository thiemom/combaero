from __future__ import annotations

from pathlib import Path

from extract_thermo_transport_data import extract_species_data


def _write_minimal_yaml_with_he(path: Path) -> None:
    content = """
phases:
  - name: gas
    thermo: ideal-gas
    species: [H2, HE]
    state:
      T: 300.0
      P: 101325.0

species:
  - name: H2
    composition: {H: 2}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [1, 2, 3, 4, 5, 6, 7]
        - [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1]
    transport:
      geometry: linear
      well-depth: 38.0
      diameter: 2.92
      polarizability: 0.79
      rotational-relaxation: 280.0

  - name: HE
    composition: {He: 1}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [1, 2, 3, 4, 5, 6, 7]
        - [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1]
    transport:
      geometry: atom
      well-depth: 10.0
      diameter: 2.0
      polarizability: 0.0
"""
    path.write_text(content, encoding="utf-8")


def test_extract_species_data_includes_helium(tmp_path: Path) -> None:
    yaml_path = tmp_path / "minimal_with_he.yaml"

    _write_minimal_yaml_with_he(yaml_path)

    species_data, warnings = extract_species_data(
        yaml_file=yaml_path,
        selected_species=["H2", "HE"],
    )

    # We expect two species: H2 and HE
    assert len(species_data) == 2

    species_names = {sp["species_name"] for sp in species_data}
    assert "H2" in species_names
    assert "HE" in species_names

    # Check that the molar mass for helium is computed from the atomic mass table
    he_entry = next(sp for sp in species_data if sp["species_name"] == "HE")
    # 4.0026 is the value in ATOMIC_MASSES for He
    assert abs(he_entry["molar_mass"] - 4.0026) < 1e-3

    # Basic sanity for thermo and transport fields
    for sp in species_data:
        assert sp["thermo"]["T_low"] == 200.0
        assert sp["thermo"]["T_mid"] == 1000.0
        assert sp["thermo"]["T_high"] == 3500.0
        assert sp["transport"]["geometry"] in {"linear", "atom"}

    # Check rotational_relaxation is extracted
    h2_entry = next(sp for sp in species_data if sp["species_name"] == "H2")
    assert h2_entry["transport"]["rotational_relaxation"] == 280.0
