from __future__ import annotations

import csv
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
    csv_path = tmp_path / "species_data.csv"

    _write_minimal_yaml_with_he(yaml_path)

    extract_species_data(
        yaml_file=str(yaml_path),
        output_file=str(csv_path),
        selected_species=["H2", "HE"],
    )

    rows: list[dict[str, str]] = []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        rows.extend(reader)

    # We expect two rows: H2 and HE
    assert len(rows) == 2

    species_names = {row["species_name"] for row in rows}
    assert "H2" in species_names
    assert "HE" in species_names

    # Check that the molar mass for helium is computed from the atomic mass table
    he_row = next(row for row in rows if row["species_name"] == "HE")
    # 4.002602 is the value we added for He/HE; allow small float formatting differences
    assert abs(float(he_row["molar_mass"]) - 4.002602) < 1e-6

    # Basic sanity for thermo and transport fields
    for row in rows:
        assert row["T_low"] == "200.0"
        assert row["T_mid"] == "1000.0"
        assert row["T_high"] == "3500.0"
        assert row["geometry"] in {"linear", "atom"}
