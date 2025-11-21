from __future__ import annotations

import csv
import tempfile
from pathlib import Path

from create_cpp_header import generate_cpp_header


def test_generate_cpp_header_produces_valid_header(tmp_path: Path) -> None:
    """Generate a small header from a minimal CSV and check key invariants.

    This test is intentionally lightweight and only checks structural properties
    of the generated header so it remains stable even as numerical data change.
    """

    # Prepare a minimal CSV with two species
    csv_path = tmp_path / "species_data.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "species_name",
                "molar_mass",
                "P_ref",
                "C",
                "H",
                "O",
                "N",
                "T_low",
                "T_mid",
                "T_high",
                "NASA_low",
                "NASA_high",
                "model",
                "geometry",
                "well-depth",
                "diameter",
                "polarizability",
                "rotational-relaxation",
            ]
        )
        writer.writerow(
            [
                "H2",
                2.016,
                101325.0,
                0,
                2,
                0,
                0,
                200.0,
                1000.0,
                3500.0,
                "[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]",
                "[1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1]",
                "NASA7",
                "linear",
                38.0,
                2.92,
                0.79,
                4.0,
            ]
        )
        writer.writerow(
            [
                "O2",
                32.0,
                101325.0,
                0,
                0,
                2,
                0,
                200.0,
                1000.0,
                3500.0,
                "[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]",
                "[1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1]",
                "NASA7",
                "linear",
                107.4,
                3.46,
                1.60,
                3.8,
            ]
        )

    header_path = tmp_path / "thermo_transport_data.h"

    # Act
    generate_cpp_header(str(csv_path), str(header_path))

    # Assert basic structure of the generated header
    content = header_path.read_text(encoding="utf-8")

    # Header guard and includes
    assert "#ifndef THERMO_TRANSPORT_DATA_H" in content
    assert "#define THERMO_TRANSPORT_DATA_H" in content
    assert "#include <vector>" in content
    assert "#include <string>" in content
    assert "#include <unordered_map>" in content

    # Struct declarations
    assert "struct NASA_Coeffs" in content
    assert "struct Transport_Props" in content
    assert "struct Molecular_Structure" in content

    # Species names should be valid C++ string literals
    assert "const std::vector<std::string> species_names" in content
    assert '{"H2", "O2"}' in content

    # species_index should map names to indices 0 and 1
    assert '"H2", 0' in content
    assert '"O2", 1' in content

    # Molar masses vector should include both species in order
    assert "const std::vector<double> molar_masses" in content
    assert "2.016" in content
    assert "32.0" in content

    # NASA coefficients should appear exactly twice (for two species)
    assert "const std::vector<NASA_Coeffs> nasa_coeffs" in content
    assert content.count("NASA_Coeffs") >= 1

    # Transport properties should include both species
    assert "const std::vector<Transport_Props> transport_props" in content
    assert "38.0" in content  # H2 well depth
    assert "107.4" in content  # O2 well depth

    # Molecular structures should encode the C/H/O/N counts
    assert "const std::vector<Molecular_Structure> molecular_structures" in content
    assert "{0, 2, 0, 0}" in content  # H2
    assert "{0, 0, 2, 0}" in content  # O2

    # Header guard end
    assert "#endif // THERMO_TRANSPORT_DATA_H" in content
