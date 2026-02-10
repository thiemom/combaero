#!/usr/bin/env python3
"""Generate Cantera-compatible NASA-9 format file from CombAero NASA-9 data.

This script converts NASA-9 polynomial data from NASA9_coeffs.json into
Cantera's Chemkin-style NASA-9 format, which can then be converted to YAML
using Cantera's ck2yaml utility.

Usage:
    python generate_cantera_nasa9.py --output combaero_nasa9.txt
    ck2yaml --input=combaero_nasa9.txt --output=combaero_nasa9.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

# Atomic masses for molar mass calculation (IUPAC 2021)
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "AR": 39.948,
}

# Species elemental composition (for CombAero's 14 species)
SPECIES_COMPOSITION = {
    "N2": {"N": 2},
    "O2": {"O": 2},
    "AR": {"AR": 1},
    "CO2": {"C": 1, "O": 2},
    "H2O": {"H": 2, "O": 1},
    "CH4": {"C": 1, "H": 4},
    "C2H6": {"C": 2, "H": 6},
    "C3H8": {"C": 3, "H": 8},
    "IC4H10": {"C": 4, "H": 10},
    "NC5H12": {"C": 5, "H": 12},
    "NC6H14": {"C": 6, "H": 14},
    "NC7H16": {"C": 7, "H": 16},
    "CO": {"C": 1, "O": 1},
    "H2": {"H": 2},
}


def calculate_molar_mass(composition: dict[str, int]) -> float:
    """Calculate molar mass from elemental composition."""
    mass = 0.0
    for element, count in composition.items():
        mass += ATOMIC_MASSES.get(element.upper(), 0.0) * count
    return mass


def format_scientific_d(value: float) -> str:
    """Format number in scientific notation with D exponent (Fortran style)."""
    if value == 0.0:
        return " 0.000000000D+00"

    # Format with 9 decimal places in scientific notation
    s = f"{value:16.9E}"
    # Replace 'E' with 'D' for Fortran compatibility
    s = s.replace("E", "D")
    return s


def write_nasa9_chemkin(
    output_path: Path,
    nasa9_data: dict[str, Any],
    species_list: list[str] | None = None,
) -> None:
    """Write NASA-9 data in Cantera's Chemkin format.

    Args:
        output_path: Path to output file
        nasa9_data: NASA-9 coefficient data from JSON
        species_list: List of species to include (None = all)
    """
    if species_list is None:
        species_list = list(nasa9_data["species"].keys())

    # Normalize species names (case-insensitive matching)
    species_map = {name.upper(): name for name in nasa9_data["species"].keys()}

    with open(output_path, "w") as f:
        # Header
        f.write("! NASA-9 polynomial data for CombAero species\n")
        f.write("! Generated from NASA CEA database\n")
        f.write("! Source: NASA9_coeffs.json\n")
        f.write("!\n")

        # Elements section
        f.write("elements\n")
        elements = set()
        for species in species_list:
            species_upper = species.upper()
            if species_upper in SPECIES_COMPOSITION:
                elements.update(SPECIES_COMPOSITION[species_upper].keys())
        f.write(" ".join(sorted(elements)) + "\n")
        f.write("end\n\n")

        # Species section
        f.write("species\n")
        f.write(" ".join(species_list) + "\n")
        f.write("end\n\n")

        # Thermo section
        f.write("thermo nasa9\n")

        # Temperature ranges (find min and max across all species)
        all_temps = set()
        for species in species_list:
            species_upper = species.upper()
            if species_upper not in species_map:
                continue

            species_key = species_map[species_upper]
            ranges = nasa9_data["species"][species_key]

            for range_data in ranges:
                all_temps.add(range_data["T_min"])
                all_temps.add(range_data["T_max"])

        temp_list = sorted(all_temps)
        # Write temperature range line (up to 4 values)
        f.write(
            f"{temp_list[0]:8.3f} {temp_list[1]:8.3f} {temp_list[2]:8.3f} {temp_list[-1]:8.3f}\n"
        )

        # Species data
        for species in species_list:
            species_upper = species.upper()

            if species_upper not in species_map:
                print(f"Warning: Species {species} not found in NASA-9 data")
                continue

            if species_upper not in SPECIES_COMPOSITION:
                print(f"Warning: Composition not defined for {species}")
                continue

            species_key = species_map[species_upper]
            ranges = nasa9_data["species"][species_key]
            composition = SPECIES_COMPOSITION[species_upper]
            molar_mass = calculate_molar_mass(composition)

            # Species header line (80 characters total)
            # Format: NAME(18) DATE(6) FORMULA(20) PHASE(1) TLOW(10) THIGH(10) COMMON(5) MOLAR_MASS(10)
            species_name = species.ljust(18)
            date = "NASA  "  # 6 chars

            # Formula section (20 chars): element counts
            formula_parts = []
            for elem in ["C", "H", "O", "N", "AR"]:
                count = composition.get(elem, 0)
                if count > 0:
                    formula_parts.append(f"{elem:2s}{count:5.2f}")
            formula = "".join(formula_parts[:4]).ljust(20)  # Max 4 elements

            phase = "G"  # Gas phase

            # Temperature range for this species
            T_low = ranges[0]["T_min"]
            T_high = ranges[-1]["T_max"]

            # Common temperature (usually mid-range)
            T_common = ranges[0]["T_max"] if len(ranges) > 1 else (T_low + T_high) / 2

            # Number of temperature ranges
            # Number of temperature ranges: len(ranges)

            # Reference enthalpy (from first range, coefficient a8)
            # Convert to formation enthalpy at 298.15 K (approximation)
            H_ref = ranges[0]["coeffs"][7] if len(ranges[0]["coeffs"]) > 7 else 0.0

            # Write species header
            f.write(
                f"{species_name}{date}{formula}{phase}{T_low:10.3f}{T_high:10.3f}{T_common:5.0f}.{molar_mass:13.7f}{H_ref:15.3f}\n"
            )

            # Write temperature ranges
            for range_data in ranges:
                T_min = range_data["T_min"]
                T_max = range_data["T_max"]
                coeffs = range_data["coeffs"]

                # Ensure we have 10 coefficients (a1-a7, a8, a9, and reference H)
                # NASA-9 format: a1-a7 are polynomial coeffs, a8 is H integration, a9 is S integration
                if len(coeffs) < 10:
                    # Pad with zeros if needed
                    coeffs = coeffs + [0.0] * (10 - len(coeffs))

                # Reference enthalpy for this range (used in exponent line)
                H_ref_range = coeffs[7] if len(coeffs) > 7 else 0.0

                # Temperature range line with exponents
                # Format from Cantera docs example:
                # "    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0  8671.104"
                # The number after T_max is the number of intervals (always 7 for NASA-9)
                # followed by exponents and reference enthalpy
                f.write(
                    f"    {T_min:8.3f}   {T_max:8.3f}7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0{H_ref_range:12.3f}\n"
                )

                # Coefficient lines (3 coefficients per line, 5 values per coefficient field)
                # Line 1: a1, a2, a3
                f.write(
                    f" {format_scientific_d(coeffs[0])}{format_scientific_d(coeffs[1])}{format_scientific_d(coeffs[2])}\n"
                )

                # Line 2: a4, a5, a6
                f.write(
                    f" {format_scientific_d(coeffs[3])}{format_scientific_d(coeffs[4])}{format_scientific_d(coeffs[5])}\n"
                )

                # Line 3: a7, a8 (H integration), a9 (S integration)
                # Note: Last line has different format - a7, then two blank spaces, then a8 and a9
                f.write(
                    f" {format_scientific_d(coeffs[6])}                 {format_scientific_d(coeffs[7])}{format_scientific_d(coeffs[8])}\n"
                )

        f.write("END\n")

    print(f"Generated Cantera NASA-9 file: {output_path}")
    print(f"Species included: {', '.join(species_list)}")
    print("\nNext step: Convert to YAML using:")
    print(f"  ck2yaml --input={output_path} --output={output_path.stem}.yaml")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Cantera NASA-9 format file from CombAero data"
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("NASA9_coeffs.json"),
        help="Input NASA-9 JSON file (default: NASA9_coeffs.json)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("combaero_nasa9.txt"),
        help="Output Chemkin-style NASA-9 file (default: combaero_nasa9.txt)",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species to include (default: all 14 CombAero species)",
    )

    args = parser.parse_args()

    # Load NASA-9 data
    with open(args.input) as f:
        nasa9_data = json.load(f)

    # Default species list (all 14 CombAero species)
    if args.species is None:
        args.species = [
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
        ]

    # Generate Chemkin file
    write_nasa9_chemkin(args.output, nasa9_data, args.species)


if __name__ == "__main__":
    main()
