#!/usr/bin/env python3
"""Generate Cantera YAML file directly with NASA-9 polynomials from CombAero data.

This script converts NASA-9 polynomial data from NASA9_coeffs.json directly into
Cantera's YAML format, bypassing the Chemkin intermediate format.

Usage:
    python generate_cantera_nasa9_yaml.py --output combaero_nasa9.yaml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import yaml

# Atomic masses for molar mass calculation (IUPAC 2021)
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "Ar": 39.948,
    "AR": 39.948,  # Alias for compatibility
}

# Species elemental composition (for CombAero's species)
# Note: Use Cantera-compatible element names (Ar not AR)
SPECIES_COMPOSITION = {
    "N2": {"N": 2},
    "O2": {"O": 2},
    "AR": {"Ar": 1},  # Cantera uses 'Ar' not 'AR'
    "CO2": {"C": 1, "O": 2},
    "H2O": {"H": 2, "O": 1},
    "CH4": {"C": 1, "H": 4},
    "C2H6": {"C": 2, "H": 6},
    "C3H8": {"C": 3, "H": 8},
    "CO": {"C": 1, "O": 1},
    "H2": {"H": 2},
}


def calculate_molar_mass(composition: dict[str, int]) -> float:
    """Calculate molar mass from elemental composition."""
    mass = 0.0
    for element, count in composition.items():
        mass += ATOMIC_MASSES.get(element.upper(), 0.0) * count
    return mass


def create_cantera_yaml(
    output_path: Path,
    nasa9_data: dict[str, Any],
    species_list: list[str] | None = None,
) -> None:
    """Create Cantera YAML file with NASA-9 polynomials.

    Args:
        output_path: Path to output YAML file
        nasa9_data: NASA-9 coefficient data from JSON
        species_list: List of species to include (None = all available)
    """
    if species_list is None:
        species_list = list(nasa9_data["species"].keys())

    # Normalize species names (case-insensitive matching)
    species_map = {name.upper(): name for name in nasa9_data["species"].keys()}

    # Collect elements
    elements = set()
    for species in species_list:
        species_upper = species.upper()
        if species_upper in SPECIES_COMPOSITION:
            elements.update(SPECIES_COMPOSITION[species_upper].keys())

    # Create phase definition
    phase_data = {
        "name": "gas",
        "thermo": "ideal-gas",
        "elements": sorted(elements),
        "species": species_list,
        "state": {"T": 300.0, "P": 101325.0},
    }

    # Create species definitions
    species_data = []

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

        # Create thermo data for each temperature range
        thermo_ranges = []
        for range_data in ranges:
            T_min = range_data["T_min"]
            T_max = range_data["T_max"]
            coeffs = range_data["coeffs"]

            # Ensure we have 9 coefficients (a1-a7, a8, a9)
            if len(coeffs) < 9:
                coeffs = coeffs + [0.0] * (9 - len(coeffs))

            # Cantera NASA-9 format uses [a1, a2, a3, a4, a5, a6, a7, a8, a9]
            # where a8 and a9 are integration constants
            thermo_range = {
                "T-min": T_min,
                "T-max": T_max,
                "data": coeffs[:9],  # a1-a7, a8, a9
            }
            thermo_ranges.append(thermo_range)

        # Create species entry
        species_entry = {
            "name": species,
            "composition": composition,
            "thermo": {
                "model": "NASA9",
                "temperature-ranges": [r["T-min"] for r in thermo_ranges]
                + [thermo_ranges[-1]["T-max"]],
                "data": [r["data"] for r in thermo_ranges],
            },
        }

        species_data.append(species_entry)

    # Create complete YAML structure
    yaml_data = {
        "description": "CombAero NASA-9 polynomial data",
        "generator": "generate_cantera_nasa9_yaml.py",
        "cantera-version": "3.0.0",
        "date": "2026-02-10",
        "phases": [phase_data],
        "species": species_data,
    }

    # Write YAML file
    with open(output_path, "w") as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Generated Cantera YAML file: {output_path}")
    print(f"Species included: {', '.join(species_list)}")
    print(f"Elements: {', '.join(sorted(elements))}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Cantera YAML file with NASA-9 polynomials"
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
        default=Path("combaero_nasa9.yaml"),
        help="Output Cantera YAML file (default: combaero_nasa9.yaml)",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Species to include (default: all available with composition data)",
    )

    args = parser.parse_args()

    # Load NASA-9 data
    with open(args.input) as f:
        nasa9_data = json.load(f)

    # Default species list (available species with composition data)
    if args.species is None:
        args.species = ["N2", "O2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8", "CO", "H2"]

    # Generate YAML file
    create_cantera_yaml(args.output, nasa9_data, args.species)


if __name__ == "__main__":
    main()
