#!/usr/bin/env python3
"""Extract NASA-7 thermo and transport data from Cantera YAML mechanism files.

Outputs JSON with species thermo polynomials and transport properties.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import yaml


# Atomic masses (IUPAC 2021 standard atomic weights)
ATOMIC_MASSES: dict[str, float] = {
    "H": 1.008,
    "HE": 4.0026,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "NE": 20.180,
    "S": 32.06,
    "CL": 35.45,
    "AR": 39.948,
    "BR": 79.904,
    "KR": 83.798,
    "I": 126.90,
    "XE": 131.29,
    "SI": 28.085,
    "P": 30.974,
}


def get_atomic_mass(element: str) -> float:
    """Return atomic mass for element symbol, or 0 with warning if unknown."""
    key = element.upper()
    if key not in ATOMIC_MASSES:
        print(f"WARNING: Unknown element '{element}', using mass=0", file=sys.stderr)
        return 0.0
    return ATOMIC_MASSES[key]


def normalize_species_name(name: str) -> str:
    """Normalize species name for matching (uppercase, strip common prefixes)."""
    return name.upper()


def find_similar_species(target: str, available: list[str]) -> list[str]:
    """Find species in available list that might match target.
    
    Handles common naming conventions:
    - C5H12 matches NC5H12, n-C5H12, nC5H12 (n-alkane)
    - iC4H10 matches IC4H10, i-C4H10, iso-C4H10, isoC4H10
    """
    target_upper = target.upper()
    similar = []
    
    # Extract the base formula (remove prefixes like n-, i-, iso-, neo-, etc.)
    prefix_pattern = re.compile(r'^(N-?|I-?|ISO-?|NEO-?|SEC-?|TERT-?|T-?|S-?|P-?|C-?)?(.+)$', re.IGNORECASE)
    
    target_match = prefix_pattern.match(target_upper)
    if not target_match:
        return similar
    target_prefix = target_match.group(1) or ""
    target_base = target_match.group(2)
    
    for avail in available:
        avail_upper = avail.upper()
        avail_match = prefix_pattern.match(avail_upper)
        if not avail_match:
            continue
        avail_prefix = avail_match.group(1) or ""
        avail_base = avail_match.group(2)
        
        # Same base formula?
        if target_base == avail_base:
            # Different prefix or one has no prefix
            if target_upper != avail_upper:
                similar.append(avail)
    
    return similar


def extract_species_data(
    yaml_file: Path,
    selected_species: list[str],
) -> tuple[list[dict], list[str]]:
    """Extract thermo and transport data for selected species.
    
    Returns:
        Tuple of (species_data list, warnings list)
    """
    with open(yaml_file) as f:
        data = yaml.safe_load(f)
    
    warnings: list[str] = []
    species_data: list[dict] = []
    
    # Build lookup of available species
    available_species = [sp["name"] for sp in data.get("species", [])]
    available_upper = {name.upper(): name for name in available_species}
    
    # Track which selected species were found
    found_species: set[str] = set()
    
    for species in data.get("species", []):
        species_name = species["name"]
        species_upper = species_name.upper()
        
        # Check if this species is in our selected list
        selected_upper = {s.upper() for s in selected_species}
        if species_upper not in selected_upper:
            continue
        
        found_species.add(species_upper)
        
        comp = species.get("composition", {})
        thermo = species.get("thermo", {})
        transport = species.get("transport", {})
        
        # Calculate molar mass
        molar_mass = sum(
            comp.get(el, 0) * get_atomic_mass(el) for el in comp
        )
        
        # Check for missing transport data
        if not transport:
            warnings.append(f"{species_name}: No transport data available")
        
        temp_ranges = thermo.get("temperature-ranges", [None, None, None])
        thermo_data = thermo.get("data", [[], []])
        
        entry = {
            "species_name": species_name.upper(),
            "molar_mass": molar_mass,
            "composition": {el.upper(): count for el, count in comp.items()},
            "thermo": {
                "model": thermo.get("model", ""),
                "T_low": temp_ranges[0] if len(temp_ranges) > 0 else None,
                "T_mid": temp_ranges[1] if len(temp_ranges) > 1 else None,
                "T_high": temp_ranges[2] if len(temp_ranges) > 2 else None,
                "coeffs_low": thermo_data[0] if len(thermo_data) > 0 else [],
                "coeffs_high": thermo_data[1] if len(thermo_data) > 1 else [],
            },
            "transport": {
                "geometry": transport.get("geometry", ""),
                "well_depth": transport.get("well-depth", None),
                "diameter": transport.get("diameter", None),
                "polarizability": transport.get("polarizability", None),
                "rotational_relaxation": transport.get("rotational-relaxation", None),
                "dipole": transport.get("dipole", None),
            },
        }
        
        species_data.append(entry)
    
    # Check for unmatched species and suggest similar ones
    for sel in selected_species:
        sel_upper = sel.upper()
        if sel_upper not in found_species:
            similar = find_similar_species(sel, available_species)
            if similar:
                warnings.append(
                    f"'{sel}' not found. Similar species available: {similar}"
                )
            else:
                warnings.append(f"'{sel}' not found in mechanism file")
    
    return species_data, warnings


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract NASA-7 thermo and transport data from Cantera YAML"
    )
    parser.add_argument(
        "yaml_file",
        type=Path,
        help="Input Cantera YAML mechanism file",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("species_data.json"),
        help="Output JSON file (default: species_data.json)",
    )
    parser.add_argument(
        "-s", "--species",
        nargs="+",
        default=[
            "O2", "N2", "AR", "CO2", "H2O", "CH4", "C2H6", "C3H8",
            "C4H10", "iC4H10", "NC5H12", "NC6H14", "NC7H16", "H2", "CO",
        ],
        help="Species to extract (default: common combustion species)",
    )
    
    args = parser.parse_args()
    
    if not args.yaml_file.exists():
        print(f"ERROR: File not found: {args.yaml_file}", file=sys.stderr)
        sys.exit(1)
    
    species_data, warnings = extract_species_data(args.yaml_file, args.species)
    
    # Print warnings
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    
    # Write output
    output = {
        "source": str(args.yaml_file),
        "species": species_data,
    }
    
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"Extracted {len(species_data)} species to {args.output}")
    if warnings:
        print(f"  ({len(warnings)} warnings, see stderr)")


if __name__ == "__main__":
    main()
