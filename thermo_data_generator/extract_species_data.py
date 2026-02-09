#!/usr/bin/env python3
"""Unified species data extraction tool.

Extracts thermodynamic and transport data from:
- Cantera YAML mechanism files (NASA-7 polynomials + transport)
- NASA CEA web output (NASA-9 polynomials, higher temperature range)

Can merge data from both sources, preferring NASA-9 for thermo (wider T range)
and YAML for transport properties.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import yaml


# =============================================================================
# Atomic masses (IUPAC 2021 standard atomic weights)
# =============================================================================

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


# =============================================================================
# Species name normalization and matching
# =============================================================================


def normalize_species_name(name: str) -> str:
    """Normalize species name for matching.

    Converts to uppercase and normalizes common naming conventions:
    - C4H10,n-butane -> NC4H10
    - C5H12,i-pentane -> IC5H12
    """
    upper = name.upper()

    # Handle CEA-style names like "C4H10,n-butane" -> "NC4H10"
    if "," in upper:
        formula, suffix = upper.split(",", 1)
        suffix = suffix.strip()
        if suffix.startswith("N-") or suffix.startswith("N "):
            return "N" + formula
        if suffix.startswith("I-") or suffix.startswith("I ") or suffix.startswith("ISO"):
            return "I" + formula
        # Keep original if we can't parse
        return upper

    return upper


def find_similar_species(target: str, available: list[str]) -> list[str]:
    """Find species in available list that might match target."""
    target_norm = normalize_species_name(target)
    similar = []

    # Extract base formula (remove prefixes)
    prefix_pattern = re.compile(
        r"^(N-?|I-?|ISO-?|NEO-?|SEC-?|TERT-?|T-?|S-?|P-?|C-?)?(.+)$", re.IGNORECASE
    )

    target_match = prefix_pattern.match(target_norm)
    if not target_match:
        return similar
    target_base = target_match.group(2)

    for avail in available:
        avail_norm = normalize_species_name(avail)
        avail_match = prefix_pattern.match(avail_norm)
        if not avail_match:
            continue
        avail_base = avail_match.group(2)

        if target_base == avail_base and target_norm != avail_norm:
            similar.append(avail)

    return similar


# =============================================================================
# Cantera YAML extraction (NASA-7 + transport)
# =============================================================================


def extract_from_yaml(
    yaml_file: Path,
    selected_species: list[str] | None = None,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """Extract NASA-7 thermo and transport data from Cantera YAML.

    Returns:
        Tuple of (species_data dict keyed by normalized name, warnings list)
    """
    with open(yaml_file) as f:
        data = yaml.safe_load(f)

    warnings: list[str] = []
    species_data: dict[str, dict[str, Any]] = {}

    available_species = [sp["name"] for sp in data.get("species", [])]

    # Build set of selected species (normalized)
    if selected_species:
        selected_norm = {normalize_species_name(s) for s in selected_species}
    else:
        selected_norm = None  # Extract all

    found_norm: set[str] = set()

    for species in data.get("species", []):
        species_name = species["name"]
        species_norm = normalize_species_name(species_name)

        if selected_norm is not None and species_norm not in selected_norm:
            continue

        found_norm.add(species_norm)

        comp = species.get("composition", {})
        thermo = species.get("thermo", {})
        transport = species.get("transport", {})

        molar_mass = sum(comp.get(el, 0) * get_atomic_mass(el) for el in comp)

        if not transport:
            warnings.append(f"{species_name}: No transport data available")

        temp_ranges = thermo.get("temperature-ranges", [None, None, None])
        thermo_data = thermo.get("data", [[], []])

        entry: dict[str, Any] = {
            "name": species_name,
            "name_normalized": species_norm,
            "molar_mass": molar_mass,
            "composition": {el.upper(): count for el, count in comp.items()},
            "thermo_nasa7": {
                "model": "NASA-7",
                "T_min": temp_ranges[0] if len(temp_ranges) > 0 else None,
                "T_mid": temp_ranges[1] if len(temp_ranges) > 1 else None,
                "T_max": temp_ranges[2] if len(temp_ranges) > 2 else None,
                "coeffs_low": thermo_data[0] if len(thermo_data) > 0 else [],
                "coeffs_high": thermo_data[1] if len(thermo_data) > 1 else [],
            },
            "transport": {
                "geometry": transport.get("geometry", None),
                "well_depth": transport.get("well-depth", None),
                "diameter": transport.get("diameter", None),
                "polarizability": transport.get("polarizability", None),
                "rotational_relaxation": transport.get("rotational-relaxation", None),
                "dipole": transport.get("dipole", None),
            },
        }

        species_data[species_norm] = entry

    # Warn about unmatched species
    if selected_species:
        for sel in selected_species:
            sel_norm = normalize_species_name(sel)
            if sel_norm not in found_norm:
                similar = find_similar_species(sel, available_species)
                if similar:
                    warnings.append(f"YAML: '{sel}' not found. Similar: {similar}")
                else:
                    warnings.append(f"YAML: '{sel}' not found")

    return species_data, warnings


# =============================================================================
# NASA CEA extraction (NASA-9)
# =============================================================================

COEF_MARK = "COEFFICIENTS FOR FITTED THERMODYNAMIC FUNCTIONS"
THERMO_MARK = "THERMODYNAMIC FUNCTIONS CALCULATED"
FLOAT_RE = re.compile(r"[+-]?\d+\.\d+D[+-]\d{2}")
SPECIES_RE = re.compile(r"FOR\s+([A-Za-z0-9\-\+\(\),]+)")
T_BREAKS = [200.0, 1000.0, 6000.0, 20000.0]


def extract_from_cea(
    cea_file: Path,
    selected_species: list[str] | None = None,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """Extract NASA-9 coefficients from CEA web output.

    Returns:
        Tuple of (species_data dict keyed by normalized name, warnings list)
    """
    text = cea_file.read_text()

    warnings: list[str] = []
    species_data: dict[str, dict[str, Any]] = {}

    if selected_species:
        selected_norm = {normalize_species_name(s) for s in selected_species}
    else:
        selected_norm = None

    segments = text.split(COEF_MARK)[1:]
    available_species: list[str] = []

    for seg in segments:
        seg_full = COEF_MARK + seg

        try:
            coeff_block, tail = seg_full.split(THERMO_MARK, 1)
        except ValueError:
            continue

        match = SPECIES_RE.search(tail)
        if not match:
            continue

        name = match.group(1).strip()
        available_species.append(name)
        name_norm = normalize_species_name(name)

        if selected_norm is not None and name_norm not in selected_norm:
            continue

        nums_str = FLOAT_RE.findall(coeff_block)
        coeff_vals = [float(s.replace("D", "E")) for s in nums_str]

        if not coeff_vals or len(coeff_vals) % 10 != 0:
            warnings.append(f"CEA: {name}: Invalid coefficient count ({len(coeff_vals)})")
            continue

        n_intervals = len(coeff_vals) // 10

        if n_intervals == 1:
            bounds = [(T_BREAKS[0], T_BREAKS[2])]
        elif n_intervals == 2:
            bounds = [(T_BREAKS[0], T_BREAKS[1]), (T_BREAKS[1], T_BREAKS[2])]
        elif n_intervals == 3:
            bounds = [
                (T_BREAKS[0], T_BREAKS[1]),
                (T_BREAKS[1], T_BREAKS[2]),
                (T_BREAKS[2], T_BREAKS[3]),
            ]
        else:
            warnings.append(f"CEA: {name}: Unexpected {n_intervals} intervals")
            continue

        intervals = []
        for i, (t_min, t_max) in enumerate(bounds):
            coeffs = coeff_vals[10 * i : 10 * (i + 1)]
            intervals.append(
                {
                    "T_min": t_min,
                    "T_max": t_max,
                    "coeffs": coeffs,
                }
            )

        entry: dict[str, Any] = {
            "name": name,
            "name_normalized": name_norm,
            "thermo_nasa9": {
                "model": "NASA-9",
                "intervals": intervals,
            },
        }

        species_data[name_norm] = entry

    # Warn about unmatched species
    if selected_species:
        found_norm = set(species_data.keys())
        for sel in selected_species:
            sel_norm = normalize_species_name(sel)
            if sel_norm not in found_norm:
                similar = find_similar_species(sel, available_species)
                if similar:
                    warnings.append(f"CEA: '{sel}' not found. Similar: {similar}")
                else:
                    warnings.append(f"CEA: '{sel}' not found")

    return species_data, warnings


# =============================================================================
# Merge data from multiple sources
# =============================================================================


def merge_species_data(
    yaml_data: dict[str, dict[str, Any]],
    cea_data: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Merge species data from YAML and CEA sources.

    Strategy:
    - Use YAML as base (has transport data)
    - Add NASA-9 thermo from CEA where available
    """
    merged: dict[str, dict[str, Any]] = {}

    # Start with YAML data
    for name_norm, entry in yaml_data.items():
        merged[name_norm] = entry.copy()

    # Add/merge CEA data
    for name_norm, cea_entry in cea_data.items():
        if name_norm in merged:
            # Add NASA-9 to existing entry
            merged[name_norm]["thermo_nasa9"] = cea_entry["thermo_nasa9"]
        else:
            # New species from CEA only
            merged[name_norm] = cea_entry

    return merged


# =============================================================================
# CLI
# =============================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract species thermodynamic and transport data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract from YAML only
  %(prog)s --yaml mechanism.yaml -o species.json

  # Extract from CEA only
  %(prog)s --cea cea_output.txt -o species.json

  # Merge both sources (NASA-9 thermo + YAML transport)
  %(prog)s --yaml mechanism.yaml --cea cea_output.txt -o species.json

  # Extract specific species
  %(prog)s --yaml mechanism.yaml -s O2 N2 H2O CO2 -o species.json
""",
    )
    parser.add_argument(
        "--yaml",
        type=Path,
        help="Cantera YAML mechanism file (NASA-7 + transport)",
    )
    parser.add_argument(
        "--cea",
        type=Path,
        help="NASA CEA web output file (NASA-9)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("species_data.json"),
        help="Output JSON file (default: species_data.json)",
    )
    parser.add_argument(
        "-s",
        "--species",
        nargs="+",
        help="Species to extract (default: all available)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available species and exit",
    )

    args = parser.parse_args()

    if not args.yaml and not args.cea:
        parser.error("At least one of --yaml or --cea is required")

    all_warnings: list[str] = []
    yaml_data: dict[str, dict[str, Any]] = {}
    cea_data: dict[str, dict[str, Any]] = {}

    # Extract from YAML
    if args.yaml:
        if not args.yaml.exists():
            print(f"ERROR: File not found: {args.yaml}", file=sys.stderr)
            sys.exit(1)
        yaml_data, warnings = extract_from_yaml(args.yaml, args.species)
        all_warnings.extend(warnings)

    # Extract from CEA
    if args.cea:
        if not args.cea.exists():
            print(f"ERROR: File not found: {args.cea}", file=sys.stderr)
            sys.exit(1)
        cea_data, warnings = extract_from_cea(args.cea, args.species)
        all_warnings.extend(warnings)

    # List mode
    if args.list:
        print("Available species:")
        if yaml_data:
            print(f"\n  YAML ({len(yaml_data)}):")
            for name in sorted(yaml_data.keys()):
                print(f"    {yaml_data[name]['name']}")
        if cea_data:
            print(f"\n  CEA ({len(cea_data)}):")
            for name in sorted(cea_data.keys()):
                print(f"    {cea_data[name]['name']}")
        return

    # Merge data
    if yaml_data and cea_data:
        merged = merge_species_data(yaml_data, cea_data)
        sources = [str(args.yaml), str(args.cea)]
    elif yaml_data:
        merged = yaml_data
        sources = [str(args.yaml)]
    else:
        merged = cea_data
        sources = [str(args.cea)]

    # Print warnings
    for warning in all_warnings:
        print(f"WARNING: {warning}", file=sys.stderr)

    # Build output
    output = {
        "sources": sources,
        "species": list(merged.values()),
    }

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    # Summary
    n_nasa7 = sum(1 for sp in merged.values() if "thermo_nasa7" in sp)
    n_nasa9 = sum(1 for sp in merged.values() if "thermo_nasa9" in sp)
    n_transport = sum(1 for sp in merged.values() if sp.get("transport", {}).get("geometry"))

    print(f"Extracted {len(merged)} species to {args.output}")
    print(f"  NASA-7 thermo: {n_nasa7}, NASA-9 thermo: {n_nasa9}, transport: {n_transport}")
    if all_warnings:
        print(f"  ({len(all_warnings)} warnings, see stderr)")


if __name__ == "__main__":
    main()
