#!/usr/bin/env python3
"""Generate thermo_transport_data.h from thermo and transport data sources.

This script supports multiple input formats:
1. Combined Cantera YAML mechanism (thermo + transport in one file)
2. Merged JSON from extract_species_data.py (NASA-9 thermo + transport)

Usage:
    # From combined YAML mechanism (NASA-7)
    python generate_thermo_data.py --mechanism gri30.yaml --species N2,O2,CH4

    # From merged JSON (NASA-9 preferred)
    python generate_thermo_data.py --json merged_species.json --species N2,O2,CH4 --prefer-nasa9

    # List available species
    python generate_thermo_data.py --mechanism gri30.yaml --list-species
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import TextIO


class NASAFormat(Enum):
    NASA7 = auto()
    NASA9 = auto()


class Geometry(Enum):
    ATOM = "atom"
    LINEAR = "linear"
    NONLINEAR = "nonlinear"


@dataclass
class NASACoeffs:
    """NASA polynomial coefficients for a species."""

    format: NASAFormat
    T_ranges: list[float]  # [T_low, T_mid, T_high] for NASA7, can have more for NASA9
    coeffs: list[list[float]]  # One set of coefficients per temperature range


@dataclass
class TransportProps:
    """Lennard-Jones transport properties."""

    geometry: Geometry
    well_depth: float  # K
    diameter: float  # Angstrom
    polarizability: float = 0.0  # Angstrom^3
    rotational_relaxation: float = 0.0


@dataclass
class MolecularStructure:
    """Atomic composition."""

    C_atoms: int = 0
    H_atoms: int = 0
    O_atoms: int = 0
    N_atoms: int = 0
    Ar_atoms: int = 0


@dataclass
class SpeciesData:
    """Complete data for a species."""

    name: str
    molar_mass: float
    nasa: NASACoeffs
    transport: TransportProps | None = None
    structure: MolecularStructure = field(default_factory=MolecularStructure)


# Atomic masses for molar mass calculation
ATOMIC_MASSES = {
    "C": 12.011,
    "H": 1.008,
    "O": 15.999,
    "N": 14.007,
    "AR": 39.948,
    "HE": 4.0026,
    "S": 32.065,
    "CL": 35.453,
    "F": 18.998,
}


def calculate_molar_mass(composition: dict[str, int]) -> float:
    """Calculate molar mass from atomic composition."""
    return sum(count * ATOMIC_MASSES.get(elem.upper(), 0.0) for elem, count in composition.items())


# -----------------------------------------------------------------------------
# Cantera YAML Parser
# -----------------------------------------------------------------------------


def _fix_yaml_species_name(name) -> str:
    """Fix species names that YAML misparses as booleans.

    YAML 1.1 parses 'NO' as boolean False and 'ON' as boolean True.
    """
    if name is False:
        return "NO"
    if name is True:
        return "ON"
    return str(name).upper()


def load_cantera_yaml(
    path: Path, species_filter: set[str] | None = None
) -> tuple[list[SpeciesData], NASAFormat]:
    """Load species data from a Cantera YAML mechanism file."""
    import yaml

    with open(path) as f:
        data = yaml.safe_load(f)

    species_list = data.get("species", [])
    if not species_list:
        raise ValueError(f"No species found in {path}")

    # Normalize filter to uppercase
    if species_filter:
        species_filter = {s.upper() for s in species_filter}

    result: list[SpeciesData] = []
    detected_format: NASAFormat | None = None

    for sp in species_list:
        # YAML may parse some names (like NO) as booleans
        name = _fix_yaml_species_name(sp["name"])

        if species_filter and name not in species_filter:
            continue

        # Composition
        comp = sp.get("composition", {})
        structure = MolecularStructure(
            C_atoms=comp.get("C", 0),
            H_atoms=comp.get("H", 0),
            O_atoms=comp.get("O", 0),
            N_atoms=comp.get("N", 0),
            Ar_atoms=comp.get("Ar", 0),
        )
        molar_mass = calculate_molar_mass(comp)

        # Thermo
        thermo = sp.get("thermo", {})
        model = thermo.get("model", "NASA7").upper()

        if model == "NASA7":
            fmt = NASAFormat.NASA7
        elif model == "NASA9":
            fmt = NASAFormat.NASA9
        else:
            print(f"Warning: Unknown thermo model '{model}' for {name}, assuming NASA7")
            fmt = NASAFormat.NASA7

        if detected_format is None:
            detected_format = fmt
        elif detected_format != fmt:
            raise ValueError(f"Mixed NASA formats detected: {detected_format} and {fmt}")

        T_ranges = thermo.get("temperature-ranges", [])
        coeffs_data = thermo.get("data", [])

        nasa = NASACoeffs(
            format=fmt,
            T_ranges=T_ranges,
            coeffs=coeffs_data,
        )

        # Transport
        transport_data = sp.get("transport", {})
        transport: TransportProps | None = None
        if transport_data:
            geom_str = transport_data.get("geometry", "nonlinear")
            try:
                geom = Geometry(geom_str)
            except ValueError:
                geom = Geometry.NONLINEAR

            transport = TransportProps(
                geometry=geom,
                well_depth=transport_data.get("well-depth", 0.0),
                diameter=transport_data.get("diameter", 0.0),
                polarizability=transport_data.get("polarizability", 0.0),
                rotational_relaxation=transport_data.get("rotational-relaxation", 0.0),
            )

        result.append(
            SpeciesData(
                name=name,
                molar_mass=molar_mass,
                nasa=nasa,
                transport=transport,
                structure=structure,
            )
        )

    return result, detected_format or NASAFormat.NASA7


def list_species_in_yaml(path: Path) -> list[str]:
    """List all species names in a Cantera YAML file."""
    import yaml

    with open(path) as f:
        data = yaml.safe_load(f)

    return [_fix_yaml_species_name(sp["name"]) for sp in data.get("species", [])]


# -----------------------------------------------------------------------------
# NASA Thermo Database Parser (for future NASA9 support)
# -----------------------------------------------------------------------------


def detect_nasa_format(path: Path) -> NASAFormat:
    """Detect whether a thermo database is NASA7 or NASA9 format."""
    with open(path) as f:
        content = f.read(2000)  # Read first 2KB

    # NASA9 typically has 9 coefficients and different header format
    # NASA7 (Chemkin) has fixed-width format with specific column positions

    if "thermo nasa9" in content.lower() or re.search(r"\d+\.\d+\s+\d+\.\d+\s+9\s+", content):
        return NASAFormat.NASA9

    return NASAFormat.NASA7


def load_nasa7_thermo(path: Path, species_filter: set[str] | None = None) -> list[SpeciesData]:
    """Load NASA7 format thermo database (Chemkin format)."""
    # TODO: Implement Chemkin NASA7 parser
    raise NotImplementedError("Standalone NASA7 parser not yet implemented. Use YAML mechanism.")


def load_nasa9_thermo(path: Path, species_filter: set[str] | None = None) -> list[SpeciesData]:
    """Load NASA9 format thermo database (NASA Glenn format)."""
    # TODO: Implement NASA9 parser
    raise NotImplementedError("NASA9 parser not yet implemented.")


# -----------------------------------------------------------------------------
# JSON Loader (from extract_species_data.py output)
# -----------------------------------------------------------------------------


def load_merged_json(
    path: Path,
    species_filter: set[str] | None = None,
    prefer_nasa9: bool = True,
) -> tuple[list[SpeciesData], NASAFormat]:
    """Load species data from merged JSON (output of extract_species_data.py).

    Args:
        path: Path to JSON file
        species_filter: Optional set of species names to include
        prefer_nasa9: If True, use NASA-9 thermo when available; else NASA-7

    Returns:
        Tuple of (species_data list, detected NASA format)
    """
    with open(path) as f:
        data = json.load(f)

    result: list[SpeciesData] = []
    has_nasa9 = False

    for sp in data.get("species", []):
        name = sp.get("name", "").upper()
        name_norm = sp.get("name_normalized", name)

        if species_filter and name_norm not in {s.upper() for s in species_filter}:
            continue

        # Determine which thermo format to use
        has_n9 = "thermo_nasa9" in sp
        has_n7 = "thermo_nasa7" in sp

        if prefer_nasa9 and has_n9:
            nasa9 = sp["thermo_nasa9"]
            intervals = nasa9.get("intervals", [])
            # Build T_ranges from intervals: [T0, T1, T2, ...]
            t_ranges = []
            for iv in intervals:
                if not t_ranges:
                    t_ranges.append(iv["T_min"])
                t_ranges.append(iv["T_max"])
            coeffs = [iv["coeffs"] for iv in intervals]
            nasa = NASACoeffs(format=NASAFormat.NASA9, T_ranges=t_ranges, coeffs=coeffs)
            has_nasa9 = True
        elif has_n7:
            nasa7 = sp["thermo_nasa7"]
            t_ranges = [nasa7.get("T_min"), nasa7.get("T_mid"), nasa7.get("T_max")]
            coeffs = [nasa7.get("coeffs_low", []), nasa7.get("coeffs_high", [])]
            nasa = NASACoeffs(format=NASAFormat.NASA7, T_ranges=t_ranges, coeffs=coeffs)
        else:
            print(f"Warning: {name} has no thermo data, skipping")
            continue

        # Transport
        transport: TransportProps | None = None
        tr = sp.get("transport", {})
        if tr.get("geometry"):
            try:
                geom = Geometry(tr["geometry"])
            except ValueError:
                geom = Geometry.NONLINEAR
            transport = TransportProps(
                geometry=geom,
                well_depth=tr.get("well_depth", 0.0) or 0.0,
                diameter=tr.get("diameter", 0.0) or 0.0,
                polarizability=tr.get("polarizability", 0.0) or 0.0,
                rotational_relaxation=tr.get("rotational_relaxation", 0.0) or 0.0,
            )

        # Composition
        comp = sp.get("composition", {})
        structure = MolecularStructure(
            C_atoms=comp.get("C", 0),
            H_atoms=comp.get("H", 0),
            O_atoms=comp.get("O", 0),
            N_atoms=comp.get("N", 0),
            Ar_atoms=comp.get("AR", 0),
        )

        molar_mass = sp.get("molar_mass", 0.0)

        result.append(
            SpeciesData(
                name=name,
                molar_mass=molar_mass,
                nasa=nasa,
                transport=transport,
                structure=structure,
            )
        )

    detected_format = NASAFormat.NASA9 if has_nasa9 else NASAFormat.NASA7
    return result, detected_format


def list_species_in_json(path: Path) -> list[str]:
    """List all species names in a merged JSON file."""
    with open(path) as f:
        data = json.load(f)
    return [sp.get("name", "") for sp in data.get("species", [])]


# -----------------------------------------------------------------------------
# C++ Header Generator
# -----------------------------------------------------------------------------


def species_sort_key(sp: SpeciesData) -> tuple:
    """Sort species into human-friendly groups."""
    name = sp.name.upper()

    # Air species first, in specific order
    air_order = ["N2", "O2", "AR", "CO2", "H2O"]
    if name in air_order:
        return (0, air_order.index(name), name)

    C_atoms, H_atoms, O_atoms, N_atoms = (
        sp.structure.C_atoms,
        sp.structure.H_atoms,
        sp.structure.O_atoms,
        sp.structure.N_atoms,
    )

    # Inert species (no C, H, O)
    if C_atoms == 0 and H_atoms == 0 and O_atoms == 0:
        return (1, name)

    # Hydrocarbons (C>0, H>0, no O or N)
    if C_atoms > 0 and H_atoms > 0 and O_atoms == 0 and N_atoms == 0:
        return (2, C_atoms, name)

    # Other carbon-containing species
    if C_atoms > 0:
        return (3, name)

    # Everything else
    return (4, name)


def format_double(value: float) -> str:
    """Format a double, handling NaN."""
    if math.isnan(value):
        return "std::numeric_limits<double>::quiet_NaN()"
    return str(value)


def _format_nasa9_intervals(nasa: NASACoeffs) -> str:
    """Format NASA-9 intervals for C++ initializer list.

    NASA-9 can have 1-3 temperature intervals. Each interval has:
    - T_min, T_max: temperature bounds
    - 10 coefficients: a1-a7 (Cp), unused, b1 (H), b2 (S)
    """
    interval_strs = []
    for i, coeffs in enumerate(nasa.coeffs):
        # Determine temperature bounds from T_ranges
        # T_ranges for NASA-9: [T0, T1, T2, ...] where intervals are [T0,T1], [T1,T2], etc.
        if len(nasa.T_ranges) > i + 1:
            t_min = nasa.T_ranges[i]
            t_max = nasa.T_ranges[i + 1]
        else:
            # Fallback: use standard CEA breakpoints
            t_breaks = [200.0, 1000.0, 6000.0, 20000.0]
            t_min = t_breaks[i] if i < len(t_breaks) else 200.0
            t_max = t_breaks[i + 1] if i + 1 < len(t_breaks) else 6000.0

        # Pad coefficients to 10 if needed
        padded = list(coeffs) + [0.0] * (10 - len(coeffs))
        coeffs_str = "{" + ", ".join(str(c) for c in padded[:10]) + "}"

        interval_strs.append(f"{{{t_min}, {t_max}, {coeffs_str}}}")

    return "{" + ", ".join(interval_strs) + "}"


def generate_cpp_header(
    species: list[SpeciesData],
    nasa_format: NASAFormat,
    output: TextIO,
) -> None:
    """Generate the C++ header file."""

    # Sort species
    species = sorted(species, key=species_sort_key)

    output.write("""#ifndef THERMO_TRANSPORT_DATA_H
#define THERMO_TRANSPORT_DATA_H

#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <limits>

// NASA polynomial format used in this build
// NASA7: Cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
// NASA9: Cp/R = a1/T^2 + a2/T + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
""")

    output.write(f"// This file uses {nasa_format.name} format\n")
    output.write(
        f"constexpr bool USE_NASA9 = {'true' if nasa_format == NASAFormat.NASA9 else 'false'};\n\n"
    )

    # NASA coefficients struct
    if nasa_format == NASAFormat.NASA7:
        output.write("""struct NASA7_Coeffs {
    double T_low;
    double T_mid;
    double T_high;
    std::vector<double> low_coeffs;   // 7 coefficients for T_low <= T < T_mid
    std::vector<double> high_coeffs;  // 7 coefficients for T_mid <= T <= T_high
};

using NASA_Coeffs = NASA7_Coeffs;

""")
    else:
        output.write("""struct NASA9_Interval {
    double T_min;
    double T_max;
    std::array<double, 10> coeffs;  // a1-a7, unused, b1 (H), b2 (S)
};

struct NASA9_Coeffs {
    std::vector<NASA9_Interval> intervals;  // 1-3 temperature intervals
};

using NASA_Coeffs = NASA9_Coeffs;

""")

    output.write("""struct Transport_Props {
    std::string geometry;
    double well_depth;
    double diameter;
    double polarizability;
};

struct Molecular_Structure {
    std::size_t C;
    std::size_t H;
    std::size_t O;
    std::size_t N;
};

""")

    # Species names
    output.write("const std::vector<std::string> species_names = {")
    output.write(", ".join(f'"{sp.name}"' for sp in species))
    output.write("};\n\n")

    # Species index map
    output.write("const std::unordered_map<std::string, int> species_index = {\n")
    output.write(",\n".join(f'    {{"{sp.name}", {i}}}' for i, sp in enumerate(species)))
    output.write("\n};\n\n")

    # Molar masses
    output.write("const std::vector<double> molar_masses = {")
    output.write(", ".join(str(sp.molar_mass) for sp in species))
    output.write("};\n\n")

    # NASA coefficients
    output.write("const std::vector<NASA_Coeffs> nasa_coeffs = {\n")
    nasa_entries = []

    if nasa_format == NASAFormat.NASA7:
        for sp in species:
            T_low = sp.nasa.T_ranges[0] if len(sp.nasa.T_ranges) > 0 else 300.0
            T_mid = sp.nasa.T_ranges[1] if len(sp.nasa.T_ranges) > 1 else 1000.0
            T_high = sp.nasa.T_ranges[2] if len(sp.nasa.T_ranges) > 2 else 5000.0

            low_coeffs = sp.nasa.coeffs[0] if len(sp.nasa.coeffs) > 0 else []
            high_coeffs = sp.nasa.coeffs[1] if len(sp.nasa.coeffs) > 1 else []

            low_str = "{" + ", ".join(str(c) for c in low_coeffs) + "}"
            high_str = "{" + ", ".join(str(c) for c in high_coeffs) + "}"

            nasa_entries.append(f"{{{T_low}, {T_mid}, {T_high}, {low_str}, {high_str}}}")
    else:
        # NASA-9: variable number of intervals (1-3)
        for sp in species:
            intervals_str = _format_nasa9_intervals(sp.nasa)
            nasa_entries.append(f"{{{intervals_str}}}")

    output.write(",\n".join(nasa_entries))
    output.write("\n};\n\n")

    # Transport properties
    output.write("const std::vector<Transport_Props> transport_props = {\n")
    transport_entries = []
    for sp in species:
        if sp.transport:
            t = sp.transport
            pol_str = format_double(t.polarizability)
            transport_entries.append(
                f'{{"{t.geometry.value}", {t.well_depth}, {t.diameter}, {pol_str}}}'
            )
        else:
            transport_entries.append('{"nonlinear", 0.0, 0.0, 0.0}')

    output.write(",\n".join(transport_entries))
    output.write("\n};\n\n")

    # Molecular structures
    output.write("const std::vector<Molecular_Structure> molecular_structures = {\n")
    struct_entries = []
    for sp in species:
        s = sp.structure
        struct_entries.append(f"{{{s.C}, {s.H}, {s.O}, {s.N}}}")

    output.write(",\n".join(struct_entries))
    output.write("\n};\n\n")

    output.write("#endif // THERMO_TRANSPORT_DATA_H\n")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def parse_species_list(s: str) -> list[str]:
    """Parse comma-separated species list."""
    return [sp.strip() for sp in s.split(",") if sp.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate thermo_transport_data.h from thermo/transport sources.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input sources
    parser.add_argument(
        "--mechanism",
        "-m",
        type=Path,
        help="Cantera YAML mechanism file (combined thermo + transport)",
    )
    parser.add_argument(
        "--json",
        "-j",
        type=Path,
        help="Merged JSON from extract_species_data.py",
    )
    parser.add_argument(
        "--prefer-nasa9",
        action="store_true",
        default=True,
        help="Prefer NASA-9 thermo when available (default: True)",
    )
    parser.add_argument(
        "--prefer-nasa7",
        action="store_true",
        help="Prefer NASA-7 thermo when available",
    )

    # Species selection
    parser.add_argument(
        "--species",
        "-s",
        type=str,
        help="Comma-separated list of species to include",
    )
    parser.add_argument(
        "--list-species",
        action="store_true",
        help="List available species and exit",
    )

    # Output
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("thermo_transport_data.h"),
        help="Output header file (default: thermo_transport_data.h)",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.mechanism and args.json:
        parser.error("Cannot use --mechanism with --json")

    if not args.mechanism and not args.json:
        parser.error("Must specify either --mechanism or --json")

    # List species mode
    if args.list_species:
        if args.mechanism:
            species_names = list_species_in_yaml(args.mechanism)
            print(f"Species in {args.mechanism} ({len(species_names)} total):")
        elif args.json:
            species_names = list_species_in_json(args.json)
            print(f"Species in {args.json} ({len(species_names)} total):")
        else:
            parser.error("--list-species requires --mechanism or --json")
            return 1
        for name in sorted(species_names):
            print(f"  {name}")
        return 0

    # Parse species filter
    species_filter: set[str] | None = None
    if args.species:
        species_filter = set(parse_species_list(args.species))
        print(f"Filtering to {len(species_filter)} species")

    # Determine NASA format preference
    prefer_nasa9 = not args.prefer_nasa7

    # Load data
    if args.mechanism:
        print(f"Loading mechanism: {args.mechanism}")
        species_data, nasa_format = load_cantera_yaml(args.mechanism, species_filter)
    elif args.json:
        print(f"Loading JSON: {args.json}")
        species_data, nasa_format = load_merged_json(args.json, species_filter, prefer_nasa9)
    else:
        parser.error("Must specify --mechanism or --json")
        return 1

    print(f"Loaded {len(species_data)} species, format: {nasa_format.name}")

    # Check for missing transport
    missing_transport = [sp.name for sp in species_data if sp.transport is None]
    if missing_transport:
        print(
            f"Warning: {len(missing_transport)} species missing transport data: "
            f"{missing_transport[:5]}..."
        )

    # Generate header
    print(f"Generating: {args.output}")
    with open(args.output, "w") as f:
        generate_cpp_header(species_data, nasa_format, f)

    print("Done!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
