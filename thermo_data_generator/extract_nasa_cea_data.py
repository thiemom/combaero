#!/usr/bin/env python3
"""Extract NASA-9 polynomial coefficients from CEA web output.

Parses the text dump from NASA CEA online tool and outputs JSON with
thermodynamic polynomial coefficients for each species.
"""

import argparse
import json
import re
import sys
from pathlib import Path


# Markers in CEA output
COEF_MARK = "COEFFICIENTS FOR FITTED THERMODYNAMIC FUNCTIONS"
THERMO_MARK = "THERMODYNAMIC FUNCTIONS CALCULATED"

# Regex for FORTRAN D-notation floats (e.g., 1.234567890D+02)
FLOAT_RE = re.compile(r"[+-]?\d+\.\d+D[+-]\d{2}")

# Regex for species name after "FOR <NAME>"
SPECIES_RE = re.compile(r"FOR\s+([A-Za-z0-9\-\+\(\),]+)")

# CEA standard temperature breakpoints
T_BREAKS = [200.0, 1000.0, 6000.0, 20000.0]


def parse_cea_output(text: str) -> tuple[dict[str, list[dict]], list[str]]:
    """Parse CEA web output text and extract NASA-9 coefficients.
    
    Returns:
        Tuple of (species_data dict, warnings list)
    """
    species_data: dict[str, list[dict]] = {}
    warnings: list[str] = []
    
    # Split into segments by coefficient block marker
    segments = text.split(COEF_MARK)[1:]
    
    for seg in segments:
        seg_full = COEF_MARK + seg
        
        # Split coefficient block from following thermo table
        try:
            coeff_block, tail = seg_full.split(THERMO_MARK, 1)
        except ValueError:
            warnings.append("Malformed segment: no thermo table found after coefficients")
            continue
        
        # Extract species name from "FOR <NAME>" line
        match = SPECIES_RE.search(tail)
        if not match:
            warnings.append("Could not extract species name from thermo table header")
            continue
        name = match.group(1).strip()
        
        # Extract all D-notation floats from coefficient block
        nums_str = FLOAT_RE.findall(coeff_block)
        coeff_vals = [float(s.replace("D", "E")) for s in nums_str]
        
        if not coeff_vals:
            warnings.append(f"{name}: No coefficients found in block")
            continue
        
        # NASA-9 uses 10 coefficients per interval
        if len(coeff_vals) % 10 != 0:
            warnings.append(
                f"{name}: Found {len(coeff_vals)} coefficients, "
                f"not a multiple of 10. Skipping."
            )
            continue
        
        n_intervals = len(coeff_vals) // 10
        
        # Assign standard CEA temperature ranges based on number of intervals
        if n_intervals == 1:
            bounds = [(T_BREAKS[0], T_BREAKS[2])]  # 200-6000
        elif n_intervals == 2:
            bounds = [
                (T_BREAKS[0], T_BREAKS[1]),  # 200-1000
                (T_BREAKS[1], T_BREAKS[2]),  # 1000-6000
            ]
        elif n_intervals == 3:
            bounds = [
                (T_BREAKS[0], T_BREAKS[1]),  # 200-1000
                (T_BREAKS[1], T_BREAKS[2]),  # 1000-6000
                (T_BREAKS[2], T_BREAKS[3]),  # 6000-20000
            ]
        else:
            warnings.append(
                f"{name}: Unexpected {n_intervals} intervals. Skipping."
            )
            continue
        
        # Slice coefficients into intervals
        intervals = []
        for i, (t_min, t_max) in enumerate(bounds):
            coeffs = coeff_vals[10 * i : 10 * (i + 1)]
            intervals.append({
                "T_min": t_min,
                "T_max": t_max,
                "coeffs": coeffs,
                # NASA-9 coefficient meanings:
                # a[0..6]: T^-2, T^-1, T^0, T^1, T^2, T^3, T^4
                # a[7]: unused (0 for most species, sometimes ln(T) term)
                # a[8]: H integration constant (b1)
                # a[9]: S integration constant (b2)
            })
        
        species_data[name] = intervals
    
    return species_data, warnings


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract NASA-9 coefficients from CEA web output"
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Input CEA web output text file",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("NASA9_coeffs.json"),
        help="Output JSON file (default: NASA9_coeffs.json)",
    )
    
    args = parser.parse_args()
    
    if not args.input_file.exists():
        print(f"ERROR: File not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    text = args.input_file.read_text()
    species_data, warnings = parse_cea_output(text)
    
    # Print warnings
    for warning in warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    
    # Write output
    output = {
        "source": str(args.input_file),
        "model": "NASA-9",
        "species": species_data,
    }
    
    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"Parsed {len(species_data)} species to {args.output}")
    if warnings:
        print(f"  ({len(warnings)} warnings, see stderr)")


if __name__ == "__main__":
    main()
