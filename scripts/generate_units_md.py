#!/usr/bin/env python3
"""Generate docs/UNITS.md from include/units_data.h.

This script parses the units_data.h file and generates a formatted Markdown
document with all function units organized by module.

Usage:
    python scripts/generate_units_md.py
"""

from __future__ import annotations

import re
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class UnitEntry:
    name: str
    input: str
    output: str
    section: str = ""


def parse_units_data(path: Path) -> list[UnitEntry]:
    """Parse units_data.h and extract all unit entries."""
    content = path.read_text()

    entries: list[UnitEntry] = []
    current_section = ""

    # Match section comments like: // thermo.h - Species Data
    section_pattern = re.compile(r'//\s*[-]+\s*\n\s*//\s*(.+?)\s*\n\s*//\s*[-]+')

    # Match entries like: {"name", "input", "output"},
    entry_pattern = re.compile(r'\{"([^"]+)",\s*"([^"]*)",\s*"([^"]*)"\}')

    # Find all sections and their positions
    sections: list[tuple[int, str]] = []
    for match in section_pattern.finditer(content):
        sections.append((match.end(), match.group(1).strip()))

    # Find all entries
    for match in entry_pattern.finditer(content):
        pos = match.start()
        name, input_units, output_units = match.groups()

        # Determine which section this entry belongs to
        section = ""
        for sec_pos, sec_name in reversed(sections):
            if pos > sec_pos:
                section = sec_name
                break

        entries.append(UnitEntry(name, input_units, output_units, section))

    return entries


def group_by_section(entries: list[UnitEntry]) -> dict[str, list[UnitEntry]]:
    """Group entries by their section."""
    grouped: dict[str, list[UnitEntry]] = defaultdict(list)
    for entry in entries:
        grouped[entry.section].append(entry)
    return grouped


def format_table(entries: list[UnitEntry]) -> str:
    """Format a list of entries as a Markdown table."""
    if not entries:
        return ""

    # Calculate column widths
    name_width = max(len(f"`{e.name}`") for e in entries)
    input_width = max(len(e.input) for e in entries)
    output_width = max(len(e.output) for e in entries)

    # Minimum widths for headers
    name_width = max(name_width, len("Function"))
    input_width = max(input_width, len("Input Units"))
    output_width = max(output_width, len("Output Unit"))

    lines = []
    # Header
    lines.append(f"| {'Function':<{name_width}} | {'Input Units':<{input_width}} | {'Output Unit':<{output_width}} |")
    lines.append(f"|{'-' * (name_width + 2)}|{'-' * (input_width + 2)}|{'-' * (output_width + 2)}|")

    # Rows
    for e in entries:
        name_col = f"`{e.name}`"
        lines.append(f"| {name_col:<{name_width}} | {e.input:<{input_width}} | {e.output:<{output_width}} |")

    return "\n".join(lines)


def generate_markdown(entries: list[UnitEntry]) -> str:
    """Generate the full UNITS.md content."""
    grouped = group_by_section(entries)

    lines = [
        "# combaero Units Reference",
        "",
        "This document defines the SI-based unit system used throughout the combaero library.",
        "All functions use consistent units to avoid conversion errors.",
        "",
        "**Auto-generated from `include/units_data.h` - do not edit manually.**",
        "",
        "---",
        "",
        "## Base SI Units",
        "",
        "| Quantity           | Unit   | Symbol |",
        "|--------------------|--------|--------|",
        "| Temperature        | Kelvin | K      |",
        "| Pressure           | Pascal | Pa     |",
        "| Mass               | kilogram | kg   |",
        "| Length             | meter  | m      |",
        "| Time               | second | s      |",
        "| Amount of substance| mole   | mol    |",
        "",
        "---",
        "",
        "## Derived Units",
        "",
        "| Quantity             | Unit        | Symbol     |",
        "|----------------------|-------------|------------|",
        "| Energy               | Joule       | J = kg*m^2/s^2 |",
        "| Power                | Watt        | W = J/s    |",
        "| Force                | Newton      | N = kg*m/s^2 |",
        "| Dynamic viscosity    | Pascal-second | Pa*s     |",
        "| Kinematic viscosity  | -           | m^2/s      |",
        "| Thermal conductivity | -           | W/(m*K)    |",
        "| Specific heat        | -           | J/(mol*K) or J/(kg*K) |",
        "| Enthalpy             | -           | J/mol or J/kg |",
        "| Entropy              | -           | J/(mol*K) or J/(kg*K) |",
        "| Density              | -           | kg/m^3     |",
        "| Velocity             | -           | m/s        |",
        "| Area                 | -           | m^2        |",
        "| Volume               | -           | m^3        |",
        "| Mass flow rate       | -           | kg/s       |",
        "| Volumetric flow rate | -           | m^3/s      |",
        "",
        "---",
        "",
        "## Physical Constants",
        "",
        "| Constant             | Value              | Unit       |",
        "|----------------------|--------------------|------------|",
        "| Universal gas constant R | 8.31446261815324 | J/(mol*K)  |",
        "| Boltzmann constant   | 1.380649e-23       | J/K        |",
        "| Avogadro's number    | 6.02214076e23      | 1/mol      |",
        "| Standard gravity g0  | 9.80665            | m/s^2      |",
        "",
        "---",
        "",
        "## Function Units by Module",
        "",
    ]

    # Section order (roughly matching the header file organization)
    section_order = [
        "thermo.h - Species Data",
        "thermo.h - Dimensionless NASA Polynomials",
        "thermo.h - Per-Species Properties",
        "thermo.h - Mixture Properties (Molar Basis)",
        "thermo.h - Mixture Properties (Mass/Other Basis)",
        "thermo.h - Inverse Solvers",
        "transport.h - Transport Properties",
        "compressible.h - Nozzle Flow",
        "compressible.h - Fanno Flow",
        "compressible.h - Thrust",
        "incompressible.h - Bernoulli & Orifice",
        "incompressible.h - Pipe Flow",
        "friction.h - Friction Factor Correlations",
        "orifice.h - Discharge Coefficients",
        "humidair.h - Humid Air Properties",
        "combustion.h - Stoichiometry",
        "combustion.h - Equivalence Ratio",
        "combustion.h - Mixture Fraction (Bilger)",
        "combustion.h - Complete Combustion",
        "combustion.h - Stream Solvers",
        "equilibrium.h - Chemical Equilibrium",
        "state.h - State Properties",
        "state.h - Stream Properties",
    ]

    # Add sections in order, tracking current header to avoid repetition
    current_header = ""
    for section in section_order:
        if section in grouped:
            # Parse header file and subsection
            if " - " in section:
                header, subsection = section.split(" - ", 1)
                if header != current_header:
                    if current_header:
                        lines.append("---")
                        lines.append("")
                    lines.append(f"### {header}")
                    lines.append("")
                    current_header = header
                lines.append(f"#### {subsection}")
            else:
                if current_header:
                    lines.append("---")
                    lines.append("")
                lines.append(f"### {section}")
                current_header = section
            lines.append("")
            lines.append(format_table(grouped[section]))
            lines.append("")

    # Add any sections not in the predefined order
    for section, section_entries in grouped.items():
        if section not in section_order and section:
            lines.append(f"### {section}")
            lines.append("")
            lines.append(format_table(section_entries))
            lines.append("")

    # Add dimensionless quantities section
    lines.extend([
        "---",
        "",
        "## Dimensionless Quantities",
        "",
        "| Quantity                | Symbol | Unit     | Definition                    |",
        "|-------------------------|--------|----------|-------------------------------|",
        "| Mach number             | M      | -        | v / a                         |",
        "| Reynolds number         | Re     | -        | rho*V*L / mu                  |",
        "| Prandtl number          | Pr     | -        | mu*Cp / k                     |",
        "| Peclet number           | Pe     | -        | V*L / alpha                   |",
        "| Isentropic exponent     | gamma  | -        | Cp / Cv                       |",
        "| Equivalence ratio       | phi    | -        | (F/A) / (F/A)_stoich          |",
        "| Mixture fraction        | Z      | -        | Bilger definition             |",
        "| Discharge coefficient   | Cd     | -        | mdot_actual / mdot_ideal      |",
        "| Friction factor (Darcy) | f      | -        | dP / (L/D * rho*v^2/2)        |",
        "| Pressure ratio          | -      | -        | P / P0                        |",
        "| Diameter ratio          | beta   | -        | d / D                         |",
        "",
        "---",
        "",
        "## Common Reference Values",
        "",
        "| Quantity                | Value          | Unit   | Notes                    |",
        "|-------------------------|----------------|--------|--------------------------|",
        "| Standard pressure       | 101325         | Pa     | 1 atm                    |",
        "| Standard temperature    | 298.15         | K      | 25 C                     |",
        "| Standard gravity        | 9.80665        | m/s^2  | Used for Isp             |",
        "| Sea level air density   | ~1.225         | kg/m^3 | At 15 C, 101325 Pa       |",
        "",
        "---",
        "",
        "## Summary: Key Unit Conventions",
        "",
        "1. **Temperature**: Always Kelvin (K)",
        "2. **Pressure**: Always Pascal (Pa)",
        "3. **Molar mass**: g/mol (historical convention; convert to kg/mol for mass-basis calculations)",
        "4. **Thermodynamic properties**: Molar basis (J/mol, J/(mol*K)) in thermo functions",
        "5. **Compressible flow**: Mass basis (J/kg, J/(kg*K)) in flow solutions",
        "6. **Fractions**: mole fractions X (mol/mol), mass fractions Y (kg/kg)",
        "7. **Relative humidity**: Fraction (0-1), not percentage",
        "8. **Angles**: Radians (rad)",
        "",
    ])

    return "\n".join(lines)


def main() -> None:
    # Determine paths relative to script location
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    units_data_path = project_root / "include" / "units_data.h"
    output_path = project_root / "docs" / "UNITS.md"

    if not units_data_path.exists():
        raise FileNotFoundError(f"units_data.h not found at {units_data_path}")

    print(f"Parsing {units_data_path}...")
    entries = parse_units_data(units_data_path)
    print(f"Found {len(entries)} unit entries")

    print(f"Generating {output_path}...")
    content = generate_markdown(entries)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content)

    print(f"Done! Generated {output_path}")


if __name__ == "__main__":
    main()
