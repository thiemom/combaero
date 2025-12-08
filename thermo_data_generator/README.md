# Thermo Data Generator

This directory contains tools to generate `thermo_transport_data.h` from thermodynamic and transport property databases.

## Quick Start

```bash
# Generate header from a Cantera YAML mechanism
python generate_thermo_data.py \
    --mechanism gri30_highT.yaml \
    --species "N2,O2,AR,CO2,H2O,CH4,C2H6,C3H8,H2,CO" \
    --output ../include/thermo_transport_data.h

# List available species in a mechanism
python generate_thermo_data.py --mechanism gri30_highT.yaml --list-species
```

## Scripts

| Script | Purpose |
|--------|---------|
| `generate_thermo_data.py` | **Main generator** - creates C++ header from YAML mechanism |
| `extract_thermo_transport_data.py` | Legacy: extracts to CSV (deprecated) |
| `create_cpp_header.py` | Legacy: CSV to header (deprecated) |

## Input Formats

### Cantera YAML Mechanism (current)

Combined thermo + transport data in Cantera's YAML format. These files contain:
- Species definitions with composition
- NASA polynomial coefficients (NASA7 format)
- Lennard-Jones transport parameters

Available mechanisms in this directory:
- `gri30_highT.yaml` - GRI-Mech 3.0 with extended temperature range
- `JetSurf2.yaml` - JetSurF 2.0 mechanism
- `aramco2.yaml` - AramcoMech 2.0
- `sandiego20161214.yaml` - San Diego mechanism

### Standalone Databases (future)

The generator is designed to support separate thermo and transport databases:

```bash
# Future usage (not yet implemented)
python generate_thermo_data.py \
    --thermo nasa9.dat \
    --transport transport.dat \
    --species "N2,O2,CH4"
```

This will enable:
- NASA9 polynomial format from NASA Glenn databases
- Mixing thermo from one source with transport from another

## NASA Polynomial Formats

### NASA7 (current)

Standard Chemkin format with 7 coefficients per temperature range:

```
Cp/R = a1 + a2*T + a3*T² + a4*T³ + a5*T⁴
H/RT = a1 + a2*T/2 + a3*T²/3 + a4*T³/4 + a5*T⁴/5 + a6/T
S/R  = a1*ln(T) + a2*T + a3*T²/2 + a4*T³/3 + a5*T⁴/4 + a7
```

### NASA9 (future)

Extended format with 9 coefficients, better accuracy at low temperatures:

```
Cp/R = a1/T² + a2/T + a3 + a4*T + a5*T² + a6*T³ + a7*T⁴
```

## Output

The generator produces `thermo_transport_data.h` containing:

- `USE_NASA9` - constexpr bool indicating polynomial format
- `species_names` - vector of species name strings
- `species_index` - map from name to index
- `molar_masses` - vector of molar masses [g/mol]
- `nasa_coeffs` - vector of NASA polynomial coefficients
- `transport_props` - vector of Lennard-Jones parameters
- `molecular_structures` - vector of atomic compositions (C, H, O, N)

## Testing

```bash
pytest test_create_cpp_header.py test_extract_thermo_transport_data.py
```
