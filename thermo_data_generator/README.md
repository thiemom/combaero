# Thermo Data Generator

Tools to extract thermodynamic and transport data from external databases and
generate `thermo_transport_data.h` for the C++ library.

## Quick Start

### Step 1: Obtain Source Data (not included)

Source data files are **not included** in this repository to avoid licensing issues.
You must obtain them yourself:

**Cantera YAML mechanisms** (NASA-7 + transport):
- Download from [Cantera](https://cantera.org/examples/input-files.html)
- Or convert from Chemkin format using `ck2yaml`

**NASA CEA data** (NASA-9, high temperature):
- Use the [NASA CEA web interface](https://cearun.grc.nasa.gov/)
- Select species, generate output, save as text file

### Step 2: Extract Species Data

```bash
# Extract from YAML only (NASA-7 + transport)
python extract_species_data.py --yaml mechanism.yaml -o species.json

# Extract from CEA only (NASA-9, up to 20000 K)
python extract_species_data.py --cea cea_output.txt -o species.json

# Merge both sources (recommended: NASA-9 thermo + YAML transport)
python extract_species_data.py \
    --yaml mechanism.yaml \
    --cea cea_output.txt \
    -s O2 N2 AR H2O CO2 CH4 H2 CO \
    -o species.json
```

### Step 3: Generate C++ Header

```bash
python generate_thermo_data.py \
    --mechanism mechanism.yaml \
    --species "N2,O2,AR,CO2,H2O,CH4,C2H6,C3H8,H2,CO" \
    --output ../include/thermo_transport_data.h
```

## Scripts

| Script | Purpose |
|--------|---------|
| `extract_species_data.py` | **Unified extractor** - extracts from YAML and/or CEA to JSON |
| `generate_thermo_data.py` | **Header generator** - creates C++ header from YAML mechanism |

### Legacy Scripts (deprecated)

| Script | Purpose |
|--------|---------|
| `extract_thermo_transport_data.py` | Standalone YAML extractor |
| `extract_nasa_cea_data.py` | Standalone CEA extractor |
| `create_cpp_header.py` | CSV to header converter |

## Data Sources

### Cantera YAML Mechanisms

Combined thermo + transport data in Cantera's YAML format:
- NASA-7 polynomial coefficients (typically 200-3500 K)
- Lennard-Jones transport parameters

Recommended sources:
- **GRI-Mech 3.0** - natural gas combustion
- **JetSurF 2.0** - jet fuel surrogate
- **San Diego mechanism** - hydrocarbon combustion

**Note: YAML 1.1 Boolean Parsing**

PyYAML parses certain species names as booleans (`NO` → `False`).
The generator handles this automatically.

### NASA CEA Database

NASA-9 polynomials with extended temperature range (200-20000 K):
- Higher accuracy at extreme temperatures
- No transport data (must be obtained separately)

## NASA Polynomial Formats

### NASA-7 (from YAML)

7 coefficients per temperature range (typically 2 ranges):

```
Cp/R = a1 + a2*T + a3*T² + a4*T³ + a5*T⁴
H/RT = a1 + a2*T/2 + a3*T²/3 + a4*T³/4 + a5*T⁴/5 + a6/T
S/R  = a1*ln(T) + a2*T + a3*T²/2 + a4*T³/3 + a5*T⁴/4 + a7
```

### NASA-9 (from CEA)

10 coefficients per interval (typically 2-3 intervals):

```
Cp/R = a1/T² + a2/T + a3 + a4*T + a5*T² + a6*T³ + a7*T⁴
```

Coefficients a8, a9 are integration constants for H and S.

## Output Formats

### JSON (from extract_species_data.py)

```json
{
  "sources": ["mechanism.yaml", "cea_output.txt"],
  "species": [
    {
      "name": "N2",
      "molar_mass": 28.014,
      "thermo_nasa7": {...},
      "thermo_nasa9": {"intervals": [...]},
      "transport": {"geometry": "linear", ...}
    }
  ]
}
```

### C++ Header (from generate_thermo_data.py)

- `species_names` - vector of species name strings
- `molar_masses` - vector of molar masses [g/mol]
- `nasa_coeffs` - vector of NASA polynomial coefficients
- `transport_props` - vector of Lennard-Jones parameters

The header uses format-specific structs with a common alias:

**NASA-7**: `NASA7_Coeffs` with fixed low/high temperature ranges (2 intervals)

**NASA-9**: `NASA9_Coeffs` with variable `intervals` vector (1-3 intervals per species)

## Testing

```bash
pytest
```

## Licensing Note

**Source data files (YAML mechanisms, CEA output) are not included** in this
repository. Many thermodynamic databases have specific licensing terms:

- Some mechanisms incorporate Burcat's database (restricted use)
- NASA CEA data is public domain
- Cantera mechanism files have varying licenses

Always verify the license of your source data before use.
