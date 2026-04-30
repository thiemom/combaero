# Development Workflow

This document outlines the development practices, testing procedures, and version synchronization rules for the CombAero project.

## Environment Setup

The root `.venv` is the single source of truth for all Python-based development tasks.

```bash
# Bootstrap the environment
./scripts/bootstrap.sh

# Verify the setup
./scripts/ensure_venv.py
```

## Testing

### C++ Tests
C++ tests use GoogleTest and are integrated with CTest.

```bash
cd build
ctest --output-on-failure
```

### Python Tests
Python tests use `pytest`.

```bash
python -m pytest python/tests
```

### Cantera Validation
We validate our physics against Cantera. These tests require the `CANTERA_DATA` environment variable.

```bash
CANTERA_DATA=$PWD/cantera_validation_tests python -m pytest cantera_validation_tests -v
```

## Version Synchronization

When updating tool versions or Python targets, ensure the following files remain synchronized:

1.  `.python-version`
2.  `ruff.toml` (`target-version`)
3.  `.github/workflows/ci.yml`
4.  `.github/workflows/validation-tests.yml`
5.  `.pre-commit-config.yaml`
6.  `Makefile` and the `.sh` scripts in `scripts/`

## Pre-commit Hooks

The project uses `pre-commit` to ensure code style and basic checks pass before every commit.

```bash
# Install hooks
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

## Documentation

### Units Reference
`docs/UNITS.md` is auto-generated from `include/units_data.h`. To update it:

```bash
python scripts/generate_units_md.py
```

## Environment Gotchas

### Here-Document Issues
Here-documents (`<<EOF`) consistently fail or get stuck in this execution environment.

**Avoid:**
```bash
cat << 'EOF' > file.txt
content
EOF
```

**Use Alternatives:**
- **`printf`**: `printf "line1\nline2\n" > file.txt`
- **`echo`**: `echo "line1" > file.txt` followed by `echo "line2" >> file.txt`
- **Python**: `python3 -c "open('file.txt', 'w').write('content')"`
- **Direct Tools**: Use the `write_to_file` tool provided by the environment.
