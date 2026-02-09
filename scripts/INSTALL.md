# Installing Script Dependencies

The Python style checker (`check-python-style.sh`) requires several Python tools.

## Option 1: Install via pip (Global)

```bash
pip install black isort flake8 mypy
```

## Option 2: Install via Poetry (Project-specific)

### For cantera_validation_tests

```bash
cd cantera_validation_tests
poetry add --group dev black isort flake8 mypy
```

### For thermo_data_generator

```bash
cd thermo_data_generator
poetry add --group dev black isort flake8 mypy
```

### For python bindings

Create a `pyproject.toml` in the `python/` directory:

```bash
cd python
poetry init --no-interaction
poetry add --group dev black isort flake8 mypy
```

## Option 3: Install via uv (Fast, Modern)

### Global installation

```bash
uv pip install black isort flake8 mypy
```

### Project-specific (using uv's virtual environment)

```bash
# Create a venv for development tools
uv venv .venv-dev
source .venv-dev/bin/activate  # or .venv-dev\Scripts\activate on Windows

# Install tools
uv pip install black isort flake8 mypy
```

## Option 4: Use Existing Poetry Environments

The validation tests already have Poetry set up:

```bash
cd cantera_validation_tests
poetry install  # Installs existing dependencies
poetry add --group dev black isort flake8 mypy  # Add style tools
poetry run black --version  # Verify installation
```

## Recommended Setup

For development, add these to the **root** `pyproject.toml` (if you create one):

```toml
[tool.poetry.group.dev.dependencies]
black = "^24.0.0"
isort = "^5.13.0"
flake8 = "^7.0.0"
mypy = "^1.8.0"

[tool.black]
line-length = 100
target-version = ['py311']

[tool.isort]
profile = "black"
line_length = 100

[tool.flake8]
max-line-length = 100
extend-ignore = ["E203", "W503"]
exclude = [".venv", "__pycache__", ".pytest_cache", "build", "dist"]
```

## Verify Installation

```bash
# Check if tools are available
black --version
isort --version
flake8 --version
mypy --version

# Run the style checker
./scripts/check-python-style.sh
```

## CI/CD Integration

For GitHub Actions or similar:

```yaml
- name: Install Python style tools
  run: |
    pip install black isort flake8 mypy
    # or: uv pip install black isort flake8 mypy

- name: Check Python style
  run: ./scripts/check-python-style.sh
```

## Optional: Pre-commit Hooks

Install pre-commit and add to `.pre-commit-config.yaml`:

```yaml
repos:
  - repo: https://github.com/psf/black
    rev: 24.1.0
    hooks:
      - id: black
        args: [--line-length=100]

  - repo: https://github.com/pycqa/isort
    rev: 5.13.2
    hooks:
      - id: isort
        args: [--profile=black, --line-length=100]

  - repo: https://github.com/pycqa/flake8
    rev: 7.0.0
    hooks:
      - id: flake8
        args: [--max-line-length=100, --extend-ignore=E203,W503]
```

Then:

```bash
pip install pre-commit
pre-commit install
```

## Quick Start (Recommended)

```bash
# Using uv (fastest)
uv pip install black isort flake8 mypy

# Run style checker
./scripts/check-python-style.sh

# Auto-fix issues
./scripts/check-python-style.sh --fix
```
