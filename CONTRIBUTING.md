# Contributing to CombAero

Thank you for your interest in contributing to the CombAero project! This document provides guidelines and instructions for contributing.

## Code Style

- Use consistent indentation (4 spaces)
- Follow C++17 best practices
- Include descriptive comments for complex algorithms
- Add documentation for public APIs
- **ASCII only**: No non-ASCII characters in code (allowed in comments and string literals)
- **Line comments only**: Use `//` comments, not `/* */` block comments

These rules are enforced by `scripts/check-source-style.sh` (C++) and `scripts/check-python-style.sh` (Python/ruff) which run as pre-commit hooks and in CI.

## Development Workflow

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature-name`)
3. Make your changes
4. Run tests to ensure they pass
5. Commit your changes (`git commit -m 'Add some feature'`)
6. Push to the branch (`git push origin feature/your-feature-name`)
7. Create a new Pull Request

## Testing

All new features should include appropriate tests:

- Unit tests for individual functions
- Integration tests for complex features
- Accuracy tests for numerical algorithms

Run the C++ test suite before submitting a pull request:

```bash
mkdir -p build
cd build
cmake ..
cmake --build .

# All C++ tests
ctest --output-on-failure

# Main test executable directly (bypassing CTest)
./tests/combaero_tests

# Accuracy tests for saturation vapor pressure
./tests/test_ice_equation_accuracy
./tests/test_water_equation_accuracy
```

If you modify the Python bindings or C++ source, rebuild and run the full Python test suite:

```bash
./scripts/bootstrap.sh
./.venv/bin/python scripts/ensure_venv.py
# Rebuild the C++ extension (required after any src/ or include/ change)
./.venv/bin/pip install -e . --no-build-isolation
./.venv/bin/python -m pytest python/tests/
```

Pre-commit hooks run automatically on `git commit`. To run them manually:

```bash
./.venv/bin/pre-commit run --all-files
```

## Pull Request Process

1. Update the README.md with details of changes if applicable
2. Update the documentation if you're changing public APIs
3. The PR should work on the GitHub Actions CI pipeline
4. PRs require review from at least one maintainer

## Reporting Bugs

When reporting bugs, please include:

- A clear description of the issue
- Steps to reproduce the behavior
- Expected behavior
- Actual behavior
- Environment details (OS, compiler version, etc.)

## Feature Requests

Feature requests are welcome. Please provide:

- A clear description of the proposed feature
- The motivation for the feature
- Any potential implementation details you have in mind
