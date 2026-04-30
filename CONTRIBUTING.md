# Contributing to CombAero

Thank you for your interest in contributing to CombAero! This project welcomes contributions in physics modeling, numerical methods, and UI/UX.

## Core Rules

1. **ASCII Only**: No non-ASCII characters in source files (allowed in comments/strings).
2. **Line Comments**: Use `//` line comments for C++.
3. **Type Safety**: Mandatory type hints for Python.
4. **Jacobian Mandate**: Physics correlations must expose `(f, J)` analytical derivatives where applicable.

For detailed coding styles and hard rules, see [GEMINI.md](GEMINI.md).

## Development Workflow

We use a standard fork-and-pull-request workflow. For detailed instructions on environment setup, building, and running tests, please refer to:

- [Building CombAero](docs/BUILDING.md)
- [Development & Testing Guide](docs/DEVELOPMENT.md)

## Pull Request Process

1. **Update Documentation**: If your change modifies a public API or a physical correlation, update the relevant markdown file in `docs/`.
2. **Automated Checks**: Ensure your code passes all CI checks (style, units sync, and unit tests).
3. **Review**: All PRs require review from at least one maintainer.

## Reporting Issues

Use GitHub Issues to report bugs or request features. Please include:
- A clear description of the problem or feature.
- Reproduction steps (for bugs) or motivation (for features).
- Environment details if relevant.
