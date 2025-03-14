# Contributing to Thermo

Thank you for your interest in contributing to the Thermo project! This document provides guidelines and instructions for contributing.

## Code Style

- Use consistent indentation (4 spaces)
- Follow C++17 best practices
- Include descriptive comments for complex algorithms
- Add documentation for public APIs

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

Run tests before submitting a pull request:

```bash
cd build
make test
./tests/thermo_tests
./tests/test_ice_equation_accuracy
./tests/test_water_equation_accuracy
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
