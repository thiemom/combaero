.PHONY: help bootstrap venv-check build test test-all clean style style-check style-fix install-deps pre-build ci

ROOT_DIR := $(abspath .)
PYTHON := $(ROOT_DIR)/.venv/bin/python
PIP := $(PYTHON) -m pip

# Default target
help:
	@echo "CombAero - Make targets:"
	@echo ""
	@echo "Environment:"
	@echo "  make bootstrap          - Create/update local .venv and install tools"
	@echo "  make venv-check         - Verify commands use local .venv"
	@echo ""
	@echo "Build targets:"
	@echo "  make build              - Build C++ library and tests (Release)"
	@echo "  make build-debug        - Build with debug symbols"
	@echo "  make clean              - Clean build artifacts"
	@echo ""
	@echo "Test targets:"
	@echo "  make test               - Run C++ tests only"
	@echo "  make test-python        - Run Python unit tests"
	@echo "  make test-validation    - Run Cantera validation tests"
	@echo "  make test-all           - Run all tests (C++, Python, validation)"
	@echo ""
	@echo "Style checking:"
	@echo "  make style-check        - Check code style (C++ and Python)"
	@echo "  make style-check-cpp    - Check C++ style only"
	@echo "  make style-check-python - Check Python style only"
	@echo "  make style-fix          - Auto-fix Python style issues"
	@echo ""
	@echo "Dependencies:"
	@echo "  make install-deps       - Install Python dependencies"
	@echo ""
	@echo "Pre-build checks:"
	@echo "  make pre-build          - Run style checks before building"
	@echo "  make ci                 - Full CI workflow (style + build + test)"

bootstrap:
	@./scripts/bootstrap.sh

venv-check:
	@$(PYTHON) scripts/ensure_venv.py

# Build targets
build:
	@echo "Building CombAero (Release)..."
	@mkdir -p build
	@cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . --config Release
	@echo "✓ Build complete"

build-debug:
	@echo "Building CombAero (Debug)..."
	@mkdir -p build
	@cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && cmake --build . --config Debug
	@echo "✓ Debug build complete"

# Test targets
test:
	@echo "Running C++ tests..."
	@cd build && ctest --output-on-failure
	@echo "✓ C++ tests complete"

test-python:
	@$(PYTHON) scripts/ensure_venv.py --quiet
	@echo "Running Python unit tests..."
	@$(PYTHON) -m pytest python/tests -v
	@echo "✓ Python tests complete"

test-validation:
	@$(PYTHON) scripts/ensure_venv.py --quiet
	@echo "Running Cantera validation tests..."
	@CANTERA_DATA=$(ROOT_DIR)/cantera_validation_tests $(PYTHON) -m pytest cantera_validation_tests -v
	@echo "✓ Validation tests complete"

test-all: test test-python test-validation
	@echo "✓ All tests complete"

# Style checking
style-check-cpp:
	@echo "Checking C++ code style..."
	@./scripts/check-source-style.sh
	@echo "✓ C++ style check complete"

style-check-python:
	@$(PYTHON) scripts/ensure_venv.py --quiet
	@echo "Checking Python code style..."
	@./scripts/check-python-style.sh
	@echo "✓ Python style check complete"

style-check: style-check-cpp style-check-python
	@echo "✓ All style checks complete"

style-fix:
	@$(PYTHON) scripts/ensure_venv.py --quiet
	@echo "Auto-fixing Python code style..."
	@./scripts/check-python-style.sh --fix
	@echo "✓ Python style fixed"

# Dependencies
install-deps:
	@$(PYTHON) scripts/ensure_venv.py --quiet
	@echo "Installing Python dependencies..."
	@$(PIP) install -e .[dev,examples] --no-build-isolation
	@$(PIP) install ruff black isort flake8 mypy pre-commit build
	@echo "✓ Dependencies installed"

# Pre-build workflow
pre-build: style-check
	@echo "✓ Pre-build checks passed"

# CI workflow
ci: style-check build test-all
	@echo "✓ CI workflow complete"

# Clean
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf build
	@find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@find . -type f -name "*.pyc" -delete
	@find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	@echo "✓ Clean complete"
