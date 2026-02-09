.PHONY: help build test test-all clean style style-check style-fix install-deps

# Default target
help:
	@echo "CombAero - Make targets:"
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
	@echo "Running Python unit tests..."
	@cd python/tests && python -m pytest -v
	@echo "✓ Python tests complete"

test-validation:
	@echo "Running Cantera validation tests..."
	@cd cantera_validation_tests && poetry run pytest -v
	@echo "✓ Validation tests complete"

test-all: test test-python test-validation
	@echo "✓ All tests complete"

# Style checking
style-check-cpp:
	@echo "Checking C++ code style..."
	@./scripts/check-source-style.sh
	@echo "✓ C++ style check complete"

style-check-python:
	@echo "Checking Python code style..."
	@./scripts/check-python-style.sh
	@echo "✓ Python style check complete"

style-check: style-check-cpp style-check-python
	@echo "✓ All style checks complete"

style-fix:
	@echo "Auto-fixing Python code style..."
	@./scripts/check-python-style.sh --fix
	@echo "✓ Python style fixed"

# Dependencies
install-deps:
	@echo "Installing Python dependencies..."
	@cd cantera_validation_tests && poetry install
	@cd thermo_data_generator && poetry install
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
