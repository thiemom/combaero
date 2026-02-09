.PHONY: install test test-verbose test-parallel test-combustion test-equilibrium test-transport test-mixing clean help

help:
	@echo "Cantera Validation Tests - Make targets:"
	@echo ""
	@echo "Setup:"
	@echo "  make install           - Install dependencies with Poetry"
	@echo ""
	@echo "Run tests:"
	@echo "  make test              - Run all tests"
	@echo "  make test-verbose      - Run tests with verbose output"
	@echo "  make test-parallel     - Run tests in parallel"
	@echo "  make test-combustion   - Run combustion tests only"
	@echo "  make test-equilibrium  - Run equilibrium tests only"
	@echo "  make test-transport    - Run transport tests only"
	@echo "  make test-mixing       - Run mixing tests only"
	@echo ""
	@echo "Maintenance:"
	@echo "  make clean             - Clean cache and temporary files"

install:
	poetry install

test:
	poetry run pytest

test-verbose:
	poetry run pytest -v -s

test-parallel:
	poetry run pytest -n auto

test-combustion:
	@echo "Running combustion validation tests..."
	poetry run pytest test_combustion_validation.py -v

test-equilibrium:
	@echo "Running equilibrium validation tests..."
	poetry run pytest test_equilibrium_validation.py -v

test-transport:
	@echo "Running transport validation tests..."
	poetry run pytest test_transport_validation.py -v

test-mixing:
	@echo "Running mixing validation tests..."
	poetry run pytest test_mixing_validation.py -v

clean:
	rm -rf .pytest_cache __pycache__
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
