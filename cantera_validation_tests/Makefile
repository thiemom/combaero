.PHONY: help sync test test-verbose test-parallel combustion-test equilibrium-test transport-test mixing-test

# Use uv run with --project . to ensure we target this workspace member
UV := uv run --project .

help:
	@echo "CombAero Cantera Validation Tests"
	@echo "Usage:"
	@echo "  make sync              - Synchronize dependencies with uv (run from root ideally)"
	@echo "  make test              - Run all validation tests"
	@echo "  make test-verbose      - Run all tests with verbose output"
	@echo "  make test-parallel     - Run tests in parallel (using pytest-xdist)"
	@echo ""
	@echo "Specific Validation Areas:"
	@echo "  make combustion-test   - Validate Tad and Complete Combustion"
	@echo "  make equilibrium-test  - Validate WGS and Reforming Equilibrium"
	@echo "  make transport-test    - Validate viscosity, conductivity, and Prandtl"
	@echo "  make mixing-test       - Validate stream enthalpy/mass balancing"

sync:
	uv sync

test:
	$(UV) pytest

test-verbose:
	$(UV) pytest -v -s

test-parallel:
	$(UV) pytest -n auto

combustion-test:
	$(UV) pytest test_combustion_validation.py -v

equilibrium-test:
	$(UV) pytest test_equilibrium_validation.py -v

transport-test:
	$(UV) pytest test_transport_validation.py -v

mixing-test:
	$(UV) pytest test_mixing_validation.py -v
