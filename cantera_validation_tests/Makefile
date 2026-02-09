.PHONY: install test test-verbose test-parallel clean help

help:
	@echo "Cantera Validation Tests - Make targets:"
	@echo "  make install        - Install dependencies with Poetry"
	@echo "  make test           - Run all tests"
	@echo "  make test-verbose   - Run tests with verbose output"
	@echo "  make test-parallel  - Run tests in parallel"
	@echo "  make clean          - Clean cache and temporary files"

install:
	poetry install

test:
	poetry run pytest

test-verbose:
	poetry run pytest -v

test-parallel:
	poetry run pytest -n auto

clean:
	rm -rf .pytest_cache __pycache__
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
