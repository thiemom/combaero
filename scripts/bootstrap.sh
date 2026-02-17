#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${ROOT_DIR}/.venv"
PYTHON_BIN="${VENV_DIR}/bin/python"

cd "${ROOT_DIR}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    echo "Creating local virtual environment at ${VENV_DIR}"
    python3 -m venv "${VENV_DIR}"
fi

echo "Using interpreter: ${PYTHON_BIN}"
"${PYTHON_BIN}" -m pip install --upgrade pip
"${PYTHON_BIN}" -m pip install -e .[dev,examples] --no-build-isolation
"${PYTHON_BIN}" -m pip install ruff black isort flake8 mypy pre-commit build

echo "Installing pre-commit hooks"
"${PYTHON_BIN}" -m pre_commit install

echo "Bootstrap complete. Activate with: source .venv/bin/activate"
