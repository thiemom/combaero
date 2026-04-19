#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "${ROOT_DIR}"

# Check for uv
if ! command -v uv &> /dev/null; then
    # If not in path, check if it was installed by astral.sh in the default location
    if [[ -x "$HOME/.local/bin/uv" ]]; then
        export PATH="$HOME/.local/bin:$PATH"
    else
        echo "uv not found. Installing uv..."
        curl -LsSf https://astral.sh/uv/install.sh | sh
        export PATH="$HOME/.local/bin:$PATH"
    fi
fi

echo "Initializing environment with uv..."
uv sync --all-extras --all-groups

echo "Installing pre-commit hooks"
uv run pre-commit install

echo "Bootstrap complete. The environment is managed by uv at .venv"
echo "To activate: source .venv/bin/activate"
