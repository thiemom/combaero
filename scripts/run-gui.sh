#!/bin/bash

# Exit on error
set -e

# Get the project root directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

# Ensure venv is activated or used directly
VENV_PYTHON="./.venv/bin/python"
if [ ! -f "$VENV_PYTHON" ]; then
    echo "Error: Virtual environment not found at ./.venv"
    echo "Please run ./scripts/bootstrap.sh first."
    exit 1
fi

# Function to cleanup background processes on exit
cleanup() {
    echo ""
    echo "Shutting down GUI services..."
    if [ ! -z "$BACKEND_PID" ]; then
        kill "$BACKEND_PID" 2>/dev/null || true
    fi
    exit 0
}

# Trap SIGINT (Ctrl+C) and SIGTERM
trap cleanup SIGINT SIGTERM

echo "Starting CombAero GUI..."

# 1. Start FastAPI Backend in background
echo "Starting Backend (FastAPI)..."
export PYTHONPATH="$PYTHONPATH:$PROJECT_ROOT"
$VENV_PYTHON -m uvicorn gui.backend.main:app --host 127.0.0.1 --port 8000 --reload &
BACKEND_PID=$!

# Wait a moment for backend to initialize
sleep 2

# 2. Start Vite Frontend
echo "Starting Frontend (Vite)..."
cd gui/frontend
npm run dev -- --port 5173 --host

# The script will wait here until Vite is stopped (Ctrl+C)
# Cleanup is handled by the trap
