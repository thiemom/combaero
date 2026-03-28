# CombAero Network Designer GUI

This is a prototype graphical user interface for building and solving fluid network simulations using the `combaero` physics core.

## Structure
- `backend/`: FastAPI server that bridges React Flow JSON with C++ CombAero objects.
- `frontend/`: React + TypeScript + React Flow interactive canvas.

## Getting Started

### 1. Install Dependencies
Ensure you have the GUI optional dependencies installed in your Python environment:
```bash
uv pip install -e ".[gui]"
```
Install frontend dependencies:
```bash
cd frontend && npm install
```

### 2. Run the Servers

**Start the Backend:**
```bash
# From the project root
export PYTHONPATH=$PYTHONPATH:.
python3 -m uvicorn gui.backend.main:app --reload
```

**Start the Frontend:**
```bash
# In a new terminal
cd gui/frontend
npm run dev
```

## Features
- Drag-and-drop network construction.
- Live telemetry indicators on solved nodes (Green borders).
- Detailed result inspection in the sidebar.
- Export solved results to CSV/Pandas DataFrame.
