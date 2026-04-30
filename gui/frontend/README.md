# CombAero GUI Frontend

This directory contains the React + TypeScript + Vite frontend for the CombAero network designer.

## Tech Stack
- **Framework**: React 18+
- **Build Tool**: Vite
- **Language**: TypeScript
- **State Management**: React Context / Hooks
- **Styling**: CSS Modules (Vanilla CSS)
- **Linting/Formatting**: Biome

## Getting Started

1. **Install Dependencies**:
   ```bash
   pnpm install
   ```

2. **Run Development Server**:
   ```bash
   pnpm run dev
   ```

3. **Build for Production**:
   ```bash
   pnpm run build
   ```

## Development Rules

- **Biome**: We use Biome for linting and formatting. Run `./scripts/check-gui-style.sh` from the repo root before committing.
- **Components**: Keep components functional and focused. Use CSS modules for scoped styling.
- **API**: The frontend communicates with the FastAPI backend (located in `gui/backend/`). Ensure schemas match `gui/backend/schemas.py`.

## Agent Guidance
- If adding a new network element, you must update the node palette in `src/components/Palette.tsx` and the inspector logic in `src/components/Inspector.tsx`.
- Refer to [docs/GUI_TECHNICAL.md](../../docs/GUI_TECHNICAL.md) for the JSON schemas expected by the backend.
