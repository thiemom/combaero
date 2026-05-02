# CombAero Network Designer GUI

**Interactive design and diagnostic environment for fluid-thermal networks.**

The CombAero GUI allows for rapid construction of complex aerospace networks using a drag-and-drop interface powered by React Flow and a high-performance C++ backend.

## Quick Start

### 1. Install

Create a virtual environment (Python 3.12+; we recommend [uv](https://docs.astral.sh/uv/)) and install:

```bash
uv pip install combaero-gui
```

### 2. Launch

```bash
combaero-gui
```

Open http://127.0.0.1:8000 in your browser.

## Key Workflows

- **Build**: Drag nodes (Plenums, Boundaries) and elements (Orifices, Channels) from the sidebar.
- **Connect**: Click and drag from node handles to element ports.
- **Solve**: Click the **Solve** button in the header. Successful solves turn node borders green.
- **Inspect**: Click on any node or element to view real-time diagnostics in the Inspector sidebar.

---

## Agent Reference (UI Automation Map)

For agents with screen control (multimodal LLMs), the following selectors and IDs are provided to assist in navigation:

| UI Component | Selector / ID | Description |
| :--- | :--- | :--- |
| **Solve Button** | `#solve-btn` | Triggers the network solve. |
| **Inspector Sidebar** | `.inspector-panel` | Properties for selected node/element. |
| **Network Canvas** | `#rf-canvas` | The main React Flow workspace. |
| **Node (General)** | `.react-flow__node` | Any node on the canvas. |
| **Solved Indicator** | `.solved-status-ok` | Appears on nodes after a successful solve. |
| **Error Console** | `#error-output` | Displays solver convergence failures. |

> [!TIP]
> **Automation Hint**: To add a node via script, use the `gui.backend` API or generate a network JSON and load it via `File -> Import`.

---

## Related

- **[combaero](https://pypi.org/project/combaero/)** — the C++ physics core and Python API powering this GUI.
- **[Source & docs](https://github.com/thiemom/combaero)** — monorepo with full API reference and technical documentation.
