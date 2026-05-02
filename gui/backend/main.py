import asyncio
import io
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

import combaero as cb
from combaero.network import NetworkSolver

from .graph_builder import build_network_from_schema
from .schemas import (
    ElementResult,
    NetworkGraphSchema,
    NodeResult,
    SolveResponse,
    StateResult,
)

_FRONTEND_DIST = Path(__file__).parent.parent / "frontend" / "dist"

app = FastAPI(title="CombAero Network GUI API")

# Configure CORS for local development (Vite typically runs on 5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/favicon.svg")
async def favicon_svg():
    return FileResponse(_FRONTEND_DIST / "favicon.svg")


@app.get("/icons.svg")
async def icons_svg():
    return FileResponse(_FRONTEND_DIST / "icons.svg")


@app.get("/")
async def root():
    if _FRONTEND_DIST.is_dir():
        return FileResponse(_FRONTEND_DIST / "index.html")
    return {
        "message": "CombAero Network GUI API is running",
        "docs": "/docs",
        "health": "/health",
        "metadata": ["/metadata/species", "/metadata/presets", "/metadata/materials"],
    }


@app.get("/metadata/species")
async def get_species():
    """Returns available species metadata from the C++ core."""

    return {
        "names": cb.species.names,
        "molar_masses": cb.species._default_locator.molar_masses.tolist(),
    }


@app.get("/metadata/presets")
async def get_presets():
    """Returns available composition presets."""
    return {
        "sources": ["dry_air", "humid_air", "fuel", "custom"],
        "modes": ["mole", "mass"],
    }


@app.get("/metadata/materials")
async def get_materials():
    """Returns available wall materials from the database."""
    return {
        "names": cb.list_materials(),
    }


# Module-level persistence for warm-starting the solver across requests.
_continuation_state: dict | None = None


def _solve_sync(schema: NetworkGraphSchema):
    """
    Synchronous helper to build and solve the network.
    Runs in a separate thread to avoid blocking the FastAPI event loop.
    """
    global _continuation_state
    net = build_network_from_schema(schema)
    solver = NetworkSolver(net)
    x0 = None
    if schema.solver_settings.init_strategy == "continuation" and _continuation_state:
        # Check topology fingerprint
        solver._build_x0()
        if _continuation_state["unknown_names"] != list(solver.unknown_names):
            raise ValueError(
                "Topology mismatch: the current network unknowns do not match "
                "the saved continuation state. Please run a 'Default' solve first."
            )

        # Safety Reset: Check if initial guesses have changed
        current_guesses = {n.id: n.data.get("initial_guess", {}) for n in schema.nodes}
        if _continuation_state.get("initial_guesses") != current_guesses:
            raise ValueError(
                "Initial guess mismatch: user overrides have changed since the last solve. "
                "Please run a 'Default' solve to incorporate these changes."
            )

        x0 = _continuation_state["x"]

    result = solver.solve(
        method=schema.solver_settings.method,
        init_strategy=schema.solver_settings.init_strategy,
        timeout=schema.solver_settings.timeout,
        x0=x0,
    )
    success = bool(result.get("__success__", False))

    if success:
        _continuation_state = {
            "x": result["__x_solution__"],
            "unknown_names": result["__unknown_names__"],
            "initial_guesses": {n.id: n.data.get("initial_guess", {}) for n in schema.nodes},
        }

    node_results = {}
    for node_id in net.nodes:
        # Pull all node diagnostics and state variables directly from sol_dict
        # flat keys are like node_id.T, node_id.P, node_id.mach, node_id.rho, etc.
        prefix = f"{node_id}."
        node_vals = {k[len(prefix) :]: v for k, v in result.items() if k.startswith(prefix)}

        # Handle composition lists (Y and X)
        y_vals = []
        x_vals = []
        i = 0
        while f"Y[{i}]" in node_vals:
            y_vals.append(float(node_vals.pop(f"Y[{i}]")))
            if f"X[{i}]" in node_vals:
                x_vals.append(float(node_vals.pop(f"X[{i}]")))
            i += 1

        state_res = StateResult(
            T=float(node_vals.pop("T", 300.0)),
            P=float(node_vals.pop("P", 101325.0)),
            Pt=node_vals.pop("Pt", None),
            Tt=node_vals.pop("Tt", None),
            Y=y_vals,
            X=x_vals or None,
            **node_vals,  # All remaining diagnostics (rho, mach, mu, k, etc.)
        )
        node_results[node_id] = NodeResult(state=state_res, success=success)

    element_results = {}
    elem_diags = result.get("__element_diag__", {})
    for elem_id in net.elements:
        m_dot = float(result.get(f"{elem_id}.m_dot", 0.0))
        # Pull element diagnostics from the structured dict aggregated by solve()
        diag = elem_diags.get(elem_id, {}).copy()
        diag.pop("m_dot", None)
        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success, **diag)

    edge_results = {}

    def _safe_float(v):
        try:
            if isinstance(v, (list, tuple, dict)):
                return v
            return float(v)
        except (TypeError, ValueError):
            return v

    for edge_id in net.walls:
        edge_keys = {
            k.split(".", 1)[1]: _safe_float(v)
            for k, v in result.items()
            if k.startswith(f"{edge_id}.")
        }
        if edge_keys:
            edge_results[edge_id] = edge_keys

    return result, node_results, element_results, edge_results, net


@app.post("/solve", response_model=SolveResponse)
async def solve(schema: NetworkGraphSchema):
    # Hard async-level deadline: even if a single C++ Fanno/RK4 step blocks,
    # the event loop stays alive for health checks and subsequent requests.
    # Add a 10 s buffer on top of the solver's own soft timeout.
    soft_timeout = schema.solver_settings.timeout or 180.0
    hard_timeout = soft_timeout + 10.0
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, element_results, edge_results, _ = await asyncio.wait_for(
            asyncio.to_thread(_solve_sync, schema),
            timeout=hard_timeout,
        )

        return SolveResponse(
            success=result.get("__success__", False),
            message=result.get("__message__", "Solve completed"),
            final_norm=result.get("__final_norm__"),
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=f"Solver exceeded the hard timeout of {hard_timeout:.0f}s. "
            "Try a shorter timeout, a simpler network, or the homotopy init strategy.",
        ) from None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.post("/export")
async def export_results(schema: NetworkGraphSchema):
    soft_timeout = schema.solver_settings.timeout or 180.0
    hard_timeout = soft_timeout + 10.0
    try:
        # Solve quickly to get states
        _, node_results, element_results, edge_results, _ = await asyncio.wait_for(
            asyncio.to_thread(_solve_sync, schema),
            timeout=hard_timeout,
        )

        # Convert to DataFrame
        import pandas as pd

        # Build unified id_to_label map from nodes, edges, and auto-links
        id_to_label = {}
        # Map each schema node/edge id to its declared ``type`` so the CSV can
        # preserve specific-kind information (plenum vs combustor vs
        # discrete_loss, …) needed for later round-trips or sweep workflows.
        id_to_kind: dict[str, str] = {}
        for n in schema.nodes:
            id_to_label[n.id] = n.data.get("label") or ""
            id_to_kind[n.id] = n.type
        for e in schema.edges:
            lbl = (e.data or {}).get("label") or ""
            id_to_label[e.id] = lbl
            id_to_kind[e.id] = (e.data or {}).get("type") or "flow"
            # Also map solver auto-link IDs
            auto_id = f"__auto_link__{e.source}__{e.target}"
            id_to_label[auto_id] = lbl
            id_to_kind[auto_id] = "lossless_connection"

        _BOUNDARY_KINDS = {"mass_boundary", "pressure_boundary"}

        # Per-node combustion method (empty string for non-combustor rows).
        id_to_method: dict[str, str] = {}
        for n in schema.nodes:
            if n.type == "combustor":
                id_to_method[n.id] = n.data.get("method", "complete")

        # Species labels for composition column headers (``Y[N2]`` not ``Y[0]``).
        import combaero as cb

        species_labels: list[str] = list(cb.species.names)

        data = []

        # 1. Process Nodes
        for node_id, res in node_results.items():
            # Include all fields from StateResult (rho, mach, k, etc.)
            state_data = res.state.model_dump()
            # Handle Y/X list expansion
            Y = state_data.pop("Y", [])
            X = state_data.pop("X", None)

            kind = id_to_kind.get(node_id, "")
            base_data = {
                "type": "node",
                "kind": kind,
                "is_boundary": kind in _BOUNDARY_KINDS,
                "combustion_method": id_to_method.get(node_id, ""),
                "id": node_id,
                "label": id_to_label.get(node_id, ""),
                **state_data,
            }
            # Flatten species fractions with species-name indexing.
            for i, y_val in enumerate(Y):
                name = species_labels[i] if i < len(species_labels) else str(i)
                base_data[f"Y[{name}]"] = y_val
            if X:
                for i, x_val in enumerate(X):
                    name = species_labels[i] if i < len(species_labels) else str(i)
                    base_data[f"X[{name}]"] = x_val

            data.append(base_data)

        # 2. Process Flow Elements
        for elem_id, res in element_results.items():
            # Include all analytic/diagnostic fields (Re, Nu, h, q_dot, f, etc.)
            elem_data = res.model_dump()
            data.append(
                {
                    "type": "element",
                    "kind": id_to_kind.get(elem_id, ""),
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": elem_id,
                    "label": id_to_label.get(elem_id, ""),
                    **elem_data,
                }
            )

        # 3. Process Thermal Walls (edge_results)
        for wall_id, res in edge_results.items():
            # thermal wall results are already a dict from _solve_sync
            data.append(
                {
                    "type": "thermal_wall",
                    "kind": "thermal_wall",
                    "is_boundary": False,
                    "combustion_method": "",
                    "id": wall_id,
                    "label": id_to_label.get(wall_id, ""),
                    **res,
                }
            )

        df = pd.DataFrame(data)

        # Rename columns with unit annotations (e.g. ``P`` -> ``P [Pa]``).
        # Meta columns (``id``, ``type``, ``kind``, ``is_boundary``, ...) are
        # passed through unchanged; unknown columns also remain untouched.
        from gui.backend.units import label_with_unit

        df = df.rename(columns={col: label_with_unit(col) for col in df.columns})

        # Stream as CSV
        stream = io.StringIO()
        df.to_csv(stream, index=False)

        return StreamingResponse(
            iter([stream.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=combaero_results.csv"},
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=f"Export solve exceeded the hard timeout of {hard_timeout:.0f}s.",
        ) from None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.get("/solver/continuation_available")
async def get_continuation_available():
    """Returns whether a converged state exists for continuation."""
    return {"available": _continuation_state is not None}


@app.get("/health")
async def health():
    return {"status": "ok"}




if (_FRONTEND_DIST / "assets").is_dir():
    app.mount("/assets", StaticFiles(directory=str(_FRONTEND_DIST / "assets")), name="assets")


@app.get("/{path:path}")
async def serve_spa_route(path: str):
    if _FRONTEND_DIST.is_dir():
        return FileResponse(_FRONTEND_DIST / "index.html")
    raise HTTPException(status_code=404, detail=f"/{path} not found")


def run_server(host: str = "127.0.0.1", port: int = 8000) -> None:
    import uvicorn

    uvicorn.run(app, host=host, port=port)
