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
from .runner import NetworkResult, _build_result_objects, _schema_maps
from .schemas import (
    NetworkGraphSchema,
    SolveResponse,
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

    node_results, element_results, edge_results = _build_result_objects(result, net)
    diag = getattr(solver, "_diagnostic_data", {})
    return result, node_results, element_results, edge_results, net, diag


@app.post("/solve", response_model=SolveResponse)
async def solve(schema: NetworkGraphSchema):
    # Hard async-level deadline: even if a single C++ Fanno/RK4 step blocks,
    # the event loop stays alive for health checks and subsequent requests.
    # Add a 10 s buffer on top of the solver's own soft timeout.
    soft_timeout = schema.solver_settings.timeout or 180.0
    hard_timeout = soft_timeout + 10.0
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, element_results, edge_results, _, diag = await asyncio.wait_for(
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
            convergence_history=diag.get("convergence_history", []),
            worst_residuals=diag.get("worst_residuals", []),
            solver_settings_used=diag.get("solver_settings_used", {}),
            lm_started_at_eval=diag.get("lm_started_at_eval"),
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=f"Solver exceeded the hard timeout of {hard_timeout:.0f}s. "
            "Try raising the solver timeout, simplifying the network, or a "
            "different init strategy.",
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
        raw, node_results, element_results, edge_results, net = await asyncio.wait_for(
            asyncio.to_thread(_solve_sync, schema),
            timeout=hard_timeout,
        )

        id_to_label, id_to_kind, id_to_method = _schema_maps(schema)
        result_obj = NetworkResult(
            raw=raw,
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
            net=net,
            schema=schema,
            id_to_label=id_to_label,
            id_to_kind=id_to_kind,
            id_to_method=id_to_method,
        )
        df = result_obj.to_dataframe()

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
