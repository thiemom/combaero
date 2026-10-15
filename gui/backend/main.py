import asyncio
import io
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

import combaero as cb
from combaero.network import NetworkSolver

from . import second_law_view, validity_view
from .graph_builder import build_network_from_schema
from .messages import humanize
from .resolve_preview import resolve_preview
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
    if schema.solver_settings.incompressible_mach_limit is not None:
        solver.incompressible_mach_limit = schema.solver_settings.incompressible_mach_limit
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
        result, node_results, element_results, edge_results, net, diag = await asyncio.wait_for(
            asyncio.to_thread(_solve_sync, schema),
            timeout=hard_timeout,
        )

        gui_ids = {n.id for n in schema.nodes}
        return SolveResponse(
            success=result.get("__success__", False),
            message=humanize(result.get("__message__", "Solve completed"), schema),
            final_norm=result.get("__final_norm__"),
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
            convergence_history=diag.get("convergence_history", []),
            worst_residuals=diag.get("worst_residuals", []),
            solver_settings_used=diag.get("solver_settings_used", {}),
            lm_started_at_eval=diag.get("lm_started_at_eval"),
            second_law_failures=second_law_view.failures(result, net, gui_ids),
            incompressible_warnings=validity_view.warnings(result, net, gui_ids),
            ledger=second_law_view.ledger(result, net, gui_ids),
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=f"Solver exceeded the hard timeout of {hard_timeout:.0f}s. "
            "Try raising the solver timeout, simplifying the network, or a "
            "different init strategy.",
        ) from None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=humanize(str(e), schema)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=humanize(str(e), schema)) from e


def _resolve_preview_sync(schema: NetworkGraphSchema) -> dict:
    net = build_network_from_schema(schema)
    report = resolve_preview(net)
    gui_ids = {n.id for n in schema.nodes}
    items: dict[str, dict] = {}
    for item_id, res in report.items():
        if item_id not in gui_ids:
            continue  # synthetic builder items (auto-links, tee steps) are not drawn
        items[item_id] = {
            "status": res.status,
            "reason": res.reason,
            "fields": {k: {"value": fr.value, "source": fr.source} for k, fr in res.fields.items()},
        }
    return {"items": items}


@app.post("/resolve_preview")
async def resolve_preview_endpoint(schema: NetworkGraphSchema):
    """Resolve geometry inheritance WITHOUT solving (#559).

    Per drawn element/node, returns each geometry field's resolved value and
    provenance ("user" / "inherited from ..." / "default"), and whether the item
    is unresolved -- the network is not yet solvable, and why. The GUI reads this
    for inherited placeholders and to flag required values before Solve.
    """
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(_resolve_preview_sync, schema), timeout=30.0
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504, detail="resolve_preview exceeded its timeout."
        ) from None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=humanize(str(e), schema)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=humanize(str(e), schema)) from e


@app.post("/export")
async def export_results(schema: NetworkGraphSchema):
    soft_timeout = schema.solver_settings.timeout or 180.0
    hard_timeout = soft_timeout + 10.0
    try:
        # Solve quickly to get states
        # Six values, not five. `_solve_sync` gained a `diag` return in #168
        # (2026-05-27) and this call site was not updated, so every export that
        # got as far as solving died with "too many values to unpack
        # (expected 5)" -- for three and a half months, because `to_dataframe`
        # is unit-tested but the route that calls it was not. Named rather
        # than discarded with `*_` so the next addition breaks loudly here too.
        raw, node_results, element_results, edge_results, net, _diag = await asyncio.wait_for(
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
        raise HTTPException(status_code=400, detail=humanize(str(e), schema)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=humanize(str(e), schema)) from e


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
