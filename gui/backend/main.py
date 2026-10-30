import asyncio
import io
import itertools
import json
import math
import threading
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles

import combaero as cb

from . import second_law_view, validity_view
from .graph_builder import build_network_from_schema
from .messages import humanize
from .resolve_preview import resolve_preview
from .runner import NetworkResult, _build_result_objects, _schema_maps
from .schemas import (
    NetworkGraphSchema,
    SolveResponse,
)
from .solver_setup import make_solver

_FRONTEND_DIST = Path(__file__).parent.parent / "frontend" / "dist"
# Worked examples from papers (#793): package data of combaero-gui
# (gui/pyproject.toml), so they load from an installed wheel too.
_EXAMPLES_DIR = Path(__file__).parent.parent / "examples"

app = FastAPI(title="CombAero Network GUI API")

# Configure CORS for local development (Vite typically runs on 5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["X-Combaero-Success", "X-Combaero-Outcome"],
)


# A spreadsheet runs a cell that starts with one of these as a formula. Labels
# are user text and travel with a shared network file, so a label such as
# =HYPERLINK("http://...") would run when a colleague opens the export. A
# leading apostrophe makes the cell text (OWASP "CSV injection").
_FORMULA_PREFIXES = ("=", "+", "-", "@", "\t", "\r")


def _neutralise_formulas(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        # Not `dtype == object`: pandas 3 gives text its own string dtype.
        if pd.api.types.is_numeric_dtype(out[col]) or pd.api.types.is_bool_dtype(out[col]):
            continue
        out[col] = out[col].map(
            lambda v: "'" + v if isinstance(v, str) and v.startswith(_FORMULA_PREFIXES) else v
        )
    return out


def _dist_file(name: str) -> FileResponse:
    # Without a built frontend (a repo checkout, the Vite dev server) the file
    # is missing, and a FileResponse on it is a bare 500 from the server.
    path = _FRONTEND_DIST / name
    if not path.is_file():
        raise HTTPException(status_code=404, detail=f"/{name} not found (frontend not built)")
    return FileResponse(path)


@app.get("/favicon.svg")
async def favicon_svg():
    return _dist_file("favicon.svg")


@app.get("/icons.svg")
async def icons_svg():
    return _dist_file("icons.svg")


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

#: Seconds /solve and /export wait beyond the solver's own (soft) timeout
#: before answering 504. The worker thread cannot be cancelled; this only
#: frees the request.
HARD_TIMEOUT_MARGIN = 10.0

# A worker thread outlives a request that answered 504 (#762): it cannot be
# cancelled, and it may still converge. Each request holds a token while it
# waits; the worker writes `_continuation_state` only while its token is live,
# under the same lock the request takes to give the token up. A result that
# arrives after its request timed out is discarded, so the next continuation
# solve never warm-starts from a state the user was told had failed.
_state_lock = threading.Lock()
_live_tokens: set[int] = set()
_token_counter = itertools.count(1)


def _open_token() -> int:
    token = next(_token_counter)
    with _state_lock:
        _live_tokens.add(token)
    return token


def _close_token(token: int) -> None:
    with _state_lock:
        _live_tokens.discard(token)


def _commit_continuation(state: dict, token: int | None) -> bool:
    """Store `state` for continuation unless `token`'s request has gone.

    `token=None` (a direct, synchronous caller) always commits.
    """
    global _continuation_state
    with _state_lock:
        if token is not None and token not in _live_tokens:
            return False
        _continuation_state = state
        return True


def _hard_timeout(timeout: float | None) -> float | None:
    """The request-level deadline for a solver `timeout`.

    None is no limit for the solver, so none for the request either (#762):
    the solve runs until it finishes. Any number, 0 included, is what the
    solver itself gets (0 times out at its first residual) plus the margin;
    a negative one counts as 0.
    """
    if timeout is None:
        return None
    return max(timeout, 0.0) + HARD_TIMEOUT_MARGIN


def _timeout_detail(what: str, timeout: float | None) -> str:
    # Reached with a finite timeout only: a null one sets no deadline (#762),
    # so "raise it, or clear it for no limit" is always something the user can do.
    if timeout is None:
        return f"{what} raised a timeout with no Max Solver Time set."
    hard_timeout = _hard_timeout(timeout)
    return (
        f"{what} did not return within {hard_timeout:g} s (Max Solver Time "
        f"{max(timeout, 0.0):g} s plus a {HARD_TIMEOUT_MARGIN:g} s grace): one step could "
        "not be interrupted. Raise Max Solver Time or clear it for no limit, simplify "
        "the network, or try a different init strategy."
    )


def _solve_sync(schema: NetworkGraphSchema, token: int | None = None):
    """
    Synchronous helper to build and solve the network.
    Runs in a separate thread to avoid blocking the FastAPI event loop.

    `token` is the request's (see `_open_token`); a converged state is kept
    for continuation only while it is live.
    """
    state = _continuation_state
    net = build_network_from_schema(schema)
    solver = make_solver(net, schema)
    x0 = None
    if schema.solver_settings.init_strategy == "continuation" and state:
        # Check topology fingerprint
        solver._build_x0()
        if state["unknown_names"] != list(solver.unknown_names):
            raise ValueError(
                "Topology mismatch: the current network unknowns do not match "
                "the saved continuation state. Please run a 'Default' solve first."
            )

        # Safety Reset: Check if initial guesses have changed
        current_guesses = {n.id: n.data.get("initial_guess", {}) for n in schema.nodes}
        if state.get("initial_guesses") != current_guesses:
            raise ValueError(
                "Initial guess mismatch: user overrides have changed since the last solve. "
                "Please run a 'Default' solve to incorporate these changes."
            )

        x0 = state["x"]

    result = solver.solve(
        method=schema.solver_settings.method,
        init_strategy=schema.solver_settings.init_strategy,
        timeout=schema.solver_settings.timeout,
        x0=x0,
        stall_handover=schema.solver_settings.stall_handover,
    )
    success = bool(result.get("__success__", False))

    if success:
        _commit_continuation(
            {
                "x": result["__x_solution__"],
                "unknown_names": result["__unknown_names__"],
                "initial_guesses": {n.id: n.data.get("initial_guess", {}) for n in schema.nodes},
            },
            token,
        )

    node_results, element_results, edge_results = _build_result_objects(result, net)
    diag = getattr(solver, "_diagnostic_data", {})
    return result, node_results, element_results, edge_results, net, diag


async def _run_solve(schema: NetworkGraphSchema):
    """`_solve_sync` in a worker thread, bounded by `_hard_timeout`.

    The token is given up however the wait ends, so a worker still running
    after a 504 (or an error) cannot write the continuation state.
    """
    token = _open_token()
    try:
        return await asyncio.wait_for(
            asyncio.to_thread(_solve_sync, schema, token),
            timeout=_hard_timeout(schema.solver_settings.timeout),
        )
    finally:
        _close_token(token)


@app.post("/solve", response_model=SolveResponse)
async def solve(schema: NetworkGraphSchema):
    # Hard async-level deadline: even if a single C++ Fanno/RK4 step blocks,
    # the event loop stays alive for health checks and subsequent requests.
    # HARD_TIMEOUT_MARGIN on top of the solver's own soft timeout; none when
    # the solver has none (#762).
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, element_results, edge_results, net, diag = await _run_solve(schema)

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
            detail=_timeout_detail("The solver", schema.solver_settings.timeout),
        ) from None
    except ValueError as e:
        raise HTTPException(status_code=400, detail=humanize(str(e), schema)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=humanize(str(e), schema)) from e


def _json_number(v: object) -> object:
    # The route has no response_model, so a NaN / inf reaches the JSON encoder
    # and fails AFTER the route returned: a bare 500. /solve's response model
    # writes non-finite floats as null; do the same here.
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


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
            "fields": {
                k: {"value": _json_number(fr.value), "source": fr.source}
                for k, fr in res.fields.items()
            },
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
    try:
        # Solve quickly to get states
        # Six values, not five. `_solve_sync` gained a `diag` return in #168
        # (2026-05-27) and this call site was not updated, so every export that
        # got as far as solving died with "too many values to unpack
        # (expected 5)" -- for three and a half months, because `to_dataframe`
        # is unit-tested but the route that calls it was not. Named rather
        # than discarded with `*_` so the next addition breaks loudly here too.
        raw, node_results, element_results, edge_results, net, _diag = await _run_solve(schema)

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
        _neutralise_formulas(df).to_csv(stream, index=False)

        return StreamingResponse(
            iter([stream.getvalue()]),
            media_type="text/csv",
            headers={
                "Content-Disposition": "attachment; filename=combaero_results.csv",
                # The solve this CSV came from, for the GUI to warn on without
                # parsing the file; every row carries the same columns.
                "X-Combaero-Success": "true" if result_obj.success else "false",
                "X-Combaero-Outcome": result_obj.outcome,
            },
        )
    except TimeoutError:
        raise HTTPException(
            status_code=504,
            detail=_timeout_detail("The export solve", schema.solver_settings.timeout),
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


def _example_files() -> dict[str, Path]:
    """The worked examples shipped with the package (#793), by id."""
    if not _EXAMPLES_DIR.is_dir():
        return {}
    return {p.stem: p for p in sorted(_EXAMPLES_DIR.glob("*.json"))}


@app.get("/examples")
async def list_examples():
    """The Examples menu: id, title and validation label of each shipped example."""
    out = []
    for eid, path in _example_files().items():
        info = json.loads(path.read_text(encoding="utf-8")).get("example", {})
        out.append({"id": eid, "title": info.get("title", eid), "label": info.get("label")})
    return {"examples": out}


@app.get("/examples/{example_id}")
async def get_example(example_id: str):
    """One example as a network file, with its description under ``example``.
    Only ids from the listing are served, never a path."""
    path = _example_files().get(example_id)
    if path is None:
        raise HTTPException(status_code=404, detail=f"No example {example_id!r}")
    return json.loads(path.read_text(encoding="utf-8"))


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
