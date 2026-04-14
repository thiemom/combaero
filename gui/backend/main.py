import asyncio
import io

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

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

app = FastAPI(title="CombAero Network GUI API")

# Configure CORS for local development (Vite typically runs on 5173)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    return {
        "message": "CombAero Network GUI API is running",
        "docs": "/docs",
        "health": "/health",
        "metadata": ["/metadata/species", "/metadata/presets"],
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


def _solve_sync(schema: NetworkGraphSchema):
    """
    Synchronous helper to build and solve the network.
    Runs in a separate thread to avoid blocking the FastAPI event loop.
    """
    net = build_network_from_schema(schema)
    solver = NetworkSolver(net)
    result = solver.solve(init_strategy=schema.solver_settings.init_strategy)
    success = bool(result.get("__success__", False))

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
            P_total=node_vals.pop("P_total", None),
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
        diag = elem_diags.get(elem_id, {})
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
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, element_results, edge_results, _ = await asyncio.to_thread(
            _solve_sync, schema
        )

        return SolveResponse(
            success=result.get("__success__", False),
            message=result.get("__message__", "Solve completed"),
            final_norm=result.get("__final_norm__"),
            node_results=node_results,
            element_results=element_results,
            edge_results=edge_results,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.post("/export")
async def export_results(schema: NetworkGraphSchema):
    try:
        # Solve quickly to get states
        _, node_results, element_results, _, _ = await asyncio.to_thread(_solve_sync, schema)

        # Convert to DataFrame
        import pandas as pd

        data = []
        for node_id, res in node_results.items():
            base_data = {
                "type": "node",
                "id": node_id,
                "T": res.state.T,
                "P": res.state.P,
                "P_total": res.state.P_total,
                "m_dot": getattr(res.state, "m_dot", None),
                "mach": res.state.mach,
            }
            if hasattr(res.state, "Y") and res.state.Y:
                for i, y_val in enumerate(res.state.Y):
                    base_data[f"Y[{i}]"] = y_val
            data.append(base_data)

        for elem_id, res in element_results.items():
            base_data = {
                "type": "element",
                "id": elem_id,
                "m_dot": res.m_dot,
            }
            if hasattr(res, "mach"):
                base_data["mach"] = res.mach
            if hasattr(res, "p_ratio"):
                base_data["p_ratio"] = res.p_ratio
            data.append(base_data)

        df = pd.DataFrame(data)

        # Stream as CSV
        stream = io.StringIO()
        df.to_csv(stream, index=False)

        return StreamingResponse(
            iter([stream.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=combaero_results.csv"},
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.get("/health")
async def health():
    return {"status": "ok"}
