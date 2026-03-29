import asyncio
import io

import combaero as cb
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

from combaero.network import NetworkSolver, results_to_dataframe

from .graph_builder import build_network_from_schema
from .schemas import ElementResult, NetworkGraphSchema, NodeResult, SolveResponse

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
    }


def _solve_sync(schema: NetworkGraphSchema):
    """
    Synchronous helper to build and solve the network.
    Runs in a separate thread to avoid blocking the FastAPI event loop.
    """
    net = build_network_from_schema(schema)
    solver = NetworkSolver(net)
    result = solver.solve()
    success = bool(result.get("__success__", False))

    node_results = {}
    default_y = list(cb.species.dry_air_mass())
    for node_id in net.nodes:
        t_val = float(result.get(f"{node_id}.T", result.get(f"{node_id}.T_total", 300.0)))
        p_val = float(result.get(f"{node_id}.P", result.get(f"{node_id}.P_total", 101325.0)))

        y_vals = []
        y_idx = 0
        while f"{node_id}.Y[{y_idx}]" in result:
            y_vals.append(float(result[f"{node_id}.Y[{y_idx}]"]))
            y_idx += 1

        if not y_vals:
            y_vals = default_y

        node_results[node_id] = NodeResult(
            T=t_val,
            P=p_val,
            Y=y_vals,
            success=success,
        )

    element_results = {}
    for elem_id in net.elements:
        m_dot = float(result.get(f"{elem_id}.m_dot", 0.0))
        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success)

    return result, node_results, element_results, net


@app.post("/solve", response_model=SolveResponse)
async def solve(schema: NetworkGraphSchema):
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, element_results, _ = await asyncio.to_thread(_solve_sync, schema)

        return SolveResponse(
            success=result.get("__success__", False),
            message=result.get("__message__", "Solve completed"),
            node_results=node_results,
            element_results=element_results,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.post("/export")
async def export_results(schema: NetworkGraphSchema):
    try:
        # Solve quickly to get states
        _, _, _, net = await asyncio.to_thread(_solve_sync, schema)

        # Convert to DataFrame
        df = results_to_dataframe(net)

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
