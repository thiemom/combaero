import asyncio
import io

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

from combaero.network import NetworkSolver, results_to_dataframe

from .graph_builder import build_network_from_schema
from .schemas import NetworkGraphSchema, NodeResult, SolveResponse

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

    node_results = {}
    for node_id, node in net.nodes.items():
        state = node.get_state()
        node_results[node_id] = NodeResult(
            T=state.T, P=state.P, Y=state.Y, success=result.get("__success__", False)
        )

    return result, node_results, net


@app.post("/solve", response_model=SolveResponse)
async def solve(schema: NetworkGraphSchema):
    try:
        # Offload CPU-bound C++ solver to a worker thread
        result, node_results, _ = await asyncio.to_thread(_solve_sync, schema)

        return SolveResponse(
            success=result.get("__success__", False),
            message=result.get("__message__", "Solve completed"),
            node_results=node_results,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.post("/export")
async def export_results(schema: NetworkGraphSchema):
    try:
        # Solve quickly to get states
        _, _, net = await asyncio.to_thread(_solve_sync, schema)

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
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.get("/health")
async def health():
    return {"status": "ok"}
