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
    import combaero as cb

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
    for node_id, node in net.nodes.items():
        # P and T extraction from raw solver result
        t_val = float(result.get(f"{node_id}.T", result.get(f"{node_id}.T_total", 300.0)))
        p_val = float(result.get(f"{node_id}.P", result.get(f"{node_id}.P_total", 101325.0)))
        p_total_val = result.get(f"{node_id}.P_total")
        if p_total_val is not None:
            p_total_val = float(p_total_val)

        # Composition extraction
        y_vals = []
        y_idx = 0
        while f"{node_id}.Y[{y_idx}]" in result:
            y_vals.append(float(result[f"{node_id}.Y[{y_idx}]"]))
            y_idx += 1
        if not y_vals:
            y_vals = list(cb.species.dry_air_mass())

        # Construct a MixtureState for the diagnostics call
        # This matches the state the node saw during the last residual evaluation
        m_dot = float(result.get(f"{node_id}.m_dot", 0.0))  # Some nodes might have m_dot
        state_obj = cb.MixtureState(P=p_val, P_total=p_total_val or p_val, T=t_val, T_total=t_val, Y=y_vals, m_dot=m_dot)

        # New Single-Pass: Call diagnostics once
        # This returns all thermo/transport properties (h, s, rho, mu, k, Re, phi, etc.)
        diag = node.diagnostics(state_obj)

        x_vals = [float(x) for x in cb.mass_to_mole(y_vals)]

        state_res = StateResult(
            T=t_val,
            P=p_val,
            P_total=p_total_val,
            Y=y_vals,
            X=x_vals,
            **diag,
        )
        node_results[node_id] = NodeResult(state=state_res, success=success)

    element_results = {}
    for elem_id, elem in net.elements.items():
        m_dot = float(result.get(f"{elem_id}.m_dot", 0.0))

        # Element diagnostics (Inlet and Outlet states)
        # For now we use the states of the connected nodes
        state_in = node_results[elem.from_node].state
        state_out = node_results[elem.to_node].state

        mix_in = cb.MixtureState(P=state_in.P, P_total=state_in.P_total or state_in.P, T=state_in.T, Y=state_in.Y, m_dot=m_dot)
        mix_out = cb.MixtureState(P=state_out.P, P_total=state_out.P_total or state_out.P, T=state_out.T, Y=state_out.Y, m_dot=m_dot)

        diag = elem.diagnostics(mix_in, mix_out)

        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success, **diag)

    edge_results = {}
    for edge_id in net.walls:
        edge_keys = {
            k.split(".", 1)[1]: float(v) for k, v in result.items() if k.startswith(f"{edge_id}.")
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
