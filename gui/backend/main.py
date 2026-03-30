import asyncio
import io

import combaero as cb
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse

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
    result = solver.solve()
    success = bool(result.get("__success__", False))

    node_results = {}
    for node_id in net.nodes:
        # Extract any extra node variables for diagnostics (e.g. transport properties)
        node_keys = {
            k.split(".", 1)[1]: v
            for k, v in result.items()
            if k.startswith(f"{node_id}.") and k.split(".", 1)[1] not in ["T", "P", "T_total", "P_total", "mach"] and not k.split(".", 1)[1].startswith("Y[")
        }

        # P and T extraction
        t_val = float(result.get(f"{node_id}.T", result.get(f"{node_id}.T_total", 300.0)))
        p_val = float(result.get(f"{node_id}.P", result.get(f"{node_id}.P_total", 101325.0)))

        # Check if actual total pressure is explicitly defined
        p_total_val = result.get(f"{node_id}.P_total")
        if p_total_val is not None:
            p_total_val = float(p_total_val)
        elif "P_total" in node_keys:
            p_total_val = float(node_keys.pop("P_total"))

        # Composition extraction
        y_vals = []
        y_idx = 0
        while f"{node_id}.Y[{y_idx}]" in result:
            y_vals.append(float(result[f"{node_id}.Y[{y_idx}]"]))
            y_idx += 1

        if not y_vals:
            # Fallback to dry air if not found (shouldn't happen in solved network)
            y_vals = list(cb.species.dry_air_mass())

        # Thermodynamics extraction using combaero utilities
        x_vals = [float(x) for x in cb.mass_to_mole(y_vals)]
        rho = float(cb.density(t_val, p_val, x_vals))
        h = float(cb.h_mass(t_val, x_vals))
        s = float(cb.s_mass(t_val, x_vals, p_val))

        # Mach extraction if available (common for momentum chambers)
        mach = result.get(f"{node_id}.mach")
        if mach is not None:
            mach = float(mach)

        state = StateResult(
            T=t_val,
            P=p_val,
            P_total=p_total_val,
            rho=rho,
            h=h,
            s=s,
            mach=mach,
            Y=y_vals,
            X=x_vals,
            **node_keys,
        )
        node_results[node_id] = NodeResult(state=state, success=success)

    element_results = {}
    for elem_id in net.elements:
        elem_keys = {k.split(".", 1)[1]: v for k, v in result.items() if k.startswith(f"{elem_id}.") and k != f"{elem_id}.m_dot"}
        m_dot = float(result.get(f"{elem_id}.m_dot", 0.0))
        element_results[elem_id] = ElementResult(m_dot=m_dot, success=success, **elem_keys)

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
        _, node_results, element_results, _ = await asyncio.to_thread(_solve_sync, schema)

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
                "mach": res.state.mach
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
