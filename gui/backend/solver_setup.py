"""One solver configuration for every backend path (the /solve and /export
routes and NetworkRunner), so a saved drawing solved from a script gets the
same retries and limits as in the GUI."""

from __future__ import annotations

from combaero.network import FlowNetwork, NetworkSolver

from .schemas import NetworkGraphSchema


def make_solver(net: FlowNetwork, schema: NetworkGraphSchema) -> NetworkSolver:
    solver = NetworkSolver(net)
    # The Newton-homotopy last retry (#642), bounded only by the max solve
    # time: whatever the cold solve and its retries leave of it.
    solver.homotopy_retry = True
    solver.homotopy_budget = None
    if schema.solver_settings.incompressible_mach_limit is not None:
        solver.incompressible_mach_limit = schema.solver_settings.incompressible_mach_limit
    return solver
