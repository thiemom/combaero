import math
import time

import numpy as np

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PipeElement,
    PressureBoundary,
)

net = FlowNetwork()

# Inlet
X_air = cb.humid_air_composition(288.16, 101325, 0.6)
Y_air = cb.mole_to_mass(X_air)
m_dot_air = 1.4 * 10
T_air = 560.0
inlet = MassFlowBoundary("inlet", m_dot=m_dot_air, T_total=T_air, Y=Y_air)

# Outlet
outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=Y_air)

# Orifices
Cd = 0.8
d = 21.0  # mm
A = 0.25 * np.pi * (d / 1000) ** 2
orifice = OrificeElement("orifice", "outlet", "outlet", Cd=Cd, area=A)  # temp connection, will fix

# Pipe parameters
L_pipe = 0.1
D_pipe = 50e-3
A_pipe = 0.25 * np.pi * D_pipe**2

# Network topology
n_series = 10
n_parallel = 10

# Build pipe and chamber grids with string references
pipe_grid = {}
chamber_grid = {}
orifice_grid = {}
for i in range(n_series):
    for j in range(n_parallel):
        # For series-parallel grid: each chamber connects to the next one in same branch
        from_node = "inlet" if i == 0 else f"momentum_chamber_{i - 1}_{j}"
        to_node = f"momentum_chamber_{i}_{j}"
        pipe_grid[(i, j)] = PipeElement(
            f"pipe_{i}_{j}", from_node, to_node, length=L_pipe, diameter=D_pipe, roughness=2e-5
        )
        chamber_grid[(i, j)] = MomentumChamberNode(f"momentum_chamber_{i}_{j}", area=A_pipe)

# Create orifices - one for each parallel branch
for j in range(n_parallel):
    last_chamber_in_branch = f"momentum_chamber_{n_series - 1}_{j}"
    orifice_grid[(j)] = OrificeElement(
        f"orifice_{j}", last_chamber_in_branch, "outlet", Cd=Cd, area=A
    )

# Now create all the nodes referenced in elements
net.add_node(inlet)
net.add_node(outlet)

# Create momentum chambers
for chamber in chamber_grid.values():
    net.add_node(chamber)

# Add all elements
for pipe in pipe_grid.values():
    net.add_element(pipe)
for orifice in orifice_grid.values():
    net.add_element(orifice)

print(f"Network built: {len(net.nodes)} nodes, {len(net.elements)} elements")
print(f"Node IDs: {list(net.nodes.keys())}")
print(f"Element IDs: {list(net.elements.keys())}")

# Solver options for tweaking SciPy parameters
# Valid options vary by method:
# 'hybr': maxfev, xtol, factor
# 'lm': maxfev, xtol, ftol
# 'broyden1': maxfev, xtol, ftol, gtol, alpha, reduce_method, max_step
solver_options = {
    "maxfev": 5000,  # More evaluations for large network
    "xtol": 1e-4,  # Relaxed tolerance for faster convergence
    "ftol": 1e-8,  # Function tolerance (lm method)
}
solver_method = "lm"
solver_timeout = 120

# Test the network
try:
    solver = NetworkSolver(net)

    # Start wall time measurement
    start_time = time.time()
    result = solver.solve(method=solver_method, timeout=solver_timeout, options=solver_options)
    end_time = time.time()

    wall_time = end_time - start_time
    print(f"Wall time: {wall_time:.2f} seconds")
    print(f"Solver success: {result.get('__success__', False)}")
    print(f"Message: {result.get('__message__', '')}")
    print(f"Final norm: {result.get('__final_norm__', 'N/A')}")

    # Always print diagnostics
    print("\n=== Diagnostics ===")

    # Inlet pressure
    inlet_P = result.get("inlet.P_total", float("nan"))
    inlet_P_static = result.get("inlet.P", float("nan"))
    print(f"Inlet P_total: {inlet_P / 1e5:.3f} bar")
    print(f"Inlet P_static: {inlet_P_static / 1e5:.3f} bar")

    # Mach numbers in pipes
    mach_numbers = []
    print(f"Debug: Checking {len(pipe_grid)} pipes...")

    # Check a few sample results
    sample_keys = list(pipe_grid.keys())[:3]
    for i, j in sample_keys:
        pipe_element = pipe_grid[(i, j)]
        m_dot = result.get(f"{pipe_element.id}.m_dot", "MISSING")
        print(f"Debug: {pipe_element.id}.m_dot = {m_dot}")

    for (i, j), pipe_element in pipe_grid.items():
        # Get mass flow and properties
        m_dot = result.get(f"{pipe_element.id}.m_dot", float("nan"))
        if not math.isnan(m_dot):
            # Get upstream conditions
            from_node_id = "inlet" if i == 0 else f"momentum_chamber_{i - 1}_{j}"

            P = result.get(f"{from_node_id}.P", float("nan"))
            T = result.get(f"{from_node_id}.T", float("nan"))
            area = math.pi * (D_pipe / 2) ** 2

            if not math.isnan(P) and not math.isnan(T) and not math.isnan(area):
                # Calculate properties using combaero core functions
                rho = cb.density(T, P, Y_air)
                v = m_dot / (rho * area)
                mach = cb.mach_number(v, T, X_air)
                mach_numbers.append(mach)

    if mach_numbers:
        print(
            f"Mach in pipes - Min: {min(mach_numbers):.4f}, Mean: {sum(mach_numbers) / len(mach_numbers):.4f}, Max: {max(mach_numbers):.4f}"
        )
    else:
        print("Could not calculate Mach numbers - no valid data found")

except Exception as e:
    print(f"Error: {e}")
    import traceback

    traceback.print_exc()
