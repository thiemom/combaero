import timeit


def setup_code():
    return """
import combaero as cb
import numpy as np
X = [0.0] * 14
X[0] = 0.79  # N2
X[1] = 0.21  # O2
T = 300.0
P = 101325.0
mdot = 0.1
D = 0.05
L = 1.0
roughness = 1e-5
A = 3.14159 * D**2 / 4.0
rho = cb.density(T, P, X)
u = mdot / (rho * A)

# Pre-allocate state and output buffers
state = cb.State()
state.TPX = (T, P, X)
out_arr = np.zeros(2, dtype=np.float64)
"""


def print_result(name, code, setup, runs=100000):
    time_taken = timeit.timeit(stmt=code, setup=setup, number=runs)
    ns_per_call = (time_taken / runs) * 1e9
    us_per_call = ns_per_call / 1000.0
    print(f"{name:<40} | {us_per_call:>8.3f} μs | {ns_per_call:>8.0f} ns")


print("CombAero Function Profiling (per call)")
print("-" * 65)
print(f"{'Operation':<40} | {'Time (us)':>8} | {'Time (ns)':>8}")
print("-" * 65)

# 1. State Object & Polynomials (expected to be very fast)
print_result("State() instantiation", "s = cb.State()", setup_code())
print_result("State.TPX setter", "state.TPX = (T, P, X)", setup_code())
print_result("Polynomial: h(T)", "cb.h(T, X)", setup_code())
print_result("Polynomial: cp(T)", "cb.cp(T, X)", setup_code())

# 2. Equation of State (fast)
print_result("EOS: density(T, P, X)", "cb.density(T, P, X)", setup_code())
print_result("Acoustic: speed_of_sound(T, X)", "cb.speed_of_sound(T, X)", setup_code())

# 3. Transport Properties (expected to be slower due to collision integrals/mixing)
print_result("Transport: viscosity(T, P, X)", "cb.viscosity(T, P, X)", setup_code())
print_result(
    "Transport: thermal_conductivity(T, P, X)", "cb.thermal_conductivity(T, P, X)", setup_code()
)
print_result("Transport: transport_state(T, P, X)", "cb.transport_state(T, P, X)", setup_code())

# 4. Incompressible Elements (require thermo + correlations)
print_result(
    "Flow: orifice_flow_thermo (Cd=0.6)",
    "cb.orifice_flow_thermo(T, P, X, P-1000, A, 0.6)",
    setup_code(),
)
print_result(
    "Flow: pipe_flow_rough (Haaland)",
    "cb.pipe_flow_rough(T, P, X, u, L, D, roughness, 'haaland')",
    setup_code(),
)
print_result(
    "Flow: pipe_flow_rough (Colebrook)",
    "cb.pipe_flow_rough(T, P, X, u, L, D, roughness, 'colebrook')",
    setup_code(),
)

# 5. Compressible Elements (Iterative/ODE integration)
print_result(
    "Flow: fanno_pipe (constant f, 10 steps)",
    "cb.fanno_pipe(T, P, u, L, D, 0.02, X, 10)",
    setup_code(),
    runs=10000,
)
print_result(
    "Flow: fanno_pipe_rough (10 steps)",
    "cb.fanno_pipe_rough(T, P, u, L, D, roughness, X, 'haaland', 10)",
    setup_code(),
    runs=10000,
)
print_result(
    "Flow: fanno_pipe_rough (100 steps)",
    "cb.fanno_pipe_rough(T, P, u, L, D, roughness, X, 'haaland', 100)",
    setup_code(),
    runs=1000,
)

print("-" * 65)
print("6. PyBind11 Return Type Latency Focus")
print("-" * 65)
print_result("Return: Primitive (float)", "state.T", setup_code())
print_result("Return: Tuple (float, float, list)", "state.TPX", setup_code())
print_result("Return: List[float] (std::vector)", "cb.mole_to_mass(X)", setup_code())
print_result(
    "Return: Object (IncompressibleFlowSolution)",
    "cb.pipe_flow_rough(T, P, X, u, L, D, roughness, 'haaland')",
    setup_code(),
)
print_result(
    "Return: Tuple[float, float, float] (htc_pipe)",
    "cb.htc_pipe(T, P, X, u, D, 'gnielinski', True, 1.0, 0.0)",
    setup_code(),
)
print_result(
    "Return: Pass-by-ref (numpy array)", "cb._core.test_pass_by_ref(T, P, out_arr)", setup_code()
)
print("-" * 65)
