import numpy as np
import combaero as cba
from combaero.network import LosslessConnectionElement

# Create State
s1 = cba.State(300, 100000, np.array([0.21, 0.79]))
s1.mdot = 0.5

print("--- Lossless Connection ---")
pipe_zero = LosslessConnectionElement(L=5.0)
pipe_zero.in_nodes = [s1]
pipe_zero.out_nodes = [s1]
prof_zero = pipe_zero.get_spatial_profile(5)
for st in prof_zero:
    print(st)

print("\n--- Incompressible Pipe ---")
res = cba.pipe_flow_rough(300, 100000, np.array([0.21, 0.79]), 50.0, 10.0, 0.1, 1e-5, "haaland", 5, True)
for st in res.profile:
    print(st)
