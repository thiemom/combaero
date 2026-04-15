feat: add head-loss correlation types to GUI pressure loss dropdown

- schemas.py: ConstantHeadLossData (zeta) and LinearThetaHeadLossData (k, zeta0)
- graph_builder.py: build_pressure_loss() now accepts area kwarg and handles
  constant_head / linear_theta_head, wired with data.area from CombustorData
- Inspector.tsx: two new dropdown options with conditional zeta / zeta0 inputs
  and formula hint strings (dP = zeta * q_in, dP = (k*Theta+zeta0) * q_in)
