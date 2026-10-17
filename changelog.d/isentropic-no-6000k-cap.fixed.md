- The isentropic, energy and stagnation Newton solves (orifice and nozzle
  kernels, `T0_from_static`, `T_from_stagnation`) no longer cap T at 6000 K:
  an orifice fed by a 6872 K plenum read its density at 6000 K and passed 3.6x
  its flow. The stagnation conversions' 200 K floor is `thermo::T_FLOOR` (10 K),
  as the others', and so are the combustion closures' 134 K inlet floor and the
  quasi-1D nozzle's 150 K one.
