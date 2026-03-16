#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "state.h"
#include <vector>

namespace combaero {

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T, double well_depth);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);  // α = k/(ρ·cp) [m²/s]
double reynolds(double T, double P, const std::vector<double>& X, double V, double L);  // Re = ρVL/μ [-]
double peclet(double T, double P, const std::vector<double>& X, double V, double L);    // Pe = VL/α [-]

// -------------------------------------------------------------
// Transport State Bundle
// -------------------------------------------------------------

// Compute all transport properties at once
TransportState transport_state(double T, double P, const std::vector<double>& X);

} // namespace combaero

#endif // TRANSPORT_H
