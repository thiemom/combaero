#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "state.h"
#include <vector>

namespace combaero {

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T_star, double delta_star);

// Omega*(2,2) and dOmega*/dT*: cubic Hermite in (ln T*, ln Omega*) with
// tabulated node values and slopes (C1); see src/transport.cpp.
struct Omega22Result {
    double omega = 0.0;
    double domega_dTstar = 0.0;
};
Omega22Result omega22_and_derivative(double T_star, double delta_star);

// Mixture viscosity and thermal conductivity with analytic d/dT.
struct TransportDerivatives {
    double mu = 0.0;      // [Pa·s]
    double dmu_dT = 0.0;  // [Pa·s/K]
    double k = 0.0;       // [W/(m·K)]
    double dk_dT = 0.0;   // [W/(m·K^2)]
};
TransportDerivatives transport_and_dT(double T, double P, const std::vector<double>& X);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);  // α = k/(ρ·cp) [m²/s]
double reynolds(double T, double P, const std::vector<double>& X, double V, double L);  // Re = rho*VL/mu [-]
double peclet(double T, double P, const std::vector<double>& X, double V, double L);    // Pe = VL/alpha [-]

// Reynolds number from pre-computed density and dynamic viscosity.
// Use this when a CompleteState/TransportState is already available to avoid
// re-evaluating mixture properties.
// Re = rho * v * L / mu  [-]
inline double reynolds_from_state(double rho, double v, double L, double mu) {
    return (rho * v * L) / mu;
}

// -------------------------------------------------------------
// Transport State Bundle
// -------------------------------------------------------------

// Compute all transport properties at once
TransportState transport_state(double T, double P, const std::vector<double>& X);

} // namespace combaero

#endif // TRANSPORT_H
