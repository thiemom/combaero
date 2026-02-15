#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <vector>

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

// Bundle of transport properties for a gas mixture
// All properties computed from (T, P, X) in a single call
struct TransportState {
    double T;      // Temperature [K] (input, echoed back)
    double P;      // Pressure [Pa] (input, echoed back)
    double rho;    // Density [kg/m³]
    double mu;     // Dynamic viscosity [Pa·s]
    double k;      // Thermal conductivity [W/(m·K)]
    double nu;     // Kinematic viscosity [m²/s]
    double alpha;  // Thermal diffusivity [m²/s]
    double Pr;     // Prandtl number [-]
    double cp;     // Specific heat at constant pressure [J/(kg·K)]
};

// Compute all transport properties at once
// Parameters:
//   T : temperature [K]
//   P : pressure [Pa]
//   X : mole fractions [-]
// Returns: TransportState struct with all properties
TransportState transport_state(double T, double P, const std::vector<double>& X);

#endif // TRANSPORT_H
