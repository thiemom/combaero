#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "state.h"
#include <vector>

namespace combaero {

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T_star, double delta_star);

// Omega*(2,2) and dOmega*/dT*: quintic Hermite in (ln T*, ln Omega*) with
// tabulated node values, slopes and curvatures (C2); see src/transport.cpp.
struct Omega22Result {
    double omega = 0.0;
    double domega_dTstar = 0.0;
};
Omega22Result omega22_and_derivative(double T_star, double delta_star);

// Half-width, in units of R, of the C2 smoothing of the Eucken vibrational
// heat capacity cv_vib = max(0, cp/R - 5/2 - cv_rot/R) (#611). Outside
// |x| < w it is max(0, x) exactly; inside, its value differs by at most
// 0.156 w (at x = 0).
constexpr double CV_VIB_SMOOTH_HALF_WIDTH = 0.002;

// Mixture viscosity and thermal conductivity with analytic d/dT.
struct TransportDerivatives {
    double mu = 0.0;      // [Pa·s]
    double dmu_dT = 0.0;  // [Pa·s/K]
    double k = 0.0;       // [W/(m·K)]
    double dk_dT = 0.0;   // [W/(m·K^2)]
};
TransportDerivatives transport_and_dT(double T, double P, const std::vector<double>& X);

// d mu_mix / dX_k [Pa*s] at fixed T, analytic (Wilke's rule), for each
// species, the mole fractions taken unnormalised: chain it with dX/dY along a
// direction that keeps sum X = 1. For an absent species it is the one-sided
// derivative as that species enters.
std::vector<double> viscosity_dX(double T, double P, const std::vector<double>& X);

// Mixture viscosity and thermal conductivity with their gradients in the
// mole fractions at fixed T (#717): Wilke's rule and the mean of the series
// and parallel conductivity bounds. Same unnormalised-X convention as
// viscosity_dX; chain with dX/dY along a direction that keeps sum X = 1.
struct TransportDerivativesX {
    double mu = 0.0;                // [Pa*s]
    double k = 0.0;                 // [W/(m*K)]
    std::vector<double> dmu_dX;     // [Pa*s] per species
    std::vector<double> dk_dX;      // [W/(m*K)] per species
};
TransportDerivativesX transport_and_dX(double T, double P, const std::vector<double>& X);

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
