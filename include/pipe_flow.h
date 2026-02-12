#ifndef PIPE_FLOW_H
#define PIPE_FLOW_H

#include <string>
#include <tuple>
#include <vector>

// -------------------------------------------------------------
// Composite Pipe Flow Functions
// -------------------------------------------------------------
// High-level functions that combine thermodynamic properties,
// transport properties, friction correlations, and flow equations.
//
// These simplify common workflows by doing multiple steps in one call.

// Pressure drop in pipe with automatic property evaluation.
//
// Combines:
//   1. Density and viscosity from thermodynamic state
//   2. Reynolds number calculation
//   3. Friction factor from correlation
//   4. Darcy-Weisbach pressure drop
//
// Parameters:
//   T           : temperature [K]
//   P           : pressure [Pa]
//   X           : mole fractions [mol/mol]
//   v           : flow velocity [m/s]
//   D           : pipe diameter [m]
//   L           : pipe length [m]
//   roughness   : absolute roughness [m] (default: 0.0 = smooth)
//   correlation : friction correlation name (default: "haaland")
//                 Options: "haaland", "serghides", "colebrook", "petukhov"
//
// Returns: tuple of (dP [Pa], Re [-], f [-])
//
// Example:
//   auto [dP, Re, f] = pressure_drop_pipe(300, 101325, X_air, 10.0, 0.1, 100.0);
//
std::tuple<double, double, double> pressure_drop_pipe(
    double T,
    double P,
    const std::vector<double>& X,
    double v,
    double D,
    double L,
    double roughness = 0.0,
    const std::string& correlation = "haaland");

#endif // PIPE_FLOW_H
