#include "../include/pipe_flow.h"
#include "../include/friction.h"
#include "../include/incompressible.h"
#include "../include/thermo.h"
#include "../include/transport.h"

#include <stdexcept>
#include <string>

std::tuple<double, double, double> pressure_drop_pipe(
    double T,
    double P,
    const std::vector<double>& X,
    double v,
    double D,
    double L,
    double roughness,
    const std::string& correlation)
{
    // Input validation
    if (T <= 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: T must be positive");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: P must be positive");
    }
    if (D <= 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: D must be positive");
    }
    if (L < 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: L must be non-negative");
    }
    if (v < 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: v must be non-negative");
    }
    if (roughness < 0.0) {
        throw std::invalid_argument("pressure_drop_pipe: roughness must be non-negative");
    }

    // Step 1: Calculate fluid properties
    double rho = density(T, P, X);      // [kg/m³]
    double mu = viscosity(T, P, X);     // [Pa·s]

    // Step 2: Calculate Reynolds number
    // Re = ρ * v * D / μ
    double Re = reynolds(T, P, X, v, D);

    // Step 3: Calculate friction factor using specified correlation
    double f;
    double e_D = roughness / D;  // Relative roughness [-]

    if (correlation == "haaland") {
        f = friction_haaland(Re, e_D);
    } else if (correlation == "serghides") {
        f = friction_serghides(Re, e_D);
    } else if (correlation == "colebrook") {
        f = friction_colebrook(Re, e_D);
    } else if (correlation == "petukhov") {
        f = friction_petukhov(Re);
    } else {
        throw std::invalid_argument(
            "pressure_drop_pipe: unknown correlation '" + correlation + "'. "
            "Valid options: 'haaland', 'serghides', 'colebrook', 'petukhov'");
    }

    // Step 4: Calculate pressure drop using Darcy-Weisbach equation
    // ΔP = f * (L/D) * (ρ * v² / 2)
    double dP = pipe_dP(v, L, D, f, rho);

    return std::make_tuple(dP, Re, f);
}
