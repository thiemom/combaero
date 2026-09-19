#include "cooling_correlations.h"
#include "heat_transfer.h"
#include "math_constants.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace combaero::cooling {

double thermal_performance_factor(double Nu_ratio, double f_ratio) {
    // Webb & Eckert (1972) thermal-hydraulic performance factor.
    // eta = (Nu/Nu0) / (f/f0)^(1/3)
    // Compares heat transfer gain against pumping power penalty at equal velocity.
    // eta > 1: net improvement over smooth; eta < 1: pressure-drop generator.
    //
    // Cross-check against Singh & Ekkad tabulated data (e/D=0.0625, P/e=10, 90 deg):
    //   Re=30k:  2.45 / 6.50^(1/3) = 2.45 / 1.866 = 1.31  (table: 1.31)
    //   Re=100k: 1.95 / 6.30^(1/3) = 1.95 / 1.849 = 1.05  (table: 1.05)
    //   Re=400k: 1.55 / 6.20^(1/3) = 1.55 / 1.838 = 0.84  (table: 0.84)
    //
    // Reference: Webb, R.L. & Eckert, E.R.G. (1972)
    //   Int. J. Heat Mass Transfer 15(8), 1647-1658
    return Nu_ratio / std::pow(f_ratio, 1.0 / 3.0);
}

double adiabatic_wall_temperature(double T_hot, double T_coolant, double eta) {
    if (eta < 0.0 || eta > 1.0) {
        throw std::runtime_error(
            "Cooling effectiveness eta must be in range [0, 1], got " + std::to_string(eta)
        );
    }
    return T_hot - eta * (T_hot - T_coolant);
}

double cooled_wall_heat_flux(double T_hot, double T_coolant,
                             double h_hot, double h_coolant,
                             double eta,
                             double t_wall, double k_wall) {
    if (h_hot <= 0.0 || h_coolant <= 0.0) {
        throw std::runtime_error("Heat transfer coefficients must be positive");
    }
    if (t_wall < 0.0 || k_wall <= 0.0) {
        throw std::runtime_error("Wall thickness must be non-negative and conductivity positive");
    }

    const double T_aw = adiabatic_wall_temperature(T_hot, T_coolant, eta);
    const double U = overall_htc_wall(h_hot, h_coolant, t_wall, k_wall);
    return U * (T_aw - T_coolant);
}

}  // namespace combaero::cooling
