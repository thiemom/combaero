#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "thermo.h"
#include <cstddef>
#include <vector>

struct WgsConfig {
    std::size_t i_CO;
    std::size_t i_H2O;
    std::size_t i_CO2;
    std::size_t i_H2;
};

// Solve for adiabatic T with WGS partial equilibrium
double solve_adiabatic_T_wgs(const std::vector<double>& n_in,
                             double H_target,
                             double T_guess,
                             const WgsConfig& cfg);

#endif // EQUILIBRIUM_H