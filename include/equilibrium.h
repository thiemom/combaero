#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"
#include "thermo.h"
#include <cstddef>
#include <vector>

struct WgsConfig {
    std::size_t i_CO;
    std::size_t i_H2O;
    std::size_t i_CO2;
    std::size_t i_H2;
};

// Solve for adiabatic T with WGS partial equilibrium (legacy API)
double solve_adiabatic_T_wgs(const std::vector<double>& n_in,
                             double H_target,
                             double T_guess,
                             const WgsConfig& cfg);

// -------------------------------------------------------------
// State-based WGS equilibrium functions
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
// CO + H2O <-> CO2 + H2
State wgs_equilibrium(const State& in);

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
State wgs_equilibrium_adiabatic(const State& in);

#endif // EQUILIBRIUM_H