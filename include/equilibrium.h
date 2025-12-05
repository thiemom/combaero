#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"

// -------------------------------------------------------------
// State-based WGS equilibrium functions
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
// CO + H2O <-> CO2 + H2
State wgs_equilibrium(const State& in);

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
State wgs_equilibrium_adiabatic(const State& in);

#endif // EQUILIBRIUM_H