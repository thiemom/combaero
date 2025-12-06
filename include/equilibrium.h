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

// -------------------------------------------------------------
// SMR + WGS equilibrium (Steam Methane Reforming + Water-Gas Shift)
// -------------------------------------------------------------
// Two coupled reactions:
//   SMR: CH4 + H2O <-> CO + 3H2  (endothermic, favored at high T)
//   WGS: CO + H2O <-> CO2 + H2   (exothermic, favored at low T)
//
// These functions solve for the equilibrium composition considering
// both reactions simultaneously.

// SMR+WGS equilibrium (isothermal) - equilibrate at input temperature
State smr_wgs_equilibrium(const State& in);

// SMR+WGS equilibrium (adiabatic) - solve for equilibrium T and composition
// This is the most physically realistic model for rich combustion products
State smr_wgs_equilibrium_adiabatic(const State& in);

#endif // EQUILIBRIUM_H