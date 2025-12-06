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

// -------------------------------------------------------------
// General Reforming + WGS equilibrium (handles ALL hydrocarbons)
// -------------------------------------------------------------
// Handles all CnHm hydrocarbons via steam reforming:
//   CnHm + n*H2O <-> n*CO + (n + m/2)*H2
// Plus WGS:
//   CO + H2O <-> CO2 + H2
//
// This is the most general equilibrium model, handling:
// CH4, C2H6, C3H8, iC4H10, nC5H12, nC6H14, nC7H16, etc.

// General reforming + WGS equilibrium (isothermal)
State reforming_equilibrium(const State& in);

// General reforming + WGS equilibrium (adiabatic)
State reforming_equilibrium_adiabatic(const State& in);

#endif // EQUILIBRIUM_H