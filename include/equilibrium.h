#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"

// -------------------------------------------------------------
// State-based WGS equilibrium functions
// -------------------------------------------------------------

// WGS equilibrium (isothermal) - equilibrate at input temperature
// CO + H2O <-> CO2 + H2
EquilibriumResult wgs_equilibrium(const State& in);

// WGS equilibrium (adiabatic) - solve for equilibrium T and composition
EquilibriumResult wgs_equilibrium_adiabatic(const State& in);

// -------------------------------------------------------------
// SMR + WGS equilibrium (Steam Methane Reforming + Water-Gas Shift)
// -------------------------------------------------------------
// DEPRECATED: use reforming_equilibrium / reforming_equilibrium_adiabatic instead.
// reforming_equilibrium handles CH4 (and all other CnHm) via the same
// sequential solver and produces identical results for CH4-only inputs.

// SMR+WGS equilibrium (isothermal) - equilibrate at input temperature
EquilibriumResult smr_wgs_equilibrium(const State& in);

// SMR+WGS equilibrium (adiabatic) - solve for equilibrium T and composition
EquilibriumResult smr_wgs_equilibrium_adiabatic(const State& in);

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
EquilibriumResult reforming_equilibrium(const State& in);

// General reforming + WGS equilibrium (adiabatic)
EquilibriumResult reforming_equilibrium_adiabatic(const State& in);

// -------------------------------------------------------------
// Combustion + Equilibrium (convenience functions)
// -------------------------------------------------------------
// These functions combine complete combustion with equilibrium in one call.
// They are useful when starting from an unburned fuel+air mixture.
//
// The workflow is:
//   1. Complete combustion (fuel + O2 -> CO2 + H2O, with excess fuel remaining)
//   2. Reforming + WGS equilibrium on the products
//
// This handles the common case where users want to go directly from
// an unburned mixture to equilibrium products.

// Combustion + equilibrium (adiabatic)
// Input: unburned fuel+air mixture at inlet temperature
// Output: equilibrium products at adiabatic flame temperature
EquilibriumResult combustion_equilibrium(const State& in);

#endif // EQUILIBRIUM_H
