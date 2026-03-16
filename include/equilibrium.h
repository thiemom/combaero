#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"
#include "combustion.h"

namespace combaero {

// -------------------------------------------------------------
// State-based WGS equilibrium functions
// -------------------------------------------------------------
EquilibriumResult wgs_equilibrium(const State &in);
EquilibriumResult wgs_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// SMR + WGS equilibrium (DEPRECATED: use reforming_equilibrium)
// -------------------------------------------------------------
EquilibriumResult smr_wgs_equilibrium(const State &in);
EquilibriumResult smr_wgs_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// Reforming + WGS equilibrium
// -------------------------------------------------------------
// Higher-fidelity product model for rich hydrocarbon combustion.
// Produces CO, H2, CO2, H2O, N2 (ignoring soot/solid carbon).
EquilibriumResult reforming_equilibrium(const State &in);
EquilibriumResult reforming_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// High-level combustion interface
// -------------------------------------------------------------

// Computes product state using the specified method (Complete or Equilibrium)
EquilibriumResult combustion_equilibrium(const State &reactants,
                                         bool smooth_phi0 = true,
                                         double k0 = SMOOTHING_K_PHI0);

} // namespace combaero

#endif // EQUILIBRIUM_H
