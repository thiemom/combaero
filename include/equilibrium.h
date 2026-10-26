#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "state.h"
#include <vector>
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
//
// The reactions (steam reforming of each hydrocarbon, then WGS) are
// equilibrated in turn and the round is REPEATED until no amount moves by more
// than REFORMING_PASS_RTOL of the total: the fixed point is the simultaneous
// equilibrium. A single round leaves reforming short of equilibrium once WGS
// has moved CO, H2 and H2O -- at 650-850 K, the pre-reformer range, H2 came
// out 7-49% low against Cantera.
constexpr int REFORMING_MAX_PASSES = 2000;
constexpr double REFORMING_PASS_RTOL = 1e-14;
EquilibriumResult reforming_equilibrium(const State &in);
EquilibriumResult reforming_equilibrium_adiabatic(const State &in);

// -------------------------------------------------------------
// High-level combustion interface
// -------------------------------------------------------------

// Computes product state using the specified method (Complete or Equilibrium)
EquilibriumResult combustion_equilibrium(const State &reactants,
                                         bool smooth_phi0 = true,
                                         double k0 = SMOOTHING_K_PHI0);

// combustion_equilibrium with analytic derivatives (#716): implicit in each
// stage's extent (reforming per hydrocarbon, WGS) and in the adiabatic root,
// chained through complete_combustion_products_and_jacobian. X is taken
// unnormalised (chain with a dX that keeps sum X = 1). An extent on a bound
// takes the bound's derivative.
struct EquilibriumJacobian {
    State state;                                // as combustion_equilibrium
    double dT_dT_in = 0.0;
    double dT_dP = 0.0;
    std::vector<double> dT_dX;                  // [reactant k]
    std::vector<std::vector<double>> dXout_dX;  // [product j][reactant k]
    std::vector<double> dXout_dT_in;
    std::vector<double> dXout_dP;
};
EquilibriumJacobian combustion_equilibrium_and_jacobian(const State &reactants,
                                                        bool smooth_phi0 = true,
                                                        double k0 = SMOOTHING_K_PHI0);

} // namespace combaero

#endif // EQUILIBRIUM_H
