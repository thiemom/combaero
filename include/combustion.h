#ifndef COMBUSTION_H
#define COMBUSTION_H

#include "state.h"
#include <cstddef>
#include <vector>
#include <string>

// Combustion calculations - oxygen requirements
double oxygen_required_per_mol_fuel(std::size_t fuel_index);
double oxygen_required_per_kg_fuel(std::size_t fuel_index);
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);

// Fuel lower heating value (LHV) from complete combustion to CO2 + H2O(g)
// X_fuel must be mole fractions over the global species set.
double fuel_lhv_molar(const std::vector<double>& X_fuel, const double reference_temperature = 298.15);  // [J/mol fuel mixture]
double fuel_lhv_mass(const std::vector<double>& X_fuel, const double reference_temperature = 298.15);   // [J/kg fuel mixture]

// Combustion calculations - dry air requirements (using standard dry air composition)
double dryair_required_per_mol_fuel(std::size_t fuel_index);
double dryair_required_per_kg_fuel(std::size_t fuel_index);
double dryair_required_per_mol_mixture(const std::vector<double>& X);
double dryair_required_per_kg_mixture(const std::vector<double>& X);

// Equivalence ratio (mole basis) for multi-species fuel + oxidizer.
// X_* are mole fractions over the same species set as species_names.

// Compute φ for a given unreacted mixture X_mix that is formed only by
// mixing a fuel stream (X_fuel) and an oxidizer stream (X_ox).
double equivalence_ratio_mole(
    const std::vector<double>& X_mix,
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox);

// Given target φ, and definitions of the fuel and oxidizer streams,
// construct the unreacted mixture mole fractions X_mix.
std::vector<double> set_equivalence_ratio_mole(
    double phi,
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox);

// Equivalence ratio (mass basis) for multi-species fuel + oxidizer.
// Y_* are mass fractions over the same species set as species_names.

// Compute φ for a given unreacted mixture Y_mix that is formed only by
// mixing a fuel stream (Y_fuel) and an oxidizer stream (Y_ox).
double equivalence_ratio_mass(
    const std::vector<double>& Y_mix,
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox);

// Given target φ, and definitions of the fuel and oxidizer streams,
// construct the unreacted mixture mass fractions Y_mix.
std::vector<double> set_equivalence_ratio_mass(
    double phi,
    const std::vector<double>& Y_fuel,
    const std::vector<double>& Y_ox);

// Stoichiometric Bilger mixture fraction Z_st for given fuel & oxidizer streams
// (mass fractions Y_F, Y_O).
double bilger_stoich_mixture_fraction_mass(
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O);

// Convert Bilger mixture fraction Z -> equivalence ratio φ (mass basis)
// for given fuel & oxidizer streams.
double equivalence_ratio_from_bilger_Z_mass(
    double Z,
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O);

// Convert equivalence ratio φ (mass basis) -> Bilger mixture fraction Z
// for given fuel & oxidizer streams.
double bilger_Z_from_equivalence_ratio_mass(
    double phi,
    const std::vector<double>& Y_F,
    const std::vector<double>& Y_O);

// Complete combustion to CO2 and H2O.
// - If O2 >= stoich: all fuel burns, possible O2 left over.
// - If 0 < O2 < stoich: all fuels burn with the same fraction f of their
//   stoichiometric amount based on available O2, and O2 is fully consumed.
// - If no fuel or no O2: mixture is returned unchanged.
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X);

// Overload that also returns the fuel-burn fraction f (0 <= f <= 1).
// - f = 1 for fuel-limited cases (O2 in excess or exactly stoichiometric).
// - 0 < f < 1 for O2-limited cases.
// - f = 0 if no combustion occurs (no fuel or no O2).
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X, double& fuel_burn_fraction);

// Mixture fraction (Bilger) utilities
// All Y*, mass fractions over the same species set as species_names / molar_masses.
double bilger_beta(const std::vector<double>& Y);
double bilger_mixture_fraction(
    const std::vector<double>& Y,    // local composition
    const std::vector<double>& Y_F,  // pure fuel-stream composition
    const std::vector<double>& Y_O   // pure oxidizer-stream composition
);

// Convenience overload: Bilger mixture fraction from mole fractions X.
// Internally converts X, X_F, X_O to mass fractions and calls
// bilger_mixture_fraction(...) above.
double bilger_mixture_fraction_from_moles(
    const std::vector<double>& X,
    const std::vector<double>& X_F,
    const std::vector<double>& X_O
);

// -------------------------------------------------------------
// Stream-based equivalence ratio helpers
// -------------------------------------------------------------

// Given a target equivalence ratio phi, a fuel stream (composition only, mdot ignored),
// and an oxidizer stream (with mdot set), return a copy of the fuel stream with mdot
// set to achieve the target phi when mixed with the oxidizer.
//
// Example:
//   Stream fuel, air;
//   fuel.set_T(300).set_X(X_CH4);       // mdot not needed
//   air.set_T(298).set_X(X_air).set_mdot(10.0);
//   Stream fuel_phi = set_fuel_stream_for_phi(0.8, fuel, air);
//   Stream mixed = mix({fuel_phi, air});
Stream set_fuel_stream_for_phi(double phi, const Stream& fuel, const Stream& oxidizer);

// -------------------------------------------------------------
// Inverse solvers for fuel/oxidizer streams (complete combustion only)
// -------------------------------------------------------------
// These functions find the fuel or oxidizer mass flow rate to achieve a target
// property in the burned products (complete combustion to CO2/H2O, no WGS).
//
// For Tad solvers: The same Tad can be achieved on both lean (O2 excess) and
// rich (fuel excess) sides of stoichiometric. Use the 'lean' parameter to
// select which side to search (default: lean=true).
//
// For O2/CO2 solvers: Only lean combustion is supported (O2 > 0 in products).

// --- Find fuel stream (oxidizer mdot fixed) ---

// Find fuel mdot to achieve target adiabatic flame temperature [K].
// T_ad_target must be in (T_oxidizer, T_ad_stoich] K.
// lean: if true (default), search on lean side (O2 excess); if false, search on rich side.
// phi_max: maximum equivalence ratio for rich side search (default: 10.0).
Stream set_fuel_stream_for_Tad(double T_ad_target, const Stream& fuel, const Stream& oxidizer,
                                double tol = 1.0, std::size_t max_iter = 100, bool lean = true,
                                double phi_max = 10.0);

// Find fuel mdot to achieve target O2 mole fraction in burned products (wet basis).
Stream set_fuel_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                               double tol = 1e-6, std::size_t max_iter = 100);

// Find fuel mdot to achieve target O2 mole fraction in burned products (dry basis).
Stream set_fuel_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol = 1e-6, std::size_t max_iter = 100);

// Find fuel mdot to achieve target CO2 mole fraction in burned products (wet basis).
Stream set_fuel_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                double tol = 1e-6, std::size_t max_iter = 100);

// Find fuel mdot to achieve target CO2 mole fraction in burned products (dry basis).
Stream set_fuel_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol = 1e-6, std::size_t max_iter = 100);

// --- Find oxidizer stream (fuel mdot fixed) ---

// Find oxidizer mdot to achieve target adiabatic flame temperature [K].
// T_ad_target must be in (T_fuel, T_ad_stoich] K.
// lean: if true (default), search on lean side (O2 excess); if false, search on rich side.
// phi_max: maximum equivalence ratio for search range (default: 10.0).
Stream set_oxidizer_stream_for_Tad(double T_ad_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol = 1.0, std::size_t max_iter = 100, bool lean = true,
                                    double phi_max = 10.0);

// Find oxidizer mdot to achieve target O2 mole fraction in burned products (wet basis).
Stream set_oxidizer_stream_for_O2(double X_O2_target, const Stream& fuel, const Stream& oxidizer,
                                   double tol = 1e-6, std::size_t max_iter = 100);

// Find oxidizer mdot to achieve target O2 mole fraction in burned products (dry basis).
Stream set_oxidizer_stream_for_O2_dry(double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                       double tol = 1e-6, std::size_t max_iter = 100);

// Find oxidizer mdot to achieve target CO2 mole fraction in burned products (wet basis).
Stream set_oxidizer_stream_for_CO2(double X_CO2_target, const Stream& fuel, const Stream& oxidizer,
                                    double tol = 1e-6, std::size_t max_iter = 100);

// Find oxidizer mdot to achieve target CO2 mole fraction in burned products (dry basis).
Stream set_oxidizer_stream_for_CO2_dry(double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer,
                                        double tol = 1e-6, std::size_t max_iter = 100);

// -------------------------------------------------------------
// State-based combustion functions
// -------------------------------------------------------------

// Complete combustion to CO2 and H2O (adiabatic)
// Input: unburned state with T, P, X
// Output: burned state with adiabatic flame temperature
State complete_combustion(const State& in);

// Complete combustion to CO2 and H2O (isothermal)
// Input: unburned state with T, P, X
// Output: burned state at same temperature
State complete_combustion_isothermal(const State& in);

// -------------------------------------------------------------
// Combustion State Dataclass
// -------------------------------------------------------------
// CombustionState is defined in state.h

// Compute combustion state from equivalence ratio (typical for calculations)
// Parameters:
//   X_fuel       : fuel composition (mole fractions) [-]
//   X_ox         : oxidizer composition (mole fractions) [-]
//   phi          : equivalence ratio [-] (INPUT)
//   T_reactants  : reactant temperature [K]
//   P            : pressure [Pa]
//   fuel_name    : optional fuel label (default: "")
// Returns: CombustionState with reactants, products, phi, mixture_fraction
CombustionState combustion_state(
    const std::vector<double>& X_fuel,
    const std::vector<double>& X_ox,
    double phi,
    double T_reactants,
    double P,
    const std::string& fuel_name = ""
);

// Compute combustion state from measured streams (typical for lab data)
// Parameters:
//   fuel_stream : fuel stream with mdot, T, X
//   ox_stream   : oxidizer stream with mdot, T, X
//   fuel_name   : optional fuel label (default: "")
// Returns: CombustionState with phi COMPUTED from mass flow rates
CombustionState combustion_state_from_streams(
    const Stream& fuel_stream,
    const Stream& ox_stream,
    const std::string& fuel_name = ""
);

#endif // COMBUSTION_H
