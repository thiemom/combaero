#ifndef THERMO_H
#define THERMO_H

#include "thermo_transport_data.h"
#include <iterator>
#include <vector>
#include <string>

// Universal gas constant [J/(mol·K)]
constexpr double R_GAS = 8.31446261815324;

// Constants needed for collision integrals
constexpr double BOLTZMANN = 1.380649e-23;  // J/K, Boltzmann constant
constexpr double AVOGADRO = 6.02214076e23;  // 1/mol, Avogadro's number
constexpr double PI = 3.14159265358979323846;

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass);

// Species metadata and lookup functions
std::string species_name(std::size_t species_index);
std::size_t species_index_from_name(const std::string& name);

// Number of thermo species in the internal data tables
std::size_t num_species();

// Molar mass lookup [kg/kmol equivalent: here stored as g/mol]
double species_molar_mass(std::size_t species_index);
double species_molar_mass_from_name(const std::string& name);

// Mixture properties
double mwmix(const std::vector<double>& X);

// NASA polynomial evaluations
// cp_R      : dimensionless heat capacity Cp/R
// cv_R      : dimensionless heat capacity Cv/R
// h_RT      : dimensionless enthalpy  H/(R*T)
// s_R       : dimensionless entropy   S/R at reference pressure and pure species
// g_over_RT : dimensionless Gibbs free energy G/(R*T) = H/(R*T) - S/R
double cp_R(std::size_t species_idx, double T);
double cv_R(std::size_t species_idx, double T);
double h_RT(std::size_t species_idx, double T);
double s_R(std::size_t species_idx, double T);
double g_over_RT(std::size_t species_idx, double T);

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double cv(double T, const std::vector<double>& X);
double density(double T, double P, const std::vector<double>& X);
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Combustion calculations
double oxygen_required_per_mol_fuel(std::size_t fuel_index);
double oxygen_required_per_kg_fuel(std::size_t fuel_index);
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);

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

    // Bilger mixture fraction (mass basis) already declared:
// double bilger_beta(const std::vector<double>& Y);
// double bilger_mixture_fraction(
//     const std::vector<double>& Y,
//     const std::vector<double>& Y_F,
//     const std::vector<double>& Y_O);

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
    
// Derivatives of thermodynamic properties with respect to temperature
// dh_dT        : dH/dT at constant composition [J/(mol·K)]
// ds_dT        : dS/dT of mixture standard-state entropy [J/(mol·K²)]
// dcp_dT       : dCp/dT [J/(mol·K²)]
// dg_over_RT_dT: d(G/(R*T))/dT at constant composition and pressure [1/K]
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
/**
 * @brief Derivative of dimensionless Gibbs free energy with respect to temperature.
 * 
 * This function evaluates the derivative of the dimensionless Gibbs free energy with respect to temperature at constant composition and pressure.
 * 
 * @param T Temperature in Kelvin.
 * @param X Mole fractions of the species.
 * @return Derivative of dimensionless Gibbs free energy with respect to temperature [1/K].
 */
double dg_over_RT_dT(double T, const std::vector<double>& X);

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T, double well_depth);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);

// Inverse calculations (solving for temperature)
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// Utility functions
void print_mixture_properties(double T, double P, const std::vector<double>& X);

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions);

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions);

// Convert mole fractions X_k to mass fractions Y_k using molar_masses.
std::vector<double> mole_to_mass(const std::vector<double>& X);

// Convert mass fractions Y_k to mole fractions X_k using molar_masses.
std::vector<double> mass_to_mole(const std::vector<double>& Y);

// Combustion calculations (duplicate declarations for convenience)
double oxygen_required_per_mol_fuel(std::size_t fuel_index);
double oxygen_required_per_kg_fuel(std::size_t fuel_index);
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);

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

#endif // THERMO_H
