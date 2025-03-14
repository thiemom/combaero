#ifndef THERMO_H
#define THERMO_H

#include "thermo_transport_data.h"
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

// Species name lookup functions
std::string species_name(int species_index);
int species_index_from_name(const std::string& name);

// Mixture properties
double mwmix(const std::vector<double>& X);

// NASA polynomial evaluations
double cp_R(int species_idx, double T);
double h_RT(int species_idx, double T);
double s_R(int species_idx, double T);

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, double P, const std::vector<double>& X);
double cv(double T, const std::vector<double>& X);
double density(double T, double P, const std::vector<double>& X);
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Combustion calculations
double oxygen_required_per_mol_fuel(int fuel_index);
double oxygen_required_per_kg_fuel(int fuel_index);
double oxygen_required_per_mol_mixture(const std::vector<double>& X);
double oxygen_required_per_kg_mixture(const std::vector<double>& X);

// Derivatives
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T, double well_depth);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);

// Inverse calculations (solving for temperature)
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, int max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, int max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, int max_iter = 50);

// Utility functions
void print_mixture_properties(double T, double P, const std::vector<double>& X);

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions);

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions);

#endif // THERMO_H
