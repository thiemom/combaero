#ifndef THERMO_H
#define THERMO_H

#include "state.h"
#include <cstddef>
#include <string>
#include <vector>
#include "transport.h"

namespace combaero::thermo {

// Universal gas constant [J/(mol·K)]
constexpr double R_GAS = 8.31446261815324;

// Boltzmann constant [J/K]
constexpr double BOLTZMANN = 1.380649e-23;

// Avogadro's number [1/mol]
constexpr double AVOGADRO = 6.02214076e23;

} // namespace combaero::thermo

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass);

// Species metadata and lookup functions
std::string species_name(std::size_t species_index);
std::size_t species_index_from_name(const std::string& name);
std::size_t num_species();

// Molar mass of a species [g/mol]
double species_molar_mass(std::size_t species_index);
double species_molar_mass_from_name(const std::string& name);

// Convert mole fractions X_k to mass fractions Y_k using molar_masses.
std::vector<double> mole_to_mass(const std::vector<double>& X);

// Convert mass fractions Y_k to mole fractions X_k using molar_masses.
std::vector<double> mass_to_mole(const std::vector<double>& Y);

// Mixture properties
double mwmix(const std::vector<double>& X);

// NASA polynomial evaluations (dimensionless)
// cp_R  : dimensionless heat capacity Cp/R
// h_RT  : dimensionless enthalpy  H/(R*T)
// s_R   : dimensionless entropy   S/R at reference pressure and pure species
// g_over_RT : dimensionless Gibbs free energy G/(R*T) = H/(R*T) - S/R
double cp_R(std::size_t species_idx, double T);
double h_RT(std::size_t species_idx, double T);
double s_R(std::size_t species_idx, double T);
double g_over_RT(std::size_t species_idx, double T);

// Per-species dimensional properties
double cp_species(std::size_t species_idx, double T);  // [J/(mol·K)]
double h_species(std::size_t species_idx, double T);   // [J/mol]
double s_species(std::size_t species_idx, double T);   // [J/(mol·K)] at P_ref

// Thermodynamic properties
double cp(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
double cv(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);  // Internal energy [J/mol]
double density(double T, double P, const std::vector<double>& X);
double molar_volume(double T, double P);  // Ideal gas molar volume [m³/mol]
double specific_gas_constant(const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);

// Mass-specific thermodynamic properties (symmetric interface to molar-basis functions)
double cp_mass(double T, const std::vector<double>& X);  // [J/(kg·K)]
double cv_mass(double T, const std::vector<double>& X);  // [J/(kg·K)]
double h_mass(double T, const std::vector<double>& X);   // [J/kg]
double s_mass(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);  // [J/(kg·K)]
double u_mass(double T, const std::vector<double>& X);   // [J/kg]

// -------------------------------------------------------------
// Air Properties Bundle
// -------------------------------------------------------------

// Bundle of air properties for convenient access
// All properties computed from (T, P, humidity) in a single call
struct AirProperties {
    double rho;      // Density [kg/m³]
    double mu;       // Dynamic viscosity [Pa·s]
    double k;        // Thermal conductivity [W/(m·K)]
    double cp;       // Specific heat at constant pressure [J/(kg·K)]
    double cv;       // Specific heat at constant volume [J/(kg·K)]
    double Pr;       // Prandtl number [-]
    double nu;       // Kinematic viscosity [m²/s]
    double alpha;    // Thermal diffusivity [m²/s]
    double gamma;    // Heat capacity ratio [-]
    double a;        // Speed of sound [m/s]
};

// Compute all air properties at once
// Parameters:
//   T        : temperature [K]
//   P        : pressure [Pa]
//   humidity : relative humidity [0-1] (default: 0.0 = dry air)
// Returns: AirProperties struct with all properties
AirProperties air_properties(double T, double P, double humidity = 0.0);

// -------------------------------------------------------------
// Thermodynamic State Bundle
// -------------------------------------------------------------

// Bundle of thermodynamic properties for a gas mixture
// All properties computed from (T, P, X) in a single call
struct ThermoState {
    double T;         // Temperature [K] (input, echoed back)
    double P;         // Pressure [Pa] (input, echoed back)
    double rho;       // Density [kg/m³]
    double cp;        // Specific heat at constant pressure [J/(mol·K)]
    double cv;        // Specific heat at constant volume [J/(mol·K)]
    double h;         // Specific enthalpy [J/mol]
    double s;         // Specific entropy [J/(mol·K)]
    double u;         // Specific internal energy [J/mol]
    double gamma;     // Isentropic expansion coefficient [-]
    double a;         // Speed of sound [m/s]
    double cp_mass;   // Mass-specific cp [J/(kg·K)]
    double cv_mass;   // Mass-specific cv [J/(kg·K)]
    double h_mass;    // Mass-specific enthalpy [J/kg]
    double s_mass;    // Mass-specific entropy [J/(kg·K)]
    double u_mass;    // Mass-specific internal energy [J/kg]
    double mw;        // Molecular weight [g/mol]
};

// Compute all thermodynamic properties at once
// Parameters:
//   T     : temperature [K]
//   P     : pressure [Pa]
//   X     : mole fractions [-]
//   P_ref : reference pressure for entropy [Pa] (default: 101325.0)
// Returns: ThermoState struct with all properties
ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// -------------------------------------------------------------
// Complete State Bundle (Thermo + Transport)
// -------------------------------------------------------------

// Bundle of thermodynamic and transport properties for a gas mixture
// Combines ThermoState + TransportState in a single call
struct CompleteState {
    ThermoState thermo;       // All 16 thermodynamic properties
    TransportState transport; // All 9 transport properties
};

// Compute all thermodynamic and transport properties at once
// Parameters:
//   T     : temperature [K]
//   P     : pressure [Pa]
//   X     : mole fractions [-]
//   P_ref : reference pressure for entropy [Pa] (default: 101325.0)
// Returns: CompleteState struct with nested thermo and transport states
CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);

// Derivatives of thermodynamic properties with respect to temperature
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
double dg_over_RT_dT(double T, const std::vector<double>& X);

// Inverse calculations (solving for temperature)
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions);

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions);

// -------------------------------------------------------------
// State-based overloads
// -------------------------------------------------------------

double mwmix(const State& s);
double cp(const State& s);
double h(const State& s);
double s(const State& s);
double cv(const State& s);
double u(const State& s);
double density(const State& s);
double specific_gas_constant(const State& s);
double isentropic_expansion_coefficient(const State& s);
double speed_of_sound(const State& s);

#endif // THERMO_H
