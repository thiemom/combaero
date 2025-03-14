#ifndef HUMIDAIR_H
#define HUMIDAIR_H

#include <vector>
#include <string>
#include <unordered_map>

// Constants for dry air composition (mole fractions)
extern const std::unordered_map<std::string, double> dry_air_composition;

// Get standard dry air composition as a vector in the order defined by species_index in thermo_transport_data.h
std::vector<double> standard_dry_air_composition();

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions);

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions);

// Water vapor saturation pressure using Hyland-Wexler equation [Pa]
// Valid for 0°C ≤ T ≤ 200°C
double saturation_vapor_pressure(double T);  // T in K

// Calculate actual vapor pressure from relative humidity [Pa]
double vapor_pressure(double T, double RH);  // T in K, RH as fraction (0-1)

// Calculate humidity ratio (kg water vapor per kg dry air)
double humidity_ratio(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate mole fraction of water vapor in humid air
double water_vapor_mole_fraction(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate humid air composition (mole fractions)
// Returns a vector of mole fractions in the order defined by species_index in thermo_transport_data.h
std::vector<double> humid_air_composition(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate dewpoint temperature from T, P, RH [K]
double dewpoint(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate relative humidity from dewpoint temperature [fraction]
double relative_humidity_from_dewpoint(double T, double Tdp, double P);  // T in K, Tdp in K, P in Pa

// Calculate wet-bulb temperature [K]
double wet_bulb_temperature(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate enthalpy of humid air [J/kg]
double humid_air_enthalpy(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate density of humid air [kg/m³]
double humid_air_density(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

#endif // HUMIDAIR_H
