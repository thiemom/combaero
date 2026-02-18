#ifndef HUMIDAIR_H
#define HUMIDAIR_H

#include <vector>
#include <string>
#include <unordered_map>

// Constants for dry air composition (mole fractions)
extern const std::unordered_map<std::string, double> dry_air_composition;

// Get standard dry air composition as a vector in the order defined by species_index in thermo_transport_data.h
std::vector<double> standard_dry_air_composition();

// Water vapor saturation pressure using Huang (2018) formula [Pa]
// Validated range: -100°C ≤ T ≤ 100°C (173.15–373.15 K)
// Extrapolation allowed within 5% of range boundaries; throws beyond that.
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
// Valid range: same as saturation_vapor_pressure (-100°C to 100°C)
double relative_humidity_from_dewpoint(double T, double Tdp, double P);  // T in K, Tdp in K, P in Pa

// Calculate wet-bulb temperature [K]
double wet_bulb_temperature(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate enthalpy of humid air using psychrometric constants [J/kg dry air]
// Uses cp_air=1005, h_fg=2501000, cp_wv=1860 J/(kg·K) — fast, standard ASHRAE formula
double humid_air_enthalpy(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate enthalpy of humid air using NASA-9 thermodynamic data [J/kg dry air]
// More accurate than humid_air_enthalpy() at non-ambient temperatures.
// Requires T in NASA-9 valid range (200–6000 K).
double humid_air_enthalpy_nasa9(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

// Calculate density of humid air [kg/m³]
double humid_air_density(double T, double P, double RH);  // T in K, P in Pa, RH as fraction (0-1)

#endif // HUMIDAIR_H
