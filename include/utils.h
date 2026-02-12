#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <unordered_map>
#include <vector>

// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X);

// -------------------------------------------------------------
// Pipe Roughness Database
// -------------------------------------------------------------

// Get absolute roughness for a standard pipe material [m]
// Returns roughness value or throws std::invalid_argument if material not found
// Material names are case-insensitive
double pipe_roughness(const std::string& material);

// Get all standard pipe roughness values as a map
// Returns: map of material name -> roughness [m]
std::unordered_map<std::string, double> standard_pipe_roughness();

#endif // UTILS_H
