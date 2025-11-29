#ifndef TRANSPORT_H
#define TRANSPORT_H

// TODO: check if imports are needed and update to new refactored structure

#include "thermo_transport_data.h"
#include "thermo.h"
#include <iterator>
#include <vector>
#include <string>

// TODO: wrap this in namespace combaero::transport similar to thermo.h

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T, double well_depth);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);


#endif // TRANSPORT_H
