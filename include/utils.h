#ifndef UTILS_H
#define UTILS_H

#include <vector>

namespace combaero {
// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X);
} // namespace combaero


#endif // UTILS_H
