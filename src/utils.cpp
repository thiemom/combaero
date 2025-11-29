#include "../include/utils.h"
#include "../include/thermo.h"
#include <cmath>
#include <iostream>
#include <numeric>

// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X) {
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        std::cout << "Warning: Mole fractions sum to " << sum << ", not 1.0" << std::endl;
    }

    std::cout << "\nMixture Composition:" << std::endl;
    std::cout << "-------------------" << std::endl;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            std::cout << species_name(i) << ": " << X[i] << std::endl;
        }
    }

    std::cout << "\nThermodynamic Properties at T = " << T << " K, P = " << P << " Pa:" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Molecular Weight: " << mwmix(X) << " g/mol" << std::endl;
    std::cout << "Density: " << density(T, P, X) << " kg/m³" << std::endl;
    std::cout << "Specific Gas Constant (Rs): " << specific_gas_constant(X) << " J/(kg·K)" << std::endl;
    std::cout << "Isentropic Expansion Coefficient (gamma): " << isentropic_expansion_coefficient(T, X) << std::endl;
    std::cout << "Speed of Sound: " << speed_of_sound(T, X) << " m/s" << std::endl;
    std::cout << "Enthalpy: " << h(T, X) << " J/mol" << std::endl;
    std::cout << "Entropy: " << s(T, X, P) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cp): " << cp(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cv): " << cv(T, X) << " J/(mol·K)" << std::endl;

    std::cout << "\nDerivatives:" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "dh/dT: " << dh_dT(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "ds/dT: " << ds_dT(T, X) << " J/(mol·K²)" << std::endl;
    std::cout << "dCp/dT: " << dcp_dT(T, X) << " J/(mol·K²)" << std::endl;
}
