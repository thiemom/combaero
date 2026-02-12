#include "../include/utils.h"
#include "../include/thermo.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

// -------------------------------------------------------------
// Pipe Roughness Database
// -------------------------------------------------------------

// Absolute roughness values for common pipe materials [m]
//
// References:
// - Moody, L.F. (1944). "Friction factors for pipe flow".
//   Transactions of the ASME, 66(8), 671-684.
// - Colebrook, C.F. (1939). "Turbulent flow in pipes, with particular
//   reference to the transition region between smooth and rough pipe laws".
//   Journal of the Institution of Civil Engineers, 11(4), 133-156.
// - White, F.M. (2011). "Fluid Mechanics" (7th ed.). McGraw-Hill.
//   Table 6.1, p. 357.
// - Munson, B.R., et al. (2013). "Fundamentals of Fluid Mechanics" (7th ed.).
//   Wiley. Table 8.1, p. 428.
// - Crane Co. (2009). "Flow of Fluids Through Valves, Fittings, and Pipe".
//   Technical Paper No. 410, Table A-24.
//
// Note: Values represent typical absolute roughness (ε) for new, clean pipes.
// Actual roughness may vary with age, corrosion, and service conditions.
//
static const std::unordered_map<std::string, double> ROUGHNESS_DATA = {
    // Smooth surfaces
    {"smooth",              0.0},           // Ideally smooth (theoretical)
    {"drawn_tubing",        1.5e-6},        // Drawn brass, copper, glass, plastic (White 2011)
    {"pvc",                 1.5e-6},        // PVC, polyethylene (White 2011)
    {"plastic",             1.5e-6},        // Generic plastic pipe

    // Steel pipes
    {"commercial_steel",    4.5e-5},        // New commercial steel (Moody 1944, White 2011)
    {"new_steel",           4.5e-5},        // Alias for commercial_steel
    {"wrought_iron",        4.5e-5},        // Wrought iron (White 2011)
    {"galvanized_iron",     1.5e-4},        // Galvanized iron/steel (Moody 1944, Crane 2009)
    {"galvanized_steel",    1.5e-4},        // Alias for galvanized_iron
    {"rusted_steel",        2.5e-4},        // Moderately rusted steel (Munson 2013)

    // Cast iron
    {"cast_iron",           2.6e-4},        // Uncoated cast iron (Moody 1944, White 2011)
    {"asphalted_cast_iron", 1.2e-4},        // Asphalted cast iron (White 2011)

    // Concrete
    {"concrete",            3.0e-4},        // Concrete, well-finished (Moody 1944, Munson 2013)
    {"rough_concrete",      3.0e-3},        // Rough concrete (White 2011)

    // Other materials
    {"riveted_steel",       9.0e-4},        // Riveted steel (Moody 1944, Crane 2009)
    {"wood_stave",          1.8e-4},        // Wood stave pipe (White 2011)
    {"corrugated_metal",    4.5e-2},        // Corrugated metal pipe (Munson 2013)
};

std::unordered_map<std::string, double> standard_pipe_roughness() {
    return ROUGHNESS_DATA;
}

double pipe_roughness(const std::string& material) {
    // Convert to lowercase for case-insensitive lookup
    std::string material_lower = material;
    std::transform(material_lower.begin(), material_lower.end(),
                   material_lower.begin(), ::tolower);

    auto it = ROUGHNESS_DATA.find(material_lower);
    if (it != ROUGHNESS_DATA.end()) {
        return it->second;
    }

    // Material not found - provide helpful error message
    throw std::invalid_argument(
        "pipe_roughness: unknown material '" + material + "'. "
        "Use standard_pipe_roughness() to see available materials.");
}

// -------------------------------------------------------------
// Utility Functions
// -------------------------------------------------------------

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
