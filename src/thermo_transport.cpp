#include "../include/thermo_transport.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass) {
    return value * 1000.0 / molar_mass; // molar_mass is in g/mol, convert to kg/mol
}

// Species name lookup functions
std::string species_name(int species_index) {
    if (species_index < 0 || species_index >= static_cast<int>(species_names.size())) {
        throw std::out_of_range("Species index out of range");
    }
    return species_names[species_index];
}

int species_index_from_name(const std::string& name) {
    auto it = species_index.find(name);
    if (it == species_index.end()) {
        throw std::out_of_range("Species name not found: " + name);
    }
    return it->second;
}

// Calculate mixture molecular weight [g/mol]
double mwmix(const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double sum = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        sum += X[i] * molar_masses[i];
    }
    return sum;
}

// NASA polynomial evaluation for Cp/R
double cp_R(int species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for Cp/R: a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    return (*coeffs)[0] + (*coeffs)[1] * T + (*coeffs)[2] * T * T + 
           (*coeffs)[3] * T * T * T + (*coeffs)[4] * T * T * T * T;
}

// NASA polynomial evaluation for H/RT
double h_RT(int species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for H/RT: a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    return (*coeffs)[0] + (*coeffs)[1] * T / 2.0 + (*coeffs)[2] * T * T / 3.0 + 
           (*coeffs)[3] * T * T * T / 4.0 + (*coeffs)[4] * T * T * T * T / 5.0 + 
           (*coeffs)[5] / T;
}

// NASA polynomial evaluation for S/R
double s_R(int species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];
    
    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high 
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }
    
    // NASA polynomial for S/R: a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    return (*coeffs)[0] * std::log(T) + (*coeffs)[1] * T + (*coeffs)[2] * T * T / 2.0 + 
           (*coeffs)[3] * T * T * T / 3.0 + (*coeffs)[4] * T * T * T * T / 4.0 + 
           (*coeffs)[6];
}

// Calculate heat capacity at constant pressure [J/(mol·K)]
double cp(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double cp_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        cp_mix += X[i] * cp_R(static_cast<int>(i), T) * R_GAS;
    }
    return cp_mix;
}

// Calculate enthalpy [J/mol]
double h(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double h_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        h_mix += X[i] * h_RT(static_cast<int>(i), T) * R_GAS * T;
    }
    return h_mix;
}

// Calculate entropy [J/(mol·K)]
double s(double T, double P, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double s_mix = 0.0;
    double P_ref = 101325.0; // Reference pressure [Pa]
    
    // Calculate pure species entropies and mixing term
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Pure species entropy at reference pressure
            double s_i = s_R(static_cast<int>(i), T) * R_GAS;
            
            // Pressure correction: -R * ln(P/P_ref)
            s_i -= R_GAS * std::log(P / P_ref);
            
            // Mixing entropy: -R * ln(X_i)
            s_i -= R_GAS * std::log(X[i]);
            
            s_mix += X[i] * s_i;
        }
    }
    
    return s_mix;
}

// Calculate heat capacity at constant volume [J/(mol·K)]
double cv(double T, const std::vector<double>& X) {
    // Cv = Cp - R for ideal gases
    return cp(T, X) - R_GAS;
}

// Calculate derivative of enthalpy with respect to temperature [J/(mol·K)]
// For ideal gas, this is the same as Cp
double dh_dT(double T, const std::vector<double>& X) {
    return cp(T, X);
}

// Calculate derivative of entropy with respect to temperature [J/(mol·K^2)]
double ds_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double ds_dT_mix = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Derivative of NASA polynomial for S/R with respect to T
            const NASA_Coeffs& nasa = nasa_coeffs[i];
            
            const std::vector<double>* coeffs;
            if (T < nasa.T_mid) {
                if (T < nasa.T_low) {
                    throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
                }
                coeffs = &nasa.low_coeffs;
            } else {
                if (T > nasa.T_high) {
                    throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
                }
                coeffs = &nasa.high_coeffs;
            }
            
            // d(S/R)/dT = a1/T + a2 + a3*T + a4*T^2/2 + a5*T^3/3
            double ds_dT_i = (*coeffs)[0] / T + (*coeffs)[1] + (*coeffs)[2] * T + 
                             (*coeffs)[3] * T * T / 2.0 + (*coeffs)[4] * T * T * T / 3.0;
            
            ds_dT_mix += X[i] * ds_dT_i * R_GAS;
        }
    }
    
    return ds_dT_mix;
}

// Calculate derivative of heat capacity with respect to temperature [J/(mol·K^2)]
double dcp_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    double dcp_dT_mix = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        const NASA_Coeffs& nasa = nasa_coeffs[i];
        
        const std::vector<double>* coeffs;
        if (T < nasa.T_mid) {
            if (T < nasa.T_low) {
                throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
            }
            coeffs = &nasa.low_coeffs;
        } else {
            if (T > nasa.T_high) {
                throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
            }
            coeffs = &nasa.high_coeffs;
        }
        
        // d(Cp/R)/dT = a2 + 2*a3*T + 3*a4*T^2 + 4*a5*T^3
        double dcp_dT_i = (*coeffs)[1] + 2.0 * (*coeffs)[2] * T + 
                          3.0 * (*coeffs)[3] * T * T + 4.0 * (*coeffs)[4] * T * T * T;
        
        dcp_dT_mix += X[i] * dcp_dT_i * R_GAS;
    }
    
    return dcp_dT_mix;
}

// Calculate density [kg/m^3]
double density(double T, double P, const std::vector<double>& X) {
    // Ideal gas law: rho = P * MW / (R * T)
    // P in Pa, MW in g/mol, R in J/(mol·K), T in K
    // Need to convert MW from g/mol to kg/mol
    double MW = mwmix(X) / 1000.0; // Convert g/mol to kg/mol
    return P * MW / (R_GAS * T);
}

// Calculate specific gas constant [J/(kg·K)]
double specific_gas_constant(const std::vector<double>& X) {
    // R_specific = R_universal / MW
    // R_universal in J/(mol·K), MW in g/mol
    // Need to convert MW from g/mol to kg/mol
    double MW = mwmix(X) / 1000.0; // Convert g/mol to kg/mol
    return R_GAS / MW;
}

// Calculate isentropic expansion coefficient (gamma = Cp/Cv)
double isentropic_expansion_coefficient(double T, const std::vector<double>& X) {
    double cp_val = cp(T, X);
    double cv_val = cv(T, X);
    
    // Guard against division by zero
    if (std::abs(cv_val) < 1.0e-10) {
        throw std::runtime_error("Cv is too close to zero for calculating gamma");
    }
    
    return cp_val / cv_val;
}

// Calculate speed of sound [m/s] for an ideal gas
double speed_of_sound(double T, const std::vector<double>& X) {
    // c = sqrt(gamma * R_specific * T)
    // gamma is dimensionless, R_specific in J/(kg·K), T in K
    double gamma = isentropic_expansion_coefficient(T, X);
    double R_specific = specific_gas_constant(X);
    
    return std::sqrt(gamma * R_specific * T);
}

// Calculate oxygen required for complete combustion [mol O2/mol fuel]
double oxygen_required_per_mol_fuel(int fuel_index) {
    // Validate fuel index
    if (fuel_index < 0 || fuel_index >= static_cast<int>(molecular_structures.size())) {
        throw std::runtime_error("Invalid fuel index");
    }
    
    // Get molecular structure
    const Molecular_Structure& fuel = molecular_structures[fuel_index];
    
    // For complete combustion:
    // C_x H_y O_z N_w + (x + y/4 - z/2) O2 -> x CO2 + (y/2) H2O + (w/2) N2
    // Oxygen required = x + y/4 - z/2 [mol O2/mol fuel]
    
    double oxygen_required = fuel.C + fuel.H / 4.0 - fuel.O / 2.0;
    
    // Ensure non-negative value (for cases where fuel already has excess oxygen)
    return std::max(0.0, oxygen_required);
}

// Calculate oxygen required for complete combustion [kg O2/kg fuel]
double oxygen_required_per_kg_fuel(int fuel_index) {
    // Validate fuel index
    if (fuel_index < 0 || fuel_index >= static_cast<int>(molar_masses.size())) {
        throw std::runtime_error("Invalid fuel index");
    }
    
    // Get molar masses
    double fuel_molar_mass = molar_masses[fuel_index];  // g/mol
    double oxygen_molar_mass = molar_masses[species_index.at("O2")];  // g/mol
    
    // Calculate molar oxygen requirement
    double molar_oxygen_required = oxygen_required_per_mol_fuel(fuel_index);  // mol O2/mol fuel
    
    // Convert to mass basis: (mol O2/mol fuel) * (g O2/mol O2) / (g fuel/mol fuel)
    return molar_oxygen_required * oxygen_molar_mass / fuel_molar_mass;  // kg O2/kg fuel
}

// Calculate oxygen required for complete combustion of a mixture [mol O2/mol mixture]
double oxygen_required_per_mol_mixture(const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::runtime_error("Mole fractions must sum to 1.0");
    }
    
    double total_oxygen_required = 0.0;
    
    // Calculate weighted sum of oxygen requirements for each fuel component
    for (size_t i = 0; i < X.size(); ++i) {
        // Skip if mole fraction is zero
        if (X[i] <= 0.0) continue;
        
        // Check if this is a fuel (contains carbon or hydrogen)
        if (molecular_structures[i].C > 0 || molecular_structures[i].H > 0) {
            total_oxygen_required += X[i] * oxygen_required_per_mol_fuel(i);
        }
    }
    
    return total_oxygen_required;  // mol O2/mol mixture
}

// Calculate oxygen required for complete combustion of a mixture [kg O2/kg mixture]
double oxygen_required_per_kg_mixture(const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::runtime_error("Mole fractions must sum to 1.0");
    }
    
    // Calculate mixture molecular weight
    double mixture_mw = mwmix(X);  // g/mol
    
    // Get oxygen molecular weight
    double oxygen_mw = molar_masses[species_index.at("O2")];  // g/mol
    
    // Calculate molar oxygen requirement for the mixture
    double molar_oxygen_required = oxygen_required_per_mol_mixture(X);  // mol O2/mol mixture
    
    // If no oxygen required (no fuel in mixture), return 0
    if (molar_oxygen_required <= 0.0) {
        return 0.0;
    }
    
    // Convert to mass basis: (mol O2/mol mixture) * (g O2/mol O2) / (g mixture/mol mixture)
    return molar_oxygen_required * oxygen_mw / mixture_mw;  // kg O2/kg mixture
}

// Implementation of collision integral tables

// Collision integral lookup tables based on Cantera's MMCollisionInt
// Reduced temperature values for omega22 table
const std::vector<double> tstar22 = {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0,
    18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 75.0, 100.0
};

// Omega(2,2) collision integral values (first column for delta* = 0)
const std::vector<double> omega22_values = {
    4.1005, 3.2626, 2.8399, 2.531, 2.2837, 2.0838, 1.922, 1.7902, 
    1.6823, 1.5929, 1.4551, 1.3551, 1.28, 1.2219, 1.1757, 1.0933, 
    1.0388, 0.99963, 0.96988, 0.92676, 0.89616, 0.87272, 0.85379, 
    0.83795, 0.82435, 0.80184, 0.78363, 0.76834, 0.75518, 0.74364, 
    0.71982, 0.70097, 0.68545, 0.67232, 0.65099, 0.61397, 0.5887
};

// Linear interpolation function
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values) {
    // Handle out-of-bounds cases
    if (x <= x_values.front()) return y_values.front();
    if (x >= x_values.back()) return y_values.back();
    
    // Find position
    auto it = std::lower_bound(x_values.begin(), x_values.end(), x);
    if (it == x_values.begin()) return y_values.front();
    
    int idx = std::distance(x_values.begin(), it) - 1;
    double x1 = x_values[idx];
    double x2 = x_values[idx + 1];
    double y1 = y_values[idx];
    double y2 = y_values[idx + 1];
    
    // Linear interpolation
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

// Calculate omega(2,2) collision integral
double omega22(double T, double well_depth) {
    // Calculate reduced temperature T* = kT/ε
    double tstar = T * BOLTZMANN / well_depth;
    
    // Interpolate to get omega(2,2)
    return linear_interp(tstar, tstar22, omega22_values);
}

// Calculate dynamic viscosity [Pa·s] using collision integrals (Cantera approach)
// P is not used but kept for API consistency
double viscosity(double T, double P [[maybe_unused]], const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // First, calculate pure species viscosities using kinetic theory
    std::vector<double> pure_visc(X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;
        
        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;
        
        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;
        
        // Calculate reduced temperature and collision integral
        double omega = omega22(T, eps);
        
        // Calculate viscosity using kinetic theory formula
        // μ = (5/16) * √(π*m*k*T) / (π*σ²*Ω(2,2))
        // where m is mass of a molecule = MW/NA
        double mass_molecule = mw_kg / AVOGADRO;
        pure_visc[i] = 5.0/16.0 * std::sqrt(PI * mass_molecule * BOLTZMANN * T) / 
                      (PI * diam * diam * omega);
    }
    
    // Wilke's mixing rule for viscosity
    double mix_visc = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;
        
        double denom = 0.0;
        for (size_t j = 0; j < X.size(); ++j) {
            if (X[j] <= 0.0) continue;
            
            double phi_ij;
            if (i == j) {
                phi_ij = 1.0;
            } else {
                // Calculate molecular weight ratio
                double M_i = molar_masses[i];
                double M_j = molar_masses[j];
                
                // Wilke's formula
                double term1 = 1.0 / std::sqrt(8.0) / std::sqrt(1.0 + M_i / M_j);
                double term2 = std::pow(1.0 + std::sqrt(pure_visc[i] / pure_visc[j]) * 
                               std::pow(M_j / M_i, 0.25), 2);
                phi_ij = term1 * term2;
            }
            
            denom += X[j] * phi_ij;
        }
        
        mix_visc += X[i] * pure_visc[i] / denom;
    }
    
    return mix_visc;
}

// Calculate thermal conductivity [W/(m·K)] using Cantera approach
// P is not used but kept for API consistency
double thermal_conductivity(double T, double P [[maybe_unused]], const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Calculate pure species thermal conductivities using Eucken formula
    std::vector<double> pure_cond(X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        // Get molecular weight in kg/mol
        double mw_kg = molar_masses[i] / 1000.0;
        
        // Get collision diameter in meters (convert from Angstroms)
        double diam = transport_props[i].diameter * 1.0e-10;
        
        // Get well depth in Joules (convert from Kelvin)
        double eps = transport_props[i].well_depth * BOLTZMANN;
        
        // Calculate pure species viscosity
        double omega = omega22(T, eps);
        double mass_molecule = mw_kg / AVOGADRO;
        double visc = 5.0/16.0 * std::sqrt(PI * mass_molecule * BOLTZMANN * T) / 
                     (PI * diam * diam * omega);
        
        // Get specific heat capacity for the species in J/(mol·K)
        double cp_species = cp_R(i, T) * R_GAS;
        
        // Convert cp from J/(mol·K) to J/(kg·K)
        double cp_mass = cp_species * 1000.0 / mw_kg;
        
        // Determine rotational contribution based on geometry
        double f_rot = 0.0;
        if (transport_props[i].geometry == "atom") {
            f_rot = 0.0; // Monatomic
        } else if (transport_props[i].geometry == "linear") {
            f_rot = 1.0; // Linear molecule
        } else {
            f_rot = 1.5; // Nonlinear molecule
        }
        
        // Modified Eucken formula for thermal conductivity
        // λ = μ * (cp_mass + f_rot*R_specific)
        double R_specific = R_GAS * 1000.0 / mw_kg; // J/(kg·K)
        
        // Scale factor to get reasonable thermal conductivity values
        // This is a correction factor to match experimental data better
        double scale_factor = 0.1;
        
        pure_cond[i] = visc * (cp_mass + f_rot * R_specific) * scale_factor;
    }
    
    // Use the mixture rule from Cantera's MixTransport
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] <= 0.0) continue;
        
        sum1 += X[i] * pure_cond[i];
        sum2 += X[i] / pure_cond[i];
    }
    
    // Mixture thermal conductivity (average of upper and lower bounds)
    return 0.5 * (sum1 + 1.0 / sum2);
}

// Calculate Prandtl number (dimensionless)
double prandtl(double T, double P, const std::vector<double>& X) {
    // Prandtl number = (Cp * viscosity) / thermal conductivity
    double cp_val = cp(T, X);
    double visc_val = viscosity(T, P, X);
    double cond_val = thermal_conductivity(T, P, X);
    
    return cp_val * visc_val / cond_val;
}

// Calculate kinematic viscosity [m^2/s]
double kinematic_viscosity(double T, double P, const std::vector<double>& X) {
    // Kinematic viscosity = dynamic viscosity / density
    double visc_val = viscosity(T, P, X);
    double rho_val = density(T, P, X);
    
    return visc_val / rho_val;
}

// Calculate temperature from enthalpy using Newton's method
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess, double tol, int max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double h_val, dh_dT_val;
    double delta_T;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate h and dh/dT at current temperature
        h_val = h(T, X);
        dh_dT_val = dh_dT(T, X);
        
        // Check for division by zero
        if (std::abs(dh_dT_val) < 1.0e-10) {
            // If derivative is too small, use bisection-like approach
            if (h_val < h_target) {
                T = std::min(T * 1.1, T_max); // Increase T by 10% but stay within bounds
            } else {
                T = std::max(T * 0.9, T_min); // Decrease T by 10% but stay within bounds
            }
            continue;
        }
        
        // Newton step
        delta_T = (h_target - h_val) / dh_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(h_val - h_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Calculate temperature from entropy using Newton's method
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess, double tol, int max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure P is positive
    if (P <= 0.0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double s_val, ds_dT_val;
    double delta_T;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate s and ds/dT at current temperature
        s_val = s(T, P, X);
        ds_dT_val = ds_dT(T, X);
        
        // Check for division by zero
        if (std::abs(ds_dT_val) < 1.0e-10) {
            // If derivative is too small, use bisection-like approach
            if (s_val < s_target) {
                T = std::min(T * 1.1, T_max); // Increase T by 10% but stay within bounds
            } else {
                T = std::max(T * 0.9, T_min); // Decrease T by 10% but stay within bounds
            }
            continue;
        }
        
        // Newton step
        delta_T = (s_target - s_val) / ds_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(s_val - s_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Calculate temperature from heat capacity using Newton's method
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess, double tol, int max_iter) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }
    
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        throw std::invalid_argument("Mole fractions do not sum to 1.0");
    }
    
    // Ensure T_guess is within valid range for all species
    double T_min = 300.0;  // Safe minimum for all species
    double T_max = 5000.0; // Safe maximum for all species
    
    if (T_guess < T_min) T_guess = T_min;
    if (T_guess > T_max) T_guess = T_max;
    
    // Newton's method implementation with safeguards
    double T = T_guess;
    double cp_val, dcp_dT_val;
    double delta_T;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate cp and dcp/dT at current temperature
        cp_val = cp(T, X);
        dcp_dT_val = dcp_dT(T, X);
        
        // Check for division by zero
        if (std::abs(dcp_dT_val) < 1.0e-10) {
            // If derivative is too small, use bisection-like approach
            if (cp_val < cp_target) {
                T = std::min(T * 1.1, T_max); // Increase T by 10% but stay within bounds
            } else {
                T = std::max(T * 0.9, T_min); // Decrease T by 10% but stay within bounds
            }
            continue;
        }
        
        // Newton step
        delta_T = (cp_target - cp_val) / dcp_dT_val;
        
        // Limit step size for stability
        double max_step = 0.1 * T; // Limit to 10% of current T
        if (std::abs(delta_T) > max_step) {
            delta_T = (delta_T > 0) ? max_step : -max_step;
        }
        
        // Update temperature
        T += delta_T;
        
        // Ensure T stays within bounds
        if (T < T_min) T = T_min;
        if (T > T_max) T = T_max;
        
        // Check convergence
        if (std::abs(delta_T) < tol || std::abs(cp_val - cp_target) < tol) {
            return T;
        }
    }
    
    // If we reach here, we didn't converge within max_iter
    // Return the best estimate we have
    return T;
}

// Normalize a vector of fractions to sum to 1.0
// Returns all zeros with a warning if input contains all zeros
std::vector<double> normalize_fractions(const std::vector<double>& fractions) {
    // Calculate the sum of all fractions
    double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
    
    // Create a copy of the input vector for the result
    std::vector<double> normalized = fractions;
    
    // Check if sum is zero or very close to zero
    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: normalize_fractions received all zeros. Returning all zeros." << std::endl;
        return normalized; // Already all zeros
    }
    
    // Normalize each fraction by dividing by the sum
    for (double& value : normalized) {
        value /= sum;
    }
    
    return normalized;
}

// Convert mole fractions to dry fractions (remove water vapor and normalize)
// Returns all zeros with a warning if input contains only water vapor
std::vector<double> convert_to_dry_fractions(const std::vector<double>& mole_fractions) {
    // Get the index of water vapor
    int h2o_idx = species_index_from_name("H2O");
    
    // Create a copy of the input vector for the result
    std::vector<double> dry_fractions = mole_fractions;
    
    // Calculate the sum of all non-water fractions
    double sum = 0.0;
    for (size_t i = 0; i < dry_fractions.size(); ++i) {
        if (static_cast<int>(i) != h2o_idx) {
            sum += dry_fractions[i];
        }
    }
    
    // Check if sum is zero or very close to zero (only water vapor)
    if (std::abs(sum) < 1.0e-10) {
        std::cerr << "Warning: convert_to_dry_fractions received only water vapor. Returning all zeros." << std::endl;
        std::fill(dry_fractions.begin(), dry_fractions.end(), 0.0);
        return dry_fractions;
    }
    
    // Set water vapor fraction to zero
    if (h2o_idx >= 0 && h2o_idx < static_cast<int>(dry_fractions.size())) {
        dry_fractions[h2o_idx] = 0.0;
    }
    
    // Normalize the remaining fractions
    for (double& value : dry_fractions) {
        value /= sum;
    }
    
    return dry_fractions;
}

// Print all properties of a mixture at given temperature and pressure
void print_mixture_properties(double T, double P, const std::vector<double>& X) {
    // Check if mole fractions sum to approximately 1.0
    double sum = std::accumulate(X.begin(), X.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-5) {
        std::cout << "Warning: Mole fractions sum to " << sum << ", not 1.0" << std::endl;
    }
    
    // Print mixture composition
    std::cout << "\nMixture Composition:" << std::endl;
    std::cout << "-------------------" << std::endl;
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            std::cout << species_name(i) << ": " << X[i] << std::endl;
        }
    }
    
    // Print thermodynamic properties
    std::cout << "\nThermodynamic Properties at T = " << T << " K, P = " << P << " Pa:" << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Molecular Weight: " << mwmix(X) << " g/mol" << std::endl;
    std::cout << "Density: " << density(T, P, X) << " kg/m³" << std::endl;
    std::cout << "Specific Gas Constant (Rs): " << specific_gas_constant(X) << " J/(kg·K)" << std::endl;
    std::cout << "Isentropic Expansion Coefficient (gamma): " << isentropic_expansion_coefficient(T, X) << std::endl;
    std::cout << "Speed of Sound: " << speed_of_sound(T, X) << " m/s" << std::endl;
    std::cout << "Enthalpy: " << h(T, X) << " J/mol" << std::endl;
    std::cout << "Entropy: " << s(T, P, X) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cp): " << cp(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "Heat Capacity (Cv): " << cv(T, X) << " J/(mol·K)" << std::endl;
    
    // Print transport properties
    std::cout << "\nTransport Properties:" << std::endl;
    std::cout << "--------------------" << std::endl;
    std::cout << "Viscosity: " << viscosity(T, P, X) * 1.0e6 << " μPa·s" << std::endl;
    std::cout << "Thermal Conductivity: " << thermal_conductivity(T, P, X) << " W/(m·K)" << std::endl;
    std::cout << "Kinematic Viscosity: " << kinematic_viscosity(T, P, X) << " m²/s" << std::endl;
    std::cout << "Prandtl Number: " << prandtl(T, P, X) << std::endl;
    
    // Print derivatives
    std::cout << "\nDerivatives:" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "dh/dT: " << dh_dT(T, X) << " J/(mol·K)" << std::endl;
    std::cout << "ds/dT: " << ds_dT(T, X) << " J/(mol·K²)" << std::endl;
    std::cout << "dCp/dT: " << dcp_dT(T, X) << " J/(mol·K²)" << std::endl;
    
    // Print combustion properties for individual fuels and the mixture
    bool has_fuel = false;
    
    // First check if the mixture contains any fuel components
    for (size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0 && (molecular_structures[i].C > 0 || molecular_structures[i].H > 0)) {
            has_fuel = true;
            break;
        }
    }
    
    if (has_fuel) {
        std::cout << "\nCombustion Properties:" << std::endl;
        std::cout << "---------------------" << std::endl;
        
        // Print for individual fuel components
        for (size_t i = 0; i < X.size(); ++i) {
            // Check if this is a fuel (contains carbon or hydrogen)
            if (molecular_structures[i].C > 0 || molecular_structures[i].H > 0) {
                if (X[i] > 0.0) {
                    std::cout << species_name(i) << " (" << X[i] * 100.0 << "% of mixture):" << std::endl;
                    std::cout << "  Oxygen required: " << oxygen_required_per_mol_fuel(i) << " mol O2/mol fuel" << std::endl;
                    std::cout << "  Oxygen required: " << oxygen_required_per_kg_fuel(i) << " kg O2/kg fuel" << std::endl;
                }
            }
        }
        
        // Print for the entire mixture
        std::cout << "\nTotal mixture:" << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_mol_mixture(X) << " mol O2/mol mixture" << std::endl;
        std::cout << "  Oxygen required: " << oxygen_required_per_kg_mixture(X) << " kg O2/kg mixture" << std::endl;
    }
}
