#include "../include/thermo_transport.h"
#include "../include/thermo_transport_data.h"
#include <iostream>
#include <iomanip>
#include <vector>

void print_fractions(const std::string& title, const std::vector<double>& fractions) {
    std::cout << title << ":" << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "Species | Fraction" << std::endl;
    std::cout << "--------|-------------" << std::endl;
    
    for (size_t i = 0; i < species_names.size(); ++i) {
        std::cout << std::setw(8) << std::left << species_names[i] 
                  << "| " << std::setw(12) << std::fixed << std::setprecision(6) << fractions[i] << std::endl;
    }
    
    double sum = 0.0;
    for (double fraction : fractions) {
        sum += fraction;
    }
    
    std::cout << "\nSum of fractions: " << std::fixed << std::setprecision(6) << sum << std::endl;
    std::cout << std::endl;
}

int main() {
    // Test 1: Normalize standard fractions
    std::cout << "Test 1: Normalize standard fractions" << std::endl;
    std::cout << "===================================" << std::endl;
    
    // Create a vector with non-normalized fractions
    std::vector<double> test1 = {0.0, 0.1, 0.0, 0.4, 0.05, 0.02, 0.03};
    print_fractions("Original fractions", test1);
    
    // Normalize the fractions
    std::vector<double> normalized1 = normalize_fractions(test1);
    print_fractions("Normalized fractions", normalized1);
    
    // Test 2: Normalize all zeros
    std::cout << "Test 2: Normalize all zeros" << std::endl;
    std::cout << "==========================" << std::endl;
    
    std::vector<double> test2(species_names.size(), 0.0);
    print_fractions("Original fractions (all zeros)", test2);
    
    std::vector<double> normalized2 = normalize_fractions(test2);
    print_fractions("Normalized fractions (should remain all zeros with warning)", normalized2);
    
    // Test 3: Convert to dry fractions (standard humid air)
    std::cout << "Test 3: Convert to dry fractions (standard humid air)" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    // Create a vector with humid air (including water vapor)
    std::vector<double> test3 = {0.0, 0.2, 0.0, 0.75, 0.009, 0.001, 0.04};  // 4% water vapor
    print_fractions("Original humid air fractions", test3);
    
    // Convert to dry fractions
    std::vector<double> dry3 = convert_to_dry_fractions(test3);
    print_fractions("Dry air fractions", dry3);
    
    // Test 4: Convert to dry fractions (pure water vapor)
    std::cout << "Test 4: Convert to dry fractions (pure water vapor)" << std::endl;
    std::cout << "================================================" << std::endl;
    
    std::vector<double> test4(species_names.size(), 0.0);
    test4[species_index_from_name("H2O")] = 1.0;  // Pure water vapor
    print_fractions("Original pure water vapor", test4);
    
    std::vector<double> dry4 = convert_to_dry_fractions(test4);
    print_fractions("Dry air fractions (should be all zeros with warning)", dry4);
    
    return 0;
}
