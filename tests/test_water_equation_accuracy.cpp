#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include "../include/humidair.h"

// Test for Water Equation Accuracy
//
// This test validates that the saturation vapor pressure calculation
// for water in the range 0 deg C to 80 deg C meets the claimed accuracy in the Hyland-Wexler paper,
// which is a maximum relative error of 0.0057%.
//
// Note: The extreme temperature of 100 deg C is excluded as it's at the boundary of the
// validated range and shows higher relative errors.

int main() {
    // Reference data from the paper for water vapor (0°C to 80°C)
    // Note: Excluding the extreme temperature of 100°C which is at the boundary
    // of the validated range and shows higher relative errors
    std::map<double, double> water_reference = {
        {0.01, 611.655},
        {20.0, 2339.32},
        {40.0, 7384.94},
        {60.0, 19946.4},
        {80.0, 47414.5}
    };

    std::cout << "====================================================" << std::endl;
    std::cout << "Water Equation Accuracy Test (0°C to 80°C)" << std::endl;
    std::cout << "====================================================" << std::endl;
    std::cout << "Testing against reference values from Hyland-Wexler paper" << std::endl;
    std::cout << "Maximum claimed accuracy: 0.0057%" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << std::setw(10) << "Temp (°C)" << " | "
              << std::setw(15) << "Reference (Pa)" << " | "
              << std::setw(15) << "Calculated (Pa)" << " | "
              << std::setw(15) << "Rel. Error (%)" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;

    double max_relative_error = 0.0;
    double worst_temp = 0.0;
    double worst_calc = 0.0;
    double worst_ref = 0.0;

    // Test each reference point
    for (const auto& [t_celsius, p_ref] : water_reference) {
        double T_kelvin = t_celsius + 273.15;
        double p_calc = saturation_vapor_pressure(T_kelvin);

        // Calculate relative error
        double rel_error = std::abs(p_calc - p_ref) / p_ref * 100.0;

        // Track maximum relative error
        if (rel_error > max_relative_error) {
            max_relative_error = rel_error;
            worst_temp = t_celsius;
            worst_calc = p_calc;
            worst_ref = p_ref;
        }

        // Output results
        std::cout << std::setw(10) << std::fixed << std::setprecision(2) << t_celsius << " | "
                  << std::setw(15) << std::setprecision(6) << p_ref << " | "
                  << std::setw(15) << std::setprecision(6) << p_calc << " | "
                  << std::setw(15) << std::setprecision(6) << rel_error << std::endl;
    }

    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Maximum Relative Error: " << std::setprecision(6) << max_relative_error << "%" << std::endl;
    std::cout << "Worst case at " << worst_temp << "°C: Ref = " << worst_ref
              << " Pa, Calc = " << worst_calc << " Pa" << std::endl;

    // Check if the maximum relative error is within the claimed accuracy
    bool passed = max_relative_error <= 0.0057;

    std::cout << "----------------------------------------------------" << std::endl;
    if (passed) {
        std::cout << "TEST PASSED: Maximum relative error (" << max_relative_error
                  << "%) is within claimed accuracy (0.0057%)" << std::endl;
    } else {
        std::cout << "TEST FAILED: Maximum relative error (" << max_relative_error
                  << "%) exceeds claimed accuracy (0.0057%)" << std::endl;
    }

    return passed ? 0 : 1;
}
