#include "../include/units.h"
#include <cassert>
#include <iostream>

int main() {
    using namespace combaero::units;

    // Test basic lookup
    auto cp_u = get_units("cp");
    assert(cp_u.has_value());
    assert(cp_u->input == "T: K, X: mol/mol");
    assert(cp_u->output == "J/(mol*K)");
    (void)cp_u;

    auto density_u = get_units("density");
    assert(density_u.has_value());
    assert(density_u->output == "kg/m^3");
    (void)density_u;

    // Test convenience functions
    assert(input_units("viscosity") == "T: K, P: Pa, X: mol/mol");
    assert(output_units("viscosity") == "Pa*s");

    // Test not found
    assert(!get_units("nonexistent_function").has_value());
    assert(input_units("nonexistent") == "");
    assert(!has_units("nonexistent"));

    // Test has_units
    assert(has_units("cp"));
    assert(has_units("density"));
    assert(has_units("State::T"));

    // Test iteration
    std::size_t count = 0;
    for (auto it = begin(); it != end(); ++it) {
        ++count;
        assert(!it->name.empty());
    }
    assert(count == function_count);
    assert(count > 100);  // Sanity check: we have many entries

    std::cout << "All unit tests passed! (" << count << " functions registered)\n";

    // Demo: print a few entries
    std::cout << "\nSample entries:\n";
    const char* samples[] = {"cp", "density", "viscosity", "nozzle_flow", "State::rho"};
    for (const char* name : samples) {
        if (auto u = get_units(name)) {
            std::cout << "  " << name << "\n"
                      << "    Input:  " << u->input << "\n"
                      << "    Output: " << u->output << "\n";
        }
    }

    return 0;
}
