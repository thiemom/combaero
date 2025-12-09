#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>

#include "thermo.h"
#include "transport.h"
#include "combustion.h"
#include "equilibrium.h"
#include "humidair.h"
#include "common_names.h"
#include "compressible.h"
#include "incompressible.h"
#include "friction.h"
#include "heat_transfer.h"
#include "state.h"
#include "orifice.h"
#include "units.h"

namespace py = pybind11;

static std::vector<double> to_vec(
    py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    auto buf = arr.unchecked<1>();
    std::vector<double> v(buf.shape(0));
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        v[static_cast<std::size_t>(i)] = buf(i);
    return v;
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Python bindings for combaero core";

    // Species metadata helpers
    m.def("num_species", &num_species, "Number of thermo species in the internal tables.");

    m.def(
        "species_name",
        &species_name,
        py::arg("index"),
        "Return the canonical species name for a given index."
    );

    m.def(
        "species_index_from_name",
        &species_index_from_name,
        py::arg("name"),
        "Return the species index for a given canonical name."
    );

    m.def(
        "species_molar_mass",
        &species_molar_mass,
        py::arg("index"),
        "Return the molar mass [g/mol] for a given species index."
    );

    m.def(
        "species_molar_mass_from_name",
        &species_molar_mass_from_name,
        py::arg("name"),
        "Return the molar mass [g/mol] for a given species name."
    );

    // Common names for species symbols
    m.def(
        "formula_to_name",
        []() {
            return combaero::formula_to_name;
        },
        "Mapping from canonical species symbols (e.g. 'CH4') to human-readable common names (e.g. 'Methane')."
    );

    m.def(
        "name_to_formula",
        []() {
            return combaero::name_to_formula;
        },
        "Mapping from human-readable common names (e.g. 'Methane') to canonical species symbols (e.g. 'CH4')."
    );

    m.def(
        "common_name",
        &combaero::common_name,
        py::arg("formula"),
        "Get common name from formula (e.g., 'CH4' -> 'Methane')."
    );

    m.def(
        "formula",
        &combaero::formula,
        py::arg("name"),
        "Get formula from common name (e.g., 'Methane' -> 'CH4')."  
    );

    m.def(
        "mixture_h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy at temperature T and mole fractions X."
    );

    // Thermodynamic mixture properties
    m.def(
        "cp",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cp(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isobaric heat capacity Cp(T, X) [J/(mol·K)]."
    );

    m.def(
        "h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy h(T, X) [J/mol]."
    );

    m.def(
        "s",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double P,
           double P_ref)
        {
            auto X = to_vec(X_arr);
            return s(T, X, P, P_ref);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P"),
        py::arg("P_ref") = 101325.0,
        "Mixture entropy s(T, X, P, P_ref) [J/(mol·K)]."
    );

    m.def(
        "cv",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cv(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isochoric heat capacity Cv(T, X) [J/(mol·K)]."
    );

    m.def(
        "density",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return density(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Mixture density rho(T, P, X) [kg/m^3]."
    );

    m.def(
        "specific_gas_constant",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return specific_gas_constant(X);
        },
        py::arg("X"),
        "Mixture specific gas constant R_s(X) [J/(kg·K)]."
    );

    m.def(
        "isentropic_expansion_coefficient",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return isentropic_expansion_coefficient(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Isentropic expansion coefficient gamma(T, X) = Cp/Cv."
    );

    m.def(
        "speed_of_sound",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return speed_of_sound(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Speed of sound c(T, X) [m/s]."
    );

    m.def(
        "molar_volume",
        &molar_volume,
        py::arg("T"),
        py::arg("P"),
        "Ideal gas molar volume V_m = R*T/P [m³/mol]."
    );

    // Transport properties
    m.def(
        "viscosity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return viscosity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Dynamic viscosity mu(T, P, X) [Pa·s]."
    );

    m.def(
        "thermal_conductivity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return thermal_conductivity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Thermal conductivity lambda(T, P, X) [W/(m·K)]."
    );

    m.def(
        "kinematic_viscosity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return kinematic_viscosity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Kinematic viscosity nu(T, P, X) [m^2/s]."
    );

    m.def(
        "prandtl",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return prandtl(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Prandtl number Pr(T, P, X) [-]."
    );

    m.def(
        "thermal_diffusivity",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return thermal_diffusivity(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Thermal diffusivity alpha(T, P, X) [m^2/s]."
    );

    m.def(
        "reynolds",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double V,
           double L)
        {
            auto X = to_vec(X_arr);
            return reynolds(T, P, X, V, L);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("V"),
        py::arg("L"),
        "Reynolds number Re = rho*V*L/mu [-]."
    );

    m.def(
        "peclet",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double V,
           double L)
        {
            auto X = to_vec(X_arr);
            return peclet(T, P, X, V, L);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("V"),
        py::arg("L"),
        "Peclet number Pe = V*L/alpha (thermal) [-]."
    );

    // Humid air utilities
    m.def(
        "standard_dry_air_composition",
        &standard_dry_air_composition,
        "Standard dry air mole-fraction composition over the thermo species set."
    );

    m.def(
        "humid_air_composition",
        &humid_air_composition,
        py::arg("T"),
        py::arg("P"),
        py::arg("RH"),
        "Humid air mole-fraction composition for temperature T [K], pressure P [Pa], and relative humidity RH (0-1)."
    );

    m.def(
        "dewpoint",
        &dewpoint,
        py::arg("T"),
        py::arg("P"),
        py::arg("RH"),
        "Dewpoint temperature [K] for ambient T [K], pressure P [Pa], and RH (0-1)."
    );

    m.def(
        "relative_humidity_from_dewpoint",
        &relative_humidity_from_dewpoint,
        py::arg("T"),
        py::arg("Tdp"),
        py::arg("P"),
        "Relative humidity (0-1) from dry-bulb T [K], dewpoint Tdp [K], and pressure P [Pa]."
    );

    // Composition conversion helpers
    m.def(
        "mole_to_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return mole_to_mass(X);
        },
        py::arg("X"),
        "Convert mole fractions X to mass fractions Y.")
    ;

    m.def(
        "mass_to_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_arr)
        {
            auto Y = to_vec(Y_arr);
            return mass_to_mole(Y);
        },
        py::arg("Y"),
        "Convert mass fractions Y to mole fractions X.")
    ;

    // Fraction utilities
    m.def(
        "normalize_fractions",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> arr)
        {
            auto v = to_vec(arr);
            return normalize_fractions(v);
        },
        py::arg("fractions"),
        "Normalize a vector of fractions so it sums to 1.0 (or return all zeros if input is all zeros)."
    );

    m.def(
        "convert_to_dry_fractions",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> arr)
        {
            auto v = to_vec(arr);
            return convert_to_dry_fractions(v);
        },
        py::arg("mole_fractions"),
        "Convert mole fractions to dry fractions by removing water vapor and renormalizing."
    );

    // Equivalence ratio (mole basis)
    m.def(
        "equivalence_ratio_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_mix  = to_vec(X_mix_arr);
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return equivalence_ratio_mole(X_mix, X_fuel, X_ox);
        },
        py::arg("X_mix"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Equivalence ratio phi (mole basis) for mixture X_mix formed from fuel X_fuel and oxidizer X_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mole",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return set_equivalence_ratio_mole(phi, X_fuel, X_ox);
        },
        py::arg("phi"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Construct mole-fraction mixture X_mix with target phi from fuel X_fuel and oxidizer X_ox.")
    ;

    // Equivalence ratio (mass basis)
    m.def(
        "equivalence_ratio_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_mix  = to_vec(Y_mix_arr);
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox);
        },
        py::arg("Y_mix"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Equivalence ratio phi (mass basis) for mixture Y_mix formed from fuel Y_fuel and oxidizer Y_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mass",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return set_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
        },
        py::arg("phi"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Construct mass-fraction mixture Y_mix with target phi from fuel Y_fuel and oxidizer Y_ox.")
    ;

    // Bilger-based helpers (mass basis)
    m.def(
        "bilger_stoich_mixture_fraction_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return bilger_stoich_mixture_fraction_mass(Y_fuel, Y_ox);
        },
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Stoichiometric Bilger mixture fraction Z_st (mass basis) for fuel Y_fuel and oxidizer Y_ox.");

    m.def(
        "bilger_Z_from_equivalence_ratio_mass",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return bilger_Z_from_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
        },
        py::arg("phi"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Convert mass-basis equivalence ratio phi to Bilger mixture fraction Z for fuel/oxidizer streams.");

    m.def(
        "equivalence_ratio_from_bilger_Z_mass",
        [](double Z,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return equivalence_ratio_from_bilger_Z_mass(Z, Y_fuel, Y_ox);
        },
        py::arg("Z"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Convert Bilger mixture fraction Z to mass-basis equivalence ratio phi for fuel/oxidizer streams.");

    // Inverse temperature calculations
    m.def(
        "calc_T_from_h",
        [](double h_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_h(h_target, X, T_guess, tol, max_iter);
        },
        py::arg("h_target"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target enthalpy h_target and composition X."
    );

    m.def(
        "calc_T_from_s",
        [](double s_target,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_s(s_target, P, X, T_guess, tol, max_iter);
        },
        py::arg("s_target"),
        py::arg("P"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target entropy s_target, pressure P, and composition X."
    );

    m.def(
        "calc_T_from_cp",
        [](double cp_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double T_guess,
           double tol,
           std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return calc_T_from_cp(cp_target, X, T_guess, tol, max_iter);
        },
        py::arg("cp_target"),
        py::arg("X"),
        py::arg("T_guess") = 300.0,
        py::arg("tol") = 1.0e-6,
        py::arg("max_iter") = 50,
        "Solve for temperature T given target Cp cp_target and composition X."
    );

    // Oxygen requirement helpers (scalar fuel index and mixtures)
    m.def(
        "oxygen_required_per_mol_fuel",
        &oxygen_required_per_mol_fuel,
        py::arg("fuel_index"),
        "Moles O2 required per mole of pure fuel species (by index)."
    );

    m.def(
        "oxygen_required_per_kg_fuel",
        &oxygen_required_per_kg_fuel,
        py::arg("fuel_index"),
        "Mass of O2 [kg] required per kg of pure fuel species (by index)."
    );

    m.def(
        "oxygen_required_per_mol_mixture",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return oxygen_required_per_mol_mixture(X);
        },
        py::arg("X"),
        "Moles O2 required per mole of fuel mixture with mole fractions X."
    );

    m.def(
        "oxygen_required_per_kg_mixture",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return oxygen_required_per_kg_mixture(X);
        },
        py::arg("X"),
        "Mass of O2 [kg] required per kg of fuel mixture with mole fractions X."
    );

    // Complete combustion helper
    m.def(
        "complete_combustion_to_CO2_H2O",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return complete_combustion_to_CO2_H2O(X);
        },
        py::arg("X"),
        "Complete combustion of mixture X to CO2 and H2O (returns product mole fractions)."
    );

    m.def(
        "complete_combustion_to_CO2_H2O_with_fraction",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            double f = 0.0;
            auto X_out = complete_combustion_to_CO2_H2O(X, f);
            return py::make_tuple(X_out, f);
        },
        py::arg("X"),
        "Complete combustion of mixture X to CO2 and H2O, returning (product mole fractions, fuel burn fraction f)."
    );

    m.def(
        "adiabatic_T_wgs",
        [](double T_in,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_in_arr,
           double P)
        {
            auto X_in = to_vec(X_in_arr);
            if (X_in.empty())
                throw std::runtime_error("X_in must be non-empty");

            State in;
            in.T = T_in;
            in.P = P;
            in.X = X_in;

            State out = wgs_equilibrium_adiabatic(in);
            return out.T;
        },
        py::arg("T_in"),
        py::arg("X_in"),
        py::arg("P") = 101325.0,
        "Adiabatic flame temperature with WGS limited equilibrium.\n\n"
        "T_in : inlet temperature [K]\n"
        "X_in : inlet mole fractions (1D NumPy array)\n"
        "P    : pressure [Pa] (default 101325)"
    );

    // State struct binding - Pythonic property-based API
    py::class_<State>(m, "State")
        .def(py::init<>())
        // Mutable state variables (read/write properties)
        .def_readwrite("T", &State::T, "Temperature [K]")
        .def_readwrite("P", &State::P, "Pressure [Pa]")
        .def_readwrite("X", &State::X, "Mole fractions [-]")
        // Computed thermodynamic properties (read-only)
        .def_property_readonly("mw", &State::mw, "Molecular weight [g/mol]")
        .def_property_readonly("cp", &State::cp, "Specific heat at constant pressure [J/(mol·K)]")
        .def_property_readonly("h", &State::h, "Specific enthalpy [J/mol]")
        .def_property_readonly("s", &State::s, "Specific entropy [J/(mol·K)]")
        .def_property_readonly("cv", &State::cv, "Specific heat at constant volume [J/(mol·K)]")
        .def_property_readonly("u", &State::u, "Specific internal energy [J/mol]")
        .def_property_readonly("rho", &State::rho, "Density [kg/m³]")
        .def_property_readonly("R", &State::R, "Specific gas constant [J/(mol·K)]")
        .def_property_readonly("gamma", &State::gamma, "Isentropic expansion coefficient [-]")
        .def_property_readonly("a", &State::a, "Speed of sound [m/s]")
        // Transport properties (read-only)
        .def_property_readonly("mu", &State::mu, "Dynamic viscosity [Pa·s]")
        .def_property_readonly("k", &State::k, "Thermal conductivity [W/(m·K)]")
        .def_property_readonly("nu", &State::nu, "Kinematic viscosity [m²/s]")
        .def_property_readonly("Pr", &State::Pr, "Prandtl number [-]")
        .def_property_readonly("alpha", &State::alpha, "Thermal diffusivity [m²/s]")
        // Fluent setters (return self for chaining)
        .def("set_T", &State::set_T, py::arg("T"), "Set temperature [K], returns self")
        .def("set_P", &State::set_P, py::arg("P"), "Set pressure [Pa], returns self")
        .def("set_X", [](State& s, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) -> State& {
            s.X = to_vec(X_arr);
            return s;
        }, py::arg("X"), "Set mole fractions [-], returns self");

    // Stream struct binding - Pythonic property-based API
    py::class_<Stream>(m, "Stream")
        .def(py::init<>())
        .def_readwrite("state", &Stream::state, "Thermodynamic state")
        // Mutable state (read/write properties)
        .def_property("T",
            [](const Stream& s) { return s.T(); },
            [](Stream& s, double T) { s.set_T(T); },
            "Temperature [K]")
        .def_property("P",
            [](const Stream& s) { return s.P(); },
            [](Stream& s, double P) { s.set_P(P); },
            "Pressure [Pa]")
        .def_property("X",
            [](const Stream& s) { return s.X(); },
            [](Stream& s, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
                s.state.X = to_vec(X_arr);
            },
            "Mole fractions [-]")
        .def_readwrite("mdot", &Stream::mdot, "Mass flow rate [kg/s]")
        // Computed properties (read-only)
        .def_property_readonly("mw", &Stream::mw, "Molecular weight [g/mol]")
        .def_property_readonly("cp", &Stream::cp, "Specific heat at constant pressure [J/(mol·K)]")
        .def_property_readonly("h", &Stream::h, "Specific enthalpy [J/mol]")
        .def_property_readonly("s", &Stream::s, "Specific entropy [J/(mol·K)]")
        .def_property_readonly("rho", &Stream::rho, "Density [kg/m³]")
        // Fluent setters (return self for chaining)
        .def("set_T", &Stream::set_T, py::arg("T"), "Set temperature [K], returns self")
        .def("set_P", &Stream::set_P, py::arg("P"), "Set pressure [Pa], returns self")
        .def("set_X", [](Stream& s, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) -> Stream& {
            s.state.X = to_vec(X_arr);
            return s;
        }, py::arg("X"), "Set mole fractions [-], returns self")
        .def("set_mdot", &Stream::set_mdot, py::arg("mdot"), "Set mass flow rate [kg/s], returns self");

    // Stream mixing function
    m.def(
        "mix",
        [](const std::vector<Stream>& streams, double P_out) {
            return mix(streams, P_out);
        },
        py::arg("streams"),
        py::arg("P_out") = -1.0,
        "Mix multiple streams, solving mass and enthalpy balance.\n\n"
        "streams : list of Stream objects\n"
        "P_out   : output pressure [Pa] (default: minimum inlet pressure)\n\n"
        "Returns mixed Stream."
    );

    // -------------------------------------------------------------
    // Inverse solvers for fuel/oxidizer streams (complete combustion only)
    // -------------------------------------------------------------

    // --- Find fuel stream (oxidizer mdot fixed) ---

    m.def(
        "set_fuel_stream_for_Tad",
        [](double T_ad_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter, bool lean, double phi_max) {
            return set_fuel_stream_for_Tad(T_ad_target, fuel, oxidizer, tol, max_iter, lean, phi_max);
        },
        py::arg("T_ad_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1.0,
        py::arg("max_iter") = 100,
        py::arg("lean") = true,
        py::arg("phi_max") = 10.0,
        "Find fuel mass flow rate to achieve target adiabatic flame temperature.\n\n"
        "Uses complete combustion to CO2/H2O (no WGS).\n\n"
        "T_ad_target : target adiabatic flame temperature [K] (must be > oxidizer T)\n"
        "fuel        : fuel Stream (T, P, X set; mdot ignored)\n"
        "oxidizer    : oxidizer Stream (with mdot set)\n"
        "tol         : temperature tolerance [K] (default: 1.0)\n"
        "max_iter    : maximum iterations (default: 100)\n"
        "lean        : if True (default), search on lean side; if False, search on rich side\n"
        "phi_max     : maximum equivalence ratio for rich side search (default: 10.0)\n\n"
        "Returns fuel Stream with mdot set to achieve target T_ad."
    );

    m.def(
        "set_fuel_stream_for_O2",
        [](double X_O2_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_fuel_stream_for_O2(X_O2_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_O2_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find fuel mass flow rate to achieve target O2 mole fraction in burned products (wet basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only."
    );

    m.def(
        "set_fuel_stream_for_O2_dry",
        [](double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_fuel_stream_for_O2_dry(X_O2_dry_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_O2_dry_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find fuel mass flow rate to achieve target O2 mole fraction in burned products (dry basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
        "Dry basis: water vapor removed from products before computing mole fraction."
    );

    m.def(
        "set_fuel_stream_for_CO2",
        [](double X_CO2_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_fuel_stream_for_CO2(X_CO2_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_CO2_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find fuel mass flow rate to achieve target CO2 mole fraction in burned products (wet basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only."
    );

    m.def(
        "set_fuel_stream_for_CO2_dry",
        [](double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_fuel_stream_for_CO2_dry(X_CO2_dry_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_CO2_dry_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find fuel mass flow rate to achieve target CO2 mole fraction in burned products (dry basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
        "Dry basis: water vapor removed from products before computing mole fraction."
    );

    // --- Find oxidizer stream (fuel mdot fixed) ---

    m.def(
        "set_oxidizer_stream_for_Tad",
        [](double T_ad_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter, bool lean, double phi_max) {
            return set_oxidizer_stream_for_Tad(T_ad_target, fuel, oxidizer, tol, max_iter, lean, phi_max);
        },
        py::arg("T_ad_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1.0,
        py::arg("max_iter") = 100,
        py::arg("lean") = true,
        py::arg("phi_max") = 10.0,
        "Find oxidizer mass flow rate to achieve target adiabatic flame temperature.\n\n"
        "Uses complete combustion to CO2/H2O (no WGS).\n\n"
        "T_ad_target : target adiabatic flame temperature [K] (must be > fuel T)\n"
        "fuel        : fuel Stream (with mdot set)\n"
        "oxidizer    : oxidizer Stream (T, P, X set; mdot ignored)\n"
        "tol         : temperature tolerance [K] (default: 1.0)\n"
        "max_iter    : maximum iterations (default: 100)\n"
        "lean        : if True (default), search on lean side; if False, search on rich side\n"
        "phi_max     : maximum equivalence ratio for search range (default: 10.0)\n\n"
        "Returns oxidizer Stream with mdot set to achieve target T_ad."
    );

    m.def(
        "set_oxidizer_stream_for_O2",
        [](double X_O2_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_oxidizer_stream_for_O2(X_O2_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_O2_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find oxidizer mass flow rate to achieve target O2 mole fraction in burned products (wet basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only."
    );

    m.def(
        "set_oxidizer_stream_for_O2_dry",
        [](double X_O2_dry_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_oxidizer_stream_for_O2_dry(X_O2_dry_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_O2_dry_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find oxidizer mass flow rate to achieve target O2 mole fraction in burned products (dry basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
        "Dry basis: water vapor removed from products before computing mole fraction."
    );

    m.def(
        "set_oxidizer_stream_for_CO2",
        [](double X_CO2_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_oxidizer_stream_for_CO2(X_CO2_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_CO2_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find oxidizer mass flow rate to achieve target CO2 mole fraction in burned products (wet basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only."
    );

    m.def(
        "set_oxidizer_stream_for_CO2_dry",
        [](double X_CO2_dry_target, const Stream& fuel, const Stream& oxidizer, double tol, std::size_t max_iter) {
            return set_oxidizer_stream_for_CO2_dry(X_CO2_dry_target, fuel, oxidizer, tol, max_iter);
        },
        py::arg("X_CO2_dry_target"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Find oxidizer mass flow rate to achieve target CO2 mole fraction in burned products (dry basis).\n\n"
        "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
        "Dry basis: water vapor removed from products before computing mole fraction."
    );

    // State-based combustion functions
    m.def(
        "complete_combustion",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return complete_combustion(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic complete combustion to CO2 and H2O.\n\n"
        "Returns State with adiabatic flame temperature and burned composition."
    );

    m.def(
        "complete_combustion_isothermal",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return complete_combustion_isothermal(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal complete combustion to CO2 and H2O.\n\n"
        "Returns State with same temperature and burned composition."
    );

    // State-based WGS equilibrium functions
    m.def(
        "wgs_equilibrium",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return wgs_equilibrium(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
        "Returns State with equilibrium composition at input temperature."
    );

    m.def(
        "wgs_equilibrium_adiabatic",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return wgs_equilibrium_adiabatic(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
        "Returns State with equilibrium temperature and composition."
    );

    // SMR+WGS equilibrium (CH4 only)
    m.def(
        "smr_wgs_equilibrium",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return smr_wgs_equilibrium(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal SMR+WGS equilibrium for CH4.\n\n"
        "Steam Methane Reforming: CH4 + H2O <-> CO + 3H2\n"
        "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
        "Returns State with equilibrium composition at input temperature."
    );

    m.def(
        "smr_wgs_equilibrium_adiabatic",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return smr_wgs_equilibrium_adiabatic(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic SMR+WGS equilibrium for CH4.\n\n"
        "Steam Methane Reforming: CH4 + H2O <-> CO + 3H2\n"
        "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
        "Returns State with equilibrium temperature and composition.\n"
        "Temperature decreases due to endothermic SMR reaction."
    );

    // General reforming + WGS equilibrium (all hydrocarbons)
    m.def(
        "reforming_equilibrium",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return reforming_equilibrium(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Isothermal steam reforming + WGS equilibrium for all hydrocarbons.\n\n"
        "General reforming: CnHm + n*H2O <-> n*CO + (n + m/2)*H2\n"
        "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
        "Handles all hydrocarbons: CH4, C2H6, C3H8, iC4H10, nC5H12, etc.\n"
        "Returns State with equilibrium composition at input temperature."
    );

    m.def(
        "reforming_equilibrium_adiabatic",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return reforming_equilibrium_adiabatic(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Adiabatic steam reforming + WGS equilibrium for all hydrocarbons.\n\n"
        "General reforming: CnHm + n*H2O <-> n*CO + (n + m/2)*H2\n"
        "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
        "Handles all hydrocarbons: CH4, C2H6, C3H8, iC4H10, nC5H12, etc.\n"
        "Returns State with equilibrium temperature and composition.\n"
        "Temperature decreases due to endothermic reforming reactions."
    );

    // Combustion + Equilibrium (convenience function)
    m.def(
        "combustion_equilibrium",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, double P)
        {
            State in;
            in.T = T;
            in.P = P;
            in.X = to_vec(X_arr);
            return combustion_equilibrium(in);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P") = 101325.0,
        "Combustion + equilibrium in one step.\n\n"
        "Combines complete combustion with reforming + WGS equilibrium.\n"
        "Use this when starting from an unburned fuel+air mixture.\n\n"
        "Workflow:\n"
        "  1. Complete combustion (fuel + O2 -> CO2 + H2O)\n"
        "  2. Reforming + WGS equilibrium on products\n\n"
        "Returns State with equilibrium temperature and composition."
    );

    // -------------------------------------------------------------
    // Compressible flow
    // -------------------------------------------------------------

    // CompressibleFlowSolution struct
    py::class_<CompressibleFlowSolution>(m, "CompressibleFlowSolution")
        .def(py::init<>())
        .def_readwrite("stagnation", &CompressibleFlowSolution::stagnation, "Stagnation state")
        .def_readwrite("outlet", &CompressibleFlowSolution::outlet, "Outlet static state")
        .def_readwrite("v", &CompressibleFlowSolution::v, "Outlet velocity [m/s]")
        .def_readwrite("M", &CompressibleFlowSolution::M, "Outlet Mach number [-]")
        .def_readwrite("mdot", &CompressibleFlowSolution::mdot, "Mass flow rate [kg/s]")
        .def_readwrite("choked", &CompressibleFlowSolution::choked, "True if flow is choked");

    // NozzleStation struct
    py::class_<NozzleStation>(m, "NozzleStation")
        .def(py::init<>())
        .def_readwrite("x", &NozzleStation::x, "Axial position [m]")
        .def_readwrite("A", &NozzleStation::A, "Area [m²]")
        .def_readwrite("P", &NozzleStation::P, "Static pressure [Pa]")
        .def_readwrite("T", &NozzleStation::T, "Static temperature [K]")
        .def_readwrite("rho", &NozzleStation::rho, "Density [kg/m³]")
        .def_readwrite("u", &NozzleStation::u, "Velocity [m/s]")
        .def_readwrite("M", &NozzleStation::M, "Mach number [-]")
        .def_readwrite("h", &NozzleStation::h, "Specific enthalpy [J/kg]");

    // NozzleSolution struct
    py::class_<NozzleSolution>(m, "NozzleSolution")
        .def(py::init<>())
        .def_readwrite("inlet", &NozzleSolution::inlet, "Inlet state")
        .def_readwrite("outlet", &NozzleSolution::outlet, "Outlet state")
        .def_readwrite("mdot", &NozzleSolution::mdot, "Mass flow rate [kg/s]")
        .def_readwrite("h0", &NozzleSolution::h0, "Stagnation enthalpy [J/kg]")
        .def_readwrite("T0", &NozzleSolution::T0, "Stagnation temperature [K]")
        .def_readwrite("P0", &NozzleSolution::P0, "Stagnation pressure [Pa]")
        .def_readwrite("choked", &NozzleSolution::choked, "True if throat is sonic")
        .def_readwrite("x_throat", &NozzleSolution::x_throat, "Throat position [m]")
        .def_readwrite("A_throat", &NozzleSolution::A_throat, "Throat area [m²]")
        .def_readwrite("profile", &NozzleSolution::profile, "Axial profile");

    // ThrustResult struct
    py::class_<ThrustResult>(m, "ThrustResult")
        .def(py::init<>())
        .def_readwrite("thrust", &ThrustResult::thrust, "Thrust force [N]")
        .def_readwrite("specific_impulse", &ThrustResult::specific_impulse, "Specific impulse Isp [s]")
        .def_readwrite("thrust_coefficient", &ThrustResult::thrust_coefficient, "Thrust coefficient C_F [-]")
        .def_readwrite("mdot", &ThrustResult::mdot, "Mass flow rate [kg/s]")
        .def_readwrite("u_exit", &ThrustResult::u_exit, "Exit velocity [m/s]")
        .def_readwrite("P_exit", &ThrustResult::P_exit, "Exit pressure [Pa]");

    // FannoStation struct
    py::class_<FannoStation>(m, "FannoStation")
        .def(py::init<>())
        .def_readwrite("x", &FannoStation::x, "Position [m]")
        .def_readwrite("P", &FannoStation::P, "Static pressure [Pa]")
        .def_readwrite("T", &FannoStation::T, "Static temperature [K]")
        .def_readwrite("rho", &FannoStation::rho, "Density [kg/m³]")
        .def_readwrite("u", &FannoStation::u, "Velocity [m/s]")
        .def_readwrite("M", &FannoStation::M, "Mach number [-]")
        .def_readwrite("h", &FannoStation::h, "Specific enthalpy [J/kg]")
        .def_readwrite("s", &FannoStation::s, "Specific entropy [J/(kg·K)]");

    // FannoSolution struct
    py::class_<FannoSolution>(m, "FannoSolution")
        .def(py::init<>())
        .def_readwrite("inlet", &FannoSolution::inlet, "Inlet state")
        .def_readwrite("outlet", &FannoSolution::outlet, "Outlet state")
        .def_readwrite("mdot", &FannoSolution::mdot, "Mass flow rate [kg/s]")
        .def_readwrite("h0", &FannoSolution::h0, "Stagnation enthalpy [J/kg]")
        .def_readwrite("L", &FannoSolution::L, "Pipe length [m]")
        .def_readwrite("D", &FannoSolution::D, "Pipe diameter [m]")
        .def_readwrite("f", &FannoSolution::f, "Darcy friction factor [-]")
        .def_readwrite("choked", &FannoSolution::choked, "True if flow reached M=1")
        .def_readwrite("L_choke", &FannoSolution::L_choke, "Length to choking [m]")
        .def_readwrite("profile", &FannoSolution::profile, "Axial profile");

    // Isentropic nozzle flow
    m.def(
        "nozzle_flow",
        [](double T0, double P0, double P_back, double A_eff,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return nozzle_flow(T0, P0, P_back, A_eff, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P_back"),
        py::arg("A_eff"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Isentropic nozzle flow.\n\n"
        "T0     : stagnation temperature [K]\n"
        "P0     : stagnation pressure [Pa]\n"
        "P_back : back pressure [Pa]\n"
        "A_eff  : effective flow area [m²]\n"
        "X      : mole fractions\n\n"
        "Returns CompressibleFlowSolution."
    );

    // Inverse solvers
    m.def(
        "solve_A_eff_from_mdot",
        [](double T0, double P0, double P_back, double mdot_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return solve_A_eff_from_mdot(T0, P0, P_back, mdot_target, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P_back"),
        py::arg("mdot_target"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Find effective area for given mass flow rate."
    );

    m.def(
        "solve_P_back_from_mdot",
        [](double T0, double P0, double A_eff, double mdot_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return solve_P_back_from_mdot(T0, P0, A_eff, mdot_target, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("A_eff"),
        py::arg("mdot_target"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Find back pressure for given mass flow rate."
    );

    m.def(
        "solve_P0_from_mdot",
        [](double T0, double P_back, double A_eff, double mdot_target,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return solve_P0_from_mdot(T0, P_back, A_eff, mdot_target, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P_back"),
        py::arg("A_eff"),
        py::arg("mdot_target"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Find stagnation pressure for given mass flow rate."
    );

    // Utility functions
    m.def(
        "critical_pressure_ratio",
        [](double T0, double P0,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return critical_pressure_ratio(T0, P0, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Critical (sonic) pressure ratio P*/P0."
    );

    m.def(
        "mach_from_pressure_ratio",
        [](double T0, double P0, double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return mach_from_pressure_ratio(T0, P0, P, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Mach number from pressure ratio P/P0 (isentropic)."
    );

    m.def(
        "mass_flux_isentropic",
        [](double T0, double P0, double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return mass_flux_isentropic(T0, P0, P, X, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P"),
        py::arg("X"),
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Mass flux G = rho*v [kg/(m²·s)] at given pressure (isentropic)."
    );

    // Quasi-1D nozzle with C-D geometry
    m.def(
        "nozzle_cd",
        [](double T0, double P0, double P_exit,
           double A_inlet, double A_throat, double A_exit,
           double x_throat, double x_exit,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           std::size_t n_stations, double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return nozzle_cd(T0, P0, P_exit, A_inlet, A_throat, A_exit,
                             x_throat, x_exit, X, n_stations, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P_exit"),
        py::arg("A_inlet"),
        py::arg("A_throat"),
        py::arg("A_exit"),
        py::arg("x_throat"),
        py::arg("x_exit"),
        py::arg("X"),
        py::arg("n_stations") = 100,
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Converging-diverging nozzle flow.\n\n"
        "Uses smooth cosine area profile between inlet, throat, and exit.\n\n"
        "Returns NozzleSolution with axial profile."
    );

    // Fanno flow
    m.def(
        "fanno_pipe",
        [](double T_in, double P_in, double u_in,
           double L, double D, double f,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           std::size_t n_steps, bool store_profile)
        {
            auto X = to_vec(X_arr);
            return fanno_pipe(T_in, P_in, u_in, L, D, f, X, n_steps, store_profile);
        },
        py::arg("T_in"),
        py::arg("P_in"),
        py::arg("u_in"),
        py::arg("L"),
        py::arg("D"),
        py::arg("f"),
        py::arg("X"),
        py::arg("n_steps") = 100,
        py::arg("store_profile") = false,
        "Fanno flow (adiabatic pipe with friction).\n\n"
        "T_in : inlet temperature [K]\n"
        "P_in : inlet pressure [Pa]\n"
        "u_in : inlet velocity [m/s]\n"
        "L    : pipe length [m]\n"
        "D    : pipe diameter [m]\n"
        "f    : Darcy friction factor [-]\n\n"
        "Returns FannoSolution."
    );

    m.def(
        "fanno_max_length",
        [](double T_in, double P_in, double u_in,
           double D, double f,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return fanno_max_length(T_in, P_in, u_in, D, f, X, tol, max_iter);
        },
        py::arg("T_in"),
        py::arg("P_in"),
        py::arg("u_in"),
        py::arg("D"),
        py::arg("f"),
        py::arg("X"),
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 100,
        "Maximum pipe length before choking (L*)."
    );

    // -------------------------------------------------------------
    // Friction factor correlations
    // -------------------------------------------------------------

    m.def(
        "friction_haaland",
        &friction_haaland,
        py::arg("Re"),
        py::arg("e_D"),
        "Haaland friction factor correlation (explicit, ~2-3% accuracy).\n\n"
        "Re  : Reynolds number [-]\n"
        "e_D : relative roughness ε/D [-]"
    );

    m.def(
        "friction_serghides",
        &friction_serghides,
        py::arg("Re"),
        py::arg("e_D"),
        "Serghides friction factor correlation (explicit, <0.3% accuracy).\n\n"
        "Re  : Reynolds number [-]\n"
        "e_D : relative roughness ε/D [-]"
    );

    m.def(
        "friction_colebrook",
        &friction_colebrook,
        py::arg("Re"),
        py::arg("e_D"),
        py::arg("tol") = 1e-10,
        py::arg("max_iter") = 20,
        "Colebrook-White friction factor (implicit, reference standard).\n\n"
        "Re  : Reynolds number [-]\n"
        "e_D : relative roughness ε/D [-]"
    );

    m.def(
        "friction_petukhov",
        &friction_petukhov,
        py::arg("Re"),
        "Petukhov friction factor for smooth tubes (Re > 3000).\n\n"
        "Re : Reynolds number [-]\n\n"
        "Returns Darcy friction factor f [-]."
    );

    // -------------------------------------------------------------
    // Heat transfer correlations
    // -------------------------------------------------------------

    m.def(
        "nusselt_dittus_boelter",
        &nusselt_dittus_boelter,
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("heating") = true,
        "Dittus-Boelter correlation for turbulent pipe flow.\n\n"
        "Nu = 0.023 * Re^0.8 * Pr^n\n"
        "n = 0.4 (heating) or 0.3 (cooling)\n\n"
        "Valid: Re > 10000, 0.6 < Pr < 160, L/D > 10"
    );

    m.def(
        "nusselt_gnielinski",
        [](double Re, double Pr, double f) {
            if (f < 0) {
                return nusselt_gnielinski(Re, Pr);
            }
            return nusselt_gnielinski(Re, Pr, f);
        },
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("f") = -1.0,
        "Gnielinski correlation for transitional/turbulent pipe flow.\n\n"
        "Valid: 3000 < Re < 5e6, 0.5 < Pr < 2000\n\n"
        "f : friction factor (if < 0, uses Petukhov correlation)"
    );

    m.def(
        "nusselt_sieder_tate",
        &nusselt_sieder_tate,
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("mu_ratio"),
        "Sieder-Tate correlation with viscosity correction.\n\n"
        "Nu = 0.027 * Re^0.8 * Pr^(1/3) * (μ/μ_w)^0.14\n\n"
        "mu_ratio : μ_bulk / μ_wall"
    );

    m.def(
        "nusselt_petukhov",
        [](double Re, double Pr, double f) {
            if (f < 0) {
                return nusselt_petukhov(Re, Pr);
            }
            return nusselt_petukhov(Re, Pr, f);
        },
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("f") = -1.0,
        "Petukhov correlation for turbulent pipe flow.\n\n"
        "Valid: 1e4 < Re < 5e6, 0.5 < Pr < 2000\n\n"
        "f : friction factor (if < 0, uses Petukhov correlation)"
    );

    m.def(
        "htc_from_nusselt",
        &htc_from_nusselt,
        py::arg("Nu"),
        py::arg("k"),
        py::arg("L"),
        "Heat transfer coefficient from Nusselt number.\n\n"
        "h = Nu * k / L  [W/(m²·K)]\n\n"
        "Nu : Nusselt number [-]\n"
        "k  : thermal conductivity [W/(m·K)]\n"
        "L  : characteristic length [m]"
    );

    m.def(
        "overall_htc",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> h_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto h_values = to_vec(h_arr);
            auto t_over_k = to_vec(tk_arr);
            return overall_htc(h_values, t_over_k);
        },
        py::arg("h_values"),
        py::arg("t_over_k") = py::array_t<double>(),
        "Overall HTC from thermal resistance network.\n\n"
        "1/U = Σ(1/h_i) + Σ(t_j/k_j)\n\n"
        "h_values : convective HTCs [W/(m²·K)]\n"
        "t_over_k : thickness/conductivity ratios [m²·K/W]"
    );

    m.def(
        "overall_htc_wall",
        py::overload_cast<double, double, double, double>(&overall_htc_wall),
        py::arg("h_inner"),
        py::arg("h_outer"),
        py::arg("t_wall"),
        py::arg("k_wall"),
        "Overall HTC for single-layer wall.\n\n"
        "1/U = 1/h_inner + t_wall/k_wall + 1/h_outer"
    );

    m.def(
        "overall_htc_wall_multilayer",
        [](double h_inner, double h_outer,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr,
           double R_fouling)
        {
            auto layers = to_vec(tk_arr);
            if (R_fouling >= 0) {
                return overall_htc_wall(h_inner, h_outer, layers, R_fouling);
            }
            return overall_htc_wall(h_inner, h_outer, layers);
        },
        py::arg("h_inner"),
        py::arg("h_outer"),
        py::arg("t_over_k_layers"),
        py::arg("R_fouling") = -1.0,
        "Overall HTC for multi-layer wall.\n\n"
        "1/U = 1/h_inner + Σ(t_i/k_i) + 1/h_outer [+ R_fouling]\n\n"
        "Example: steel + insulation\n"
        "  overall_htc_wall_multilayer(500, 10, [0.003/50, 0.05/0.04])"
    );

    m.def(
        "thermal_resistance",
        &thermal_resistance,
        py::arg("h"),
        py::arg("A"),
        "Thermal resistance for convection.\n\n"
        "R = 1/(h*A) [K/W]"
    );

    m.def(
        "thermal_resistance_wall",
        &thermal_resistance_wall,
        py::arg("thickness"),
        py::arg("k"),
        py::arg("A"),
        "Thermal resistance for conduction through wall.\n\n"
        "R = t/(k*A) [K/W]"
    );

    m.def(
        "lmtd",
        &lmtd,
        py::arg("dT1"),
        py::arg("dT2"),
        "Log Mean Temperature Difference.\n\n"
        "LMTD = (ΔT1 - ΔT2) / ln(ΔT1/ΔT2)\n\n"
        "Handles equal ΔT gracefully (returns arithmetic mean)."
    );

    m.def(
        "lmtd_counterflow",
        &lmtd_counterflow,
        py::arg("T_hot_in"),
        py::arg("T_hot_out"),
        py::arg("T_cold_in"),
        py::arg("T_cold_out"),
        "LMTD for counter-flow heat exchanger.\n\n"
        "ΔT1 = T_hot_in - T_cold_out\n"
        "ΔT2 = T_hot_out - T_cold_in"
    );

    m.def(
        "lmtd_parallelflow",
        &lmtd_parallelflow,
        py::arg("T_hot_in"),
        py::arg("T_hot_out"),
        py::arg("T_cold_in"),
        py::arg("T_cold_out"),
        "LMTD for parallel-flow heat exchanger.\n\n"
        "ΔT1 = T_hot_in - T_cold_in\n"
        "ΔT2 = T_hot_out - T_cold_out"
    );

    m.def(
        "heat_rate",
        &heat_rate,
        py::arg("U"),
        py::arg("A"),
        py::arg("dT"),
        "Heat transfer rate Q = U * A * ΔT [W]"
    );

    m.def(
        "heat_flux",
        &heat_flux,
        py::arg("U"),
        py::arg("dT"),
        "Heat flux q = U * ΔT [W/m²]"
    );

    m.def(
        "heat_transfer_area",
        &heat_transfer_area,
        py::arg("Q"),
        py::arg("U"),
        py::arg("dT"),
        "Required area A = Q / (U * ΔT) [m²]"
    );

    m.def(
        "heat_transfer_dT",
        &heat_transfer_dT,
        py::arg("Q"),
        py::arg("U"),
        py::arg("A"),
        "Temperature difference ΔT = Q / (U * A) [K]"
    );

    m.def(
        "wall_temperature_profile",
        [](double T_hot, double T_cold, double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            double q;
            auto temps = wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k, q);
            return py::make_tuple(temps, q);
        },
        py::arg("T_hot"),
        py::arg("T_cold"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        "Temperature profile through multi-layer wall.\n\n"
        "Returns (temps, q) where:\n"
        "  temps : interface temperatures [K]\n"
        "  q     : heat flux [W/m²]"
    );

    m.def(
        "ntu",
        &ntu,
        py::arg("U"),
        py::arg("A"),
        py::arg("C_min"),
        "Number of Transfer Units: NTU = U * A / C_min"
    );

    m.def(
        "capacity_ratio",
        &capacity_ratio,
        py::arg("C_min"),
        py::arg("C_max"),
        "Heat capacity ratio: C_r = C_min / C_max"
    );

    m.def(
        "effectiveness_counterflow",
        &effectiveness_counterflow,
        py::arg("NTU"),
        py::arg("C_r"),
        "Effectiveness for counter-flow heat exchanger.\n\n"
        "ε = (1 - exp(-NTU*(1-C_r))) / (1 - C_r*exp(-NTU*(1-C_r)))"
    );

    m.def(
        "effectiveness_parallelflow",
        &effectiveness_parallelflow,
        py::arg("NTU"),
        py::arg("C_r"),
        "Effectiveness for parallel-flow heat exchanger.\n\n"
        "ε = (1 - exp(-NTU*(1+C_r))) / (1 + C_r)"
    );

    m.def(
        "heat_rate_from_effectiveness",
        &heat_rate_from_effectiveness,
        py::arg("epsilon"),
        py::arg("C_min"),
        py::arg("T_hot_in"),
        py::arg("T_cold_in"),
        "Heat rate from effectiveness: Q = ε * C_min * (T_hot_in - T_cold_in) [W]"
    );

    m.def(
        "heat_flux_from_T_at_edge",
        [](double T_measured, std::size_t edge_idx,
           double T_hot, double T_cold,
           double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            return heat_flux_from_T_at_edge(T_measured, edge_idx, T_hot, T_cold,
                                            h_hot, h_cold, t_over_k);
        },
        py::arg("T_measured"),
        py::arg("edge_idx"),
        py::arg("T_hot"),
        py::arg("T_cold"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        "Heat flux from temperature measured at wall edge.\n\n"
        "Edge indexing: 0 = hot surface, N = cold surface.\n"
        "Returns heat flux q [W/m²]."
    );

    m.def(
        "heat_flux_from_T_at_depth",
        [](double T_measured, double depth_from_hot,
           double T_hot, double T_cold,
           double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> t_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> k_arr)
        {
            auto thicknesses = to_vec(t_arr);
            auto conductivities = to_vec(k_arr);
            return heat_flux_from_T_at_depth(T_measured, depth_from_hot, T_hot, T_cold,
                                             h_hot, h_cold, thicknesses, conductivities);
        },
        py::arg("T_measured"),
        py::arg("depth_from_hot"),
        py::arg("T_hot"),
        py::arg("T_cold"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("thicknesses"),
        py::arg("conductivities"),
        "Heat flux from temperature measured at depth from hot surface.\n\n"
        "Useful when thermocouple is embedded slightly into wall.\n"
        "Returns heat flux q [W/m²]."
    );

    m.def(
        "bulk_T_from_edge_T_and_q",
        [](double T_measured, std::size_t edge_idx, double q,
           double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr,
           const std::string& solve_for)
        {
            auto t_over_k = to_vec(tk_arr);
            return bulk_T_from_edge_T_and_q(T_measured, edge_idx, q,
                                            h_hot, h_cold, t_over_k, solve_for);
        },
        py::arg("T_measured"),
        py::arg("edge_idx"),
        py::arg("q"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        py::arg("solve_for"),
        "Infer bulk fluid temperature from edge temperature and known heat flux.\n\n"
        "solve_for: 'hot' or 'cold'\n"
        "Returns inferred bulk temperature [K]."
    );

    m.def(
        "dT_edge_dT_hot",
        [](std::size_t edge_idx, double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            return dT_edge_dT_hot(edge_idx, h_hot, h_cold, t_over_k);
        },
        py::arg("edge_idx"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        "Sensitivity of edge temperature to hot-side bulk temperature.\n\n"
        "Returns ∂T_edge/∂T_hot [-] (0 to 1)."
    );

    m.def(
        "dT_edge_dT_cold",
        [](std::size_t edge_idx, double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            return dT_edge_dT_cold(edge_idx, h_hot, h_cold, t_over_k);
        },
        py::arg("edge_idx"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        "Sensitivity of edge temperature to cold-side bulk temperature.\n\n"
        "Returns ∂T_edge/∂T_cold [-] (0 to 1)."
    );

    m.def(
        "dT_edge_dT_bulk",
        [](std::size_t edge_idx, double h_hot, double h_cold,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            auto [dT_hot, dT_cold] = dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k);
            return py::make_tuple(dT_hot, dT_cold);
        },
        py::arg("edge_idx"),
        py::arg("h_hot"),
        py::arg("h_cold"),
        py::arg("t_over_k"),
        "Both sensitivities at once.\n\n"
        "Returns (∂T_edge/∂T_hot, ∂T_edge/∂T_cold)."
    );

    m.def(
        "dT_edge_dq",
        [](std::size_t edge_idx, double h_hot,
           py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr)
        {
            auto t_over_k = to_vec(tk_arr);
            return dT_edge_dq(edge_idx, h_hot, t_over_k);
        },
        py::arg("edge_idx"),
        py::arg("h_hot"),
        py::arg("t_over_k"),
        "Sensitivity of edge temperature to heat flux.\n\n"
        "Returns ∂T_edge/∂q [K·m²/W] (negative)."
    );

    // Nozzle thrust (from NozzleSolution)
    m.def(
        "nozzle_thrust",
        py::overload_cast<const NozzleSolution&, double>(&nozzle_thrust),
        py::arg("sol"),
        py::arg("P_amb"),
        "Calculate rocket nozzle thrust from nozzle solution.\n\n"
        "F = mdot * u_e + (P_e - P_amb) * A_e\n\n"
        "sol   : NozzleSolution from nozzle_cd()\n"
        "P_amb : Ambient pressure [Pa]\n\n"
        "Returns ThrustResult with thrust, Isp, and C_F."
    );

    // Nozzle thrust (convenience, from parameters)
    m.def(
        "nozzle_thrust_cd",
        [](double T0, double P0, double P_amb,
           double A_inlet, double A_throat, double A_exit,
           double x_throat, double x_exit,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           std::size_t n_stations, double tol, std::size_t max_iter)
        {
            auto X = to_vec(X_arr);
            return nozzle_thrust(T0, P0, P_amb, A_inlet, A_throat, A_exit,
                                 x_throat, x_exit, X, n_stations, tol, max_iter);
        },
        py::arg("T0"),
        py::arg("P0"),
        py::arg("P_amb"),
        py::arg("A_inlet"),
        py::arg("A_throat"),
        py::arg("A_exit"),
        py::arg("x_throat"),
        py::arg("x_exit"),
        py::arg("X"),
        py::arg("n_stations") = 100,
        py::arg("tol") = 1e-8,
        py::arg("max_iter") = 50,
        "Calculate rocket nozzle thrust directly from geometry.\n\n"
        "Convenience function that calls nozzle_cd() then nozzle_thrust().\n\n"
        "Returns ThrustResult with thrust, Isp, and C_F."
    );

    // =========================================================================
    // Incompressible flow
    // =========================================================================

    // Bernoulli equation
    m.def(
        "bernoulli_P2",
        &bernoulli_P2,
        py::arg("P1"),
        py::arg("v1"),
        py::arg("v2"),
        py::arg("rho"),
        py::arg("dz") = 0.0,
        py::arg("g") = 9.80665,
        "Bernoulli equation: solve for downstream pressure P2.\n\n"
        "P1 + ½ρv1² + ρgh1 = P2 + ½ρv2² + ρgh2\n\n"
        "P1  : upstream pressure [Pa]\n"
        "v1  : upstream velocity [m/s]\n"
        "v2  : downstream velocity [m/s]\n"
        "rho : fluid density [kg/m³]\n"
        "dz  : elevation change z2-z1 [m] (positive=up)\n"
        "g   : gravitational acceleration [m/s²]"
    );

    m.def(
        "bernoulli_v2",
        &bernoulli_v2,
        py::arg("P1"),
        py::arg("P2"),
        py::arg("v1"),
        py::arg("rho"),
        py::arg("dz") = 0.0,
        py::arg("g") = 9.80665,
        "Bernoulli equation: solve for downstream velocity v2.\n\n"
        "P1  : upstream pressure [Pa]\n"
        "P2  : downstream pressure [Pa]\n"
        "v1  : upstream velocity [m/s]\n"
        "rho : fluid density [kg/m³]\n"
        "dz  : elevation change z2-z1 [m]\n"
        "g   : gravitational acceleration [m/s²]"
    );

    // Orifice flow
    m.def(
        "orifice_mdot",
        static_cast<double(*)(double, double, double, double, double)>(&orifice_mdot),
        py::arg("P1"),
        py::arg("P2"),
        py::arg("A"),
        py::arg("Cd"),
        py::arg("rho"),
        "Orifice mass flow rate (incompressible).\n\n"
        "ṁ = Cd · A · √(2 · ρ · ΔP)\n\n"
        "P1  : upstream pressure [Pa]\n"
        "P2  : downstream pressure [Pa]\n"
        "A   : orifice area [m²]\n"
        "Cd  : discharge coefficient [-] (typically 0.6-0.65)\n"
        "rho : fluid density [kg/m³]\n\n"
        "Returns: mass flow rate [kg/s]"
    );

    m.def(
        "orifice_Q",
        &orifice_Q,
        py::arg("P1"),
        py::arg("P2"),
        py::arg("A"),
        py::arg("Cd"),
        py::arg("rho"),
        "Orifice volumetric flow rate.\n\n"
        "Q = Cd · A · √(2 · ΔP / ρ)\n\n"
        "Returns: volumetric flow rate [m³/s]"
    );

    m.def(
        "orifice_velocity",
        &orifice_velocity,
        py::arg("P1"),
        py::arg("P2"),
        py::arg("rho"),
        "Ideal orifice velocity (no losses).\n\n"
        "v = √(2 · ΔP / ρ)\n\n"
        "Returns: velocity [m/s]"
    );

    m.def(
        "orifice_area",
        &orifice_area,
        py::arg("mdot"),
        py::arg("P1"),
        py::arg("P2"),
        py::arg("Cd"),
        py::arg("rho"),
        "Solve for required orifice area given mass flow.\n\n"
        "A = ṁ / (Cd · √(2 · ρ · ΔP))\n\n"
        "Returns: orifice area [m²]"
    );

    m.def(
        "orifice_dP",
        static_cast<double(*)(double, double, double, double)>(&orifice_dP),
        py::arg("mdot"),
        py::arg("A"),
        py::arg("Cd"),
        py::arg("rho"),
        "Solve for pressure drop given mass flow and area.\n\n"
        "ΔP = (ṁ / (Cd · A))² / (2 · ρ)\n\n"
        "Returns: pressure drop P1-P2 [Pa]"
    );

    // Pipe flow
    m.def(
        "pipe_dP",
        &pipe_dP,
        py::arg("v"),
        py::arg("L"),
        py::arg("D"),
        py::arg("f"),
        py::arg("rho"),
        "Pipe pressure drop (Darcy-Weisbach).\n\n"
        "ΔP = f · (L/D) · (ρ · v² / 2)\n\n"
        "v   : flow velocity [m/s]\n"
        "L   : pipe length [m]\n"
        "D   : pipe diameter [m]\n"
        "f   : Darcy friction factor [-]\n"
        "rho : fluid density [kg/m³]\n\n"
        "Returns: pressure drop [Pa]"
    );

    m.def(
        "pipe_dP_mdot",
        &pipe_dP_mdot,
        py::arg("mdot"),
        py::arg("L"),
        py::arg("D"),
        py::arg("f"),
        py::arg("rho"),
        "Pipe pressure drop from mass flow rate.\n\n"
        "Returns: pressure drop [Pa]"
    );

    m.def(
        "pipe_velocity",
        &pipe_velocity,
        py::arg("mdot"),
        py::arg("D"),
        py::arg("rho"),
        "Pipe velocity from mass flow rate.\n\n"
        "v = ṁ / (ρ · π · D² / 4)\n\n"
        "Returns: velocity [m/s]"
    );

    m.def(
        "pipe_mdot",
        &pipe_mdot,
        py::arg("v"),
        py::arg("D"),
        py::arg("rho"),
        "Pipe mass flow from velocity.\n\n"
        "ṁ = ρ · v · π · D² / 4\n\n"
        "Returns: mass flow rate [kg/s]"
    );

    // Hydraulic utilities
    m.def(
        "dynamic_pressure",
        &dynamic_pressure,
        py::arg("v"),
        py::arg("rho"),
        "Dynamic pressure (velocity head).\n\n"
        "q = ½ · ρ · v²\n\n"
        "Returns: dynamic pressure [Pa]"
    );

    m.def(
        "velocity_from_q",
        &velocity_from_q,
        py::arg("q"),
        py::arg("rho"),
        "Velocity from dynamic pressure.\n\n"
        "v = √(2 · q / ρ)\n\n"
        "Returns: velocity [m/s]"
    );

    m.def(
        "hydraulic_diameter",
        &hydraulic_diameter,
        py::arg("A"),
        py::arg("P_wetted"),
        "Hydraulic diameter for non-circular cross-sections.\n\n"
        "Dh = 4 · A / P_wetted\n\n"
        "A        : cross-sectional area [m²]\n"
        "P_wetted : wetted perimeter [m]\n\n"
        "Returns: hydraulic diameter [m]"
    );

    m.def(
        "hydraulic_diameter_rect",
        &hydraulic_diameter_rect,
        py::arg("a"),
        py::arg("b"),
        "Hydraulic diameter for rectangular duct.\n\n"
        "Dh = 2·a·b / (a + b)\n\n"
        "Returns: hydraulic diameter [m]"
    );

    m.def(
        "hydraulic_diameter_annulus",
        &hydraulic_diameter_annulus,
        py::arg("D_outer"),
        py::arg("D_inner"),
        "Hydraulic diameter for annulus.\n\n"
        "Dh = D_outer - D_inner\n\n"
        "Returns: hydraulic diameter [m]"
    );

    // =========================================================================
    // Orifice Cd Correlations
    // =========================================================================

    // OrificeGeometry struct
    py::class_<OrificeGeometry>(m, "OrificeGeometry",
        "Orifice geometry for Cd correlations.\n\n"
        "Attributes:\n"
        "  d : orifice bore diameter [m]\n"
        "  D : pipe diameter [m]\n"
        "  t : plate thickness [m] (for thick plate)\n"
        "  r : inlet edge radius [m] (for rounded entry)")
        .def(py::init<>())
        .def_readwrite("d", &OrificeGeometry::d, "Orifice bore diameter [m]")
        .def_readwrite("D", &OrificeGeometry::D, "Pipe diameter [m]")
        .def_readwrite("t", &OrificeGeometry::t, "Plate thickness [m]")
        .def_readwrite("r", &OrificeGeometry::r, "Inlet edge radius [m]")
        .def_readwrite("bevel", &OrificeGeometry::bevel, "Bevel angle [rad]")
        .def_property_readonly("beta", &OrificeGeometry::beta, "Diameter ratio d/D [-]")
        .def_property_readonly("area", &OrificeGeometry::area, "Orifice area [m²]")
        .def_property_readonly("t_over_d", &OrificeGeometry::t_over_d, "Thickness ratio t/d [-]")
        .def_property_readonly("r_over_d", &OrificeGeometry::r_over_d, "Radius ratio r/d [-]")
        .def("is_valid", &OrificeGeometry::is_valid, "Check if geometry is valid");

    // OrificeState struct
    py::class_<OrificeState>(m, "OrificeState",
        "Flow state at orifice for Cd correlations.\n\n"
        "Attributes:\n"
        "  Re_D : pipe Reynolds number (based on D) [-]\n"
        "  dP   : differential pressure [Pa]\n"
        "  rho  : fluid density [kg/m³]\n"
        "  mu   : dynamic viscosity [Pa·s]")
        .def(py::init<>())
        .def_readwrite("Re_D", &OrificeState::Re_D, "Pipe Reynolds number [-]")
        .def_readwrite("dP", &OrificeState::dP, "Differential pressure [Pa]")
        .def_readwrite("rho", &OrificeState::rho, "Fluid density [kg/m³]")
        .def_readwrite("mu", &OrificeState::mu, "Dynamic viscosity [Pa·s]")
        .def("Re_d", &OrificeState::Re_d, py::arg("beta"),
             "Orifice Reynolds number (based on d)");

    // Cd correlation functions
    m.def(
        "Cd_sharp_thin_plate",
        &Cd_sharp_thin_plate,
        py::arg("geom"),
        py::arg("state"),
        "Discharge coefficient for sharp thin-plate orifice (ISO 5167-2).\n\n"
        "Uses Reader-Harris/Gallagher correlation.\n"
        "Valid for: 0.1 <= beta <= 0.75, Re_D >= 5000, D >= 50mm"
    );

    m.def(
        "Cd_thick_plate",
        &Cd_thick_plate,
        py::arg("geom"),
        py::arg("state"),
        "Discharge coefficient for thick-plate orifice.\n\n"
        "Applies Idelchik thickness correction to thin-plate Cd.\n"
        "Valid for: 0 < t/d < ~3"
    );

    m.def(
        "Cd_rounded_entry",
        &Cd_rounded_entry,
        py::arg("geom"),
        py::arg("state"),
        "Discharge coefficient for rounded-entry orifice.\n\n"
        "Based on Idelchik contraction loss coefficients.\n"
        "Valid for: 0 < r/d <= 0.2 (typical)"
    );

    m.def(
        "Cd_orifice",
        &Cd,
        py::arg("geom"),
        py::arg("state"),
        "Discharge coefficient with auto-selected correlation.\n\n"
        "Selects correlation based on geometry:\n"
        "  - r/d > 0.01: rounded entry\n"
        "  - t/d > 0.02: thick plate\n"
        "  - otherwise: sharp thin plate"
    );

    // Orifice flow calculations with Cd
    m.def(
        "orifice_mdot_Cd",
        static_cast<double(*)(const OrificeGeometry&, double, double, double)>(&orifice_mdot),
        py::arg("geom"),
        py::arg("Cd"),
        py::arg("dP"),
        py::arg("rho"),
        "Mass flow through orifice given Cd.\n\n"
        "mdot = Cd * A * sqrt(2 * rho * dP)\n\n"
        "Returns: mass flow rate [kg/s]"
    );

    m.def(
        "orifice_dP_Cd",
        static_cast<double(*)(const OrificeGeometry&, double, double, double)>(&orifice_dP),
        py::arg("geom"),
        py::arg("Cd"),
        py::arg("mdot"),
        py::arg("rho"),
        "Pressure drop for given mass flow and Cd.\n\n"
        "dP = (mdot / (Cd * A))² / (2 * rho)\n\n"
        "Returns: pressure drop [Pa]"
    );

    m.def(
        "orifice_Cd_from_measurement",
        &orifice_Cd_from_measurement,
        py::arg("geom"),
        py::arg("mdot"),
        py::arg("dP"),
        py::arg("rho"),
        "Solve for Cd given measured mass flow and pressure drop.\n\n"
        "Cd = mdot / (A * sqrt(2 * rho * dP))\n\n"
        "Returns: discharge coefficient [-]"
    );

    // Namespace functions
    m.def(
        "orifice_K_from_Cd",
        &orifice::K_from_Cd,
        py::arg("Cd"),
        py::arg("beta"),
        "Loss coefficient K from discharge coefficient Cd.\n\n"
        "K = (1/Cd² - 1) * (1 - beta⁴)"
    );

    m.def(
        "orifice_Cd_from_K",
        &orifice::Cd_from_K,
        py::arg("K"),
        py::arg("beta"),
        "Discharge coefficient Cd from loss coefficient K.\n\n"
        "Cd = 1 / sqrt(1 + K / (1 - beta⁴))"
    );

    m.def(
        "orifice_thickness_correction",
        &orifice::thickness_correction,
        py::arg("t_over_d"),
        py::arg("beta"),
        "Thickness correction factor for thick-plate orifices.\n\n"
        "Multiplies thin-plate Cd to account for flow reattachment.\n"
        "Returns 1.0 for t/d <= 0.02 (thin plate)."
    );

    // -------------------------------------------------------------------------
    // Units query API
    // -------------------------------------------------------------------------

    py::class_<combaero::units::UnitInfo>(m, "UnitInfo")
        .def_readonly("input", &combaero::units::UnitInfo::input, "Input units description")
        .def_readonly("output", &combaero::units::UnitInfo::output, "Output units description")
        .def("__repr__", [](const combaero::units::UnitInfo& u) {
            return "UnitInfo(input='" + std::string(u.input) + "', output='" + std::string(u.output) + "')";
        });

    m.def(
        "get_units",
        [](const std::string& name) -> py::object {
            auto u = combaero::units::get_units(name);
            if (u) {
                return py::cast(combaero::units::UnitInfo{u->input, u->output});
            }
            return py::none();
        },
        py::arg("function_name"),
        "Get input/output units for a function by name.\n\n"
        "Returns UnitInfo with 'input' and 'output' fields, or None if not found.\n\n"
        "Example:\n"
        "    >>> combaero.get_units('density')\n"
        "    UnitInfo(input='T: K, P: Pa, X: mol/mol', output='kg/m^3')"
    );

    m.def(
        "input_units",
        [](const std::string& name) {
            return std::string(combaero::units::input_units(name));
        },
        py::arg("function_name"),
        "Get input units string for a function (empty string if not found)."
    );

    m.def(
        "output_units",
        [](const std::string& name) {
            return std::string(combaero::units::output_units(name));
        },
        py::arg("function_name"),
        "Get output units string for a function (empty string if not found)."
    );

    m.def(
        "has_units",
        [](const std::string& name) {
            return combaero::units::has_units(name);
        },
        py::arg("function_name"),
        "Check if a function has registered units."
    );

    m.def(
        "list_functions_with_units",
        []() {
            std::vector<std::string> names;
            for (auto it = combaero::units::begin(); it != combaero::units::end(); ++it) {
                names.emplace_back(it->name);
            }
            return names;
        },
        "List all function names that have registered units."
    );

    m.def(
        "all_units",
        []() {
            py::dict result;
            for (auto it = combaero::units::begin(); it != combaero::units::end(); ++it) {
                py::dict entry;
                entry["input"] = std::string(it->input);
                entry["output"] = std::string(it->output);
                result[py::str(std::string(it->name))] = entry;
            }
            return result;
        },
        "Get all registered units as a dictionary.\n\n"
        "Returns: dict mapping function_name -> {'input': ..., 'output': ...}"
    );
}