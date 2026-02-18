#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>

#include "state.h"
#include "thermo.h"
#include "transport.h"
#include "utils.h"
#include "combustion.h"
#include "equilibrium.h"
#include "humidair.h"
#include "common_names.h"
#include "compressible.h"
#include "incompressible.h"
#include "friction.h"
#include "heat_transfer.h"
#include "acoustics.h"
#include "can_annular_solvers.h"
#include "orifice.h"
#include "pipe_flow.h"
#include "materials.h"
#include "cooling_correlations.h"
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
        "Mixture isochoric heat capacity Cv(T, X) [J/(mol*K)]."
    );

    m.def(
        "u",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return u(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture internal energy u(T, X) [J/mol].\n\n"
        "For ideal gas: u = h - RT\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: internal energy [J/mol]"
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
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
            auto X = to_vec(X_arr);
            return speed_of_sound(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Speed of sound [m/s].\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: speed of sound [m/s]"
    );

    m.def(
        "cp_mass",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
            auto X = to_vec(X_arr);
            return cp_mass(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mass-specific heat capacity at constant pressure.\n\n"
        "Cp_mass = Cp_molar / mwmix * 1000\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: heat capacity [J/(kg·K)]"
    );

    m.def(
        "cv_mass",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
            auto X = to_vec(X_arr);
            return cv_mass(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mass-specific heat capacity at constant volume.\n\n"
        "Cv_mass = Cv_molar / mwmix * 1000\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: heat capacity [J/(kg·K)]"
    );

    m.def(
        "h_mass",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
            auto X = to_vec(X_arr);
            return h_mass(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mass-specific enthalpy.\n\n"
        "h_mass = h_molar / mwmix * 1000\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: enthalpy [J/kg]"
    );

    m.def(
        "s_mass",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double P,
           double P_ref = 101325.0) {
            auto X = to_vec(X_arr);
            return s_mass(T, X, P, P_ref);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P"),
        py::arg("P_ref") = 101325.0,
        "Mass-specific entropy.\n\n"
        "s_mass = s_molar / mwmix * 1000\n\n"
        "Parameters:\n"
        "  T     : temperature [K]\n"
        "  X     : mole fractions [mol/mol]\n"
        "  P     : pressure [Pa]\n"
        "  P_ref : reference pressure [Pa] (default: 101325)\n\n"
        "Returns: entropy [J/(kg*K)]"
    );

    m.def(
        "u_mass",
        [](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
            auto X = to_vec(X_arr);
            return u_mass(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mass-specific internal energy.\n\n"
        "u_mass = u_molar / mwmix * 1000\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  X : mole fractions [mol/mol]\n\n"
        "Returns: internal energy [J/kg]"
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

    // TransportState struct binding - Bundle of transport properties
    py::class_<TransportState>(m, "TransportState",
        "Bundle of transport properties for a gas mixture.\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  T     : Temperature [K] (input, echoed back)\n"
        "  P     : Pressure [Pa] (input, echoed back)\n"
        "  rho   : Density [kg/m³]\n"
        "  mu    : Dynamic viscosity [Pa·s]\n"
        "  k     : Thermal conductivity [W/(m·K)]\n"
        "  nu    : Kinematic viscosity [m²/s]\n"
        "  alpha : Thermal diffusivity [m²/s]\n"
        "  Pr    : Prandtl number [-]\n"
        "  cp    : Specific heat at constant pressure [J/(kg·K)]")
        .def_readonly("T", &TransportState::T, "Temperature [K]")
        .def_readonly("P", &TransportState::P, "Pressure [Pa]")
        .def_readonly("rho", &TransportState::rho, "Density [kg/m³]")
        .def_readonly("mu", &TransportState::mu, "Dynamic viscosity [Pa·s]")
        .def_readonly("k", &TransportState::k, "Thermal conductivity [W/(m·K)]")
        .def_readonly("nu", &TransportState::nu, "Kinematic viscosity [m²/s]")
        .def_readonly("alpha", &TransportState::alpha, "Thermal diffusivity [m²/s]")
        .def_readonly("Pr", &TransportState::Pr, "Prandtl number [-]")
        .def_readonly("cp", &TransportState::cp, "Specific heat at constant pressure [J/(kg·K)]")
        .def("__repr__", [](const TransportState& s) {
            return "<TransportState: T=" + std::to_string(s.T) + " K, "
                   "P=" + std::to_string(s.P) + " Pa, "
                   "mu=" + std::to_string(s.mu) + " Pa·s, "
                   "k=" + std::to_string(s.k) + " W/(m·K), "
                   "Pr=" + std::to_string(s.Pr) + ">";
        });

    m.def(
        "transport_state",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return transport_state(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Compute all transport properties at once.\n\n"
        "Convenience function that computes all transport properties\n"
        "for a gas mixture in a single call. Returns TransportState struct\n"
        "with read-only attributes for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n"
        "  P : pressure [Pa]\n"
        "  X : mole fractions [-]\n\n"
        "Returns: TransportState object with attributes:\n"
        "  T, P, rho, mu, k, nu, alpha, Pr, cp\n\n"
        "Example:\n"
        "  >>> X = cb.standard_dry_air_composition()\n"
        "  >>> state = cb.transport_state(T=300, P=101325, X=X)\n"
        "  >>> print(state.mu)  # IDE autocomplete works!\n"
        "  1.85e-05\n"
        "  >>> print(state.Pr)\n"
        "  0.707"
    );

    // CompleteState struct binding - Bundle of thermo + transport properties
    py::class_<CompleteState>(m, "CompleteState",
        "Bundle of thermodynamic and transport properties for a gas mixture.\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  thermo    : ThermoState with all 16 thermodynamic properties\n"
        "  transport : TransportState with all 9 transport properties")
        .def_readonly("thermo", &CompleteState::thermo, "ThermoState with all thermodynamic properties")
        .def_readonly("transport", &CompleteState::transport, "TransportState with all transport properties")
        .def("__repr__", [](const CompleteState& s) {
            return "<CompleteState: T=" + std::to_string(s.thermo.T) + " K, "
                   "P=" + std::to_string(s.thermo.P) + " Pa, "
                   "rho=" + std::to_string(s.thermo.rho) + " kg/m³, "
                   "mu=" + std::to_string(s.transport.mu) + " Pa·s>";
        });

    m.def(
        "complete_state",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double P_ref)
        {
            auto X = to_vec(X_arr);
            return complete_state(T, P, X, P_ref);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("P_ref") = 101325.0,
        "Compute all thermodynamic and transport properties at once.\n\n"
        "Convenience function that computes all thermodynamic and transport\n"
        "properties for a gas mixture in a single call. Returns CompleteState\n"
        "struct with nested thermo and transport states.\n\n"
        "Parameters:\n"
        "  T     : temperature [K]\n"
        "  P     : pressure [Pa]\n"
        "  X     : mole fractions [-]\n"
        "  P_ref : reference pressure for entropy [Pa] (default: 101325.0)\n\n"
        "Returns: CompleteState object with nested attributes:\n"
        "  thermo.T, thermo.P, thermo.rho, thermo.cp, thermo.h, etc. (16 properties)\n"
        "  transport.mu, transport.k, transport.Pr, etc. (9 properties)\n\n"
        "Example:\n"
        "  >>> X = cb.standard_dry_air_composition()\n"
        "  >>> state = cb.complete_state(T=300, P=101325, X=X)\n"
        "  >>> print(state.thermo.h)      # Thermodynamic property\n"
        "  -103.6\n"
        "  >>> print(state.transport.mu)  # Transport property\n"
        "  1.68e-05"
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
        [](py::array_t<double> Y) { return mass_to_mole(to_vec(Y)); },
        py::arg("Y"),
        "Convert mass fractions to mole fractions."
    );

    m.def(
        "mwmix",
        [](py::array_t<double> X) { return mwmix(to_vec(X)); },
        py::arg("X"),
        "Calculate mixture molecular weight [g/mol] from mole fractions."
    );

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

    m.def(
        "fuel_lhv_molar",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           double reference_temperature)
        {
            auto X_fuel = to_vec(X_fuel_arr);
            return fuel_lhv_molar(X_fuel, reference_temperature);
        },
        py::arg("X_fuel"),
        py::arg("reference_temperature") = 298.15,
        "Fuel lower heating value (LHV) on molar basis [J/mol fuel mixture].\n\n"
        "Computed from complete combustion to CO2 + H2O(g) at reference temperature."
    );

    m.def(
        "fuel_lhv_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           double reference_temperature)
        {
            auto X_fuel = to_vec(X_fuel_arr);
            return fuel_lhv_mass(X_fuel, reference_temperature);
        },
        py::arg("X_fuel"),
        py::arg("reference_temperature") = 298.15,
        "Fuel lower heating value (LHV) on mass basis [J/kg fuel mixture].\n\n"
        "Computed from complete combustion to CO2 + H2O(g) at reference temperature."
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

    // AirProperties struct binding - Bundle of air properties
    py::class_<AirProperties>(m, "AirProperties",
        "Bundle of air properties computed from (T, P, humidity).\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  rho   : Density [kg/m³]\n"
        "  mu    : Dynamic viscosity [Pa·s]\n"
        "  k     : Thermal conductivity [W/(m·K)]\n"
        "  cp    : Specific heat at constant pressure [J/(kg·K)]\n"
        "  cv    : Specific heat at constant volume [J/(kg·K)]\n"
        "  Pr    : Prandtl number [-]\n"
        "  nu    : Kinematic viscosity [m²/s]\n"
        "  alpha : Thermal diffusivity [m²/s]\n"
        "  gamma : Heat capacity ratio [-]\n"
        "  a     : Speed of sound [m/s]")
        .def_readonly("rho", &AirProperties::rho, "Density [kg/m³]")
        .def_readonly("mu", &AirProperties::mu, "Dynamic viscosity [Pa·s]")
        .def_readonly("k", &AirProperties::k, "Thermal conductivity [W/(m·K)]")
        .def_readonly("cp", &AirProperties::cp, "Specific heat at constant pressure [J/(kg·K)]")
        .def_readonly("cv", &AirProperties::cv, "Specific heat at constant volume [J/(kg·K)]")
        .def_readonly("Pr", &AirProperties::Pr, "Prandtl number [-]")
        .def_readonly("nu", &AirProperties::nu, "Kinematic viscosity [m²/s]")
        .def_readonly("alpha", &AirProperties::alpha, "Thermal diffusivity [m²/s]")
        .def_readonly("gamma", &AirProperties::gamma, "Heat capacity ratio [-]")
        .def_readonly("a", &AirProperties::a, "Speed of sound [m/s]")
        .def("__repr__", [](const AirProperties& p) {
            return "<AirProperties: rho=" + std::to_string(p.rho) + " kg/m³, "
                   "mu=" + std::to_string(p.mu) + " Pa·s, "
                   "k=" + std::to_string(p.k) + " W/(m·K), "
                   "cp=" + std::to_string(p.cp) + " J/(kg·K), "
                   "Pr=" + std::to_string(p.Pr) + ">";
        });

    m.def(
        "air_properties",
        &air_properties,
        py::arg("T"),
        py::arg("P"),
        py::arg("humidity") = 0.0,
        "Compute all air properties at once.\n\n"
        "Convenience function that computes all thermodynamic and transport\n"
        "properties of air in a single call. Returns AirProperties struct\n"
        "with read-only attributes for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  T        : temperature [K]\n"
        "  P        : pressure [Pa]\n"
        "  humidity : relative humidity [0-1] (default: 0.0 = dry air)\n\n"
        "Returns: AirProperties object with attributes:\n"
        "  rho, mu, k, cp, cv, Pr, nu, alpha, gamma, a\n\n"
        "Example:\n"
        "  >>> props = cb.air_properties(T=300, P=101325, humidity=0.5)\n"
        "  >>> print(props.rho)  # IDE autocomplete works!\n"
        "  1.177\n"
        "  >>> print(props.Pr)\n"
        "  0.707"
    );

    // ThermoState struct binding - Bundle of thermodynamic properties
    py::class_<ThermoState>(m, "ThermoState",
        "Bundle of thermodynamic properties for a gas mixture.\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  T       : Temperature [K] (input, echoed back)\n"
        "  P       : Pressure [Pa] (input, echoed back)\n"
        "  rho     : Density [kg/m³]\n"
        "  cp      : Specific heat at constant pressure [J/(mol·K)]\n"
        "  cv      : Specific heat at constant volume [J/(mol·K)]\n"
        "  h       : Specific enthalpy [J/mol]\n"
        "  s       : Specific entropy [J/(mol·K)]\n"
        "  u       : Specific internal energy [J/mol]\n"
        "  gamma   : Isentropic expansion coefficient [-]\n"
        "  a       : Speed of sound [m/s]\n"
        "  cp_mass : Mass-specific cp [J/(kg·K)]\n"
        "  cv_mass : Mass-specific cv [J/(kg·K)]\n"
        "  h_mass  : Mass-specific enthalpy [J/kg]\n"
        "  s_mass  : Mass-specific entropy [J/(kg·K)]\n"
        "  u_mass  : Mass-specific internal energy [J/kg]\n"
        "  mw      : Molecular weight [g/mol]")
        .def_readonly("T", &ThermoState::T, "Temperature [K]")
        .def_readonly("P", &ThermoState::P, "Pressure [Pa]")
        .def_readonly("rho", &ThermoState::rho, "Density [kg/m³]")
        .def_readonly("cp", &ThermoState::cp, "Specific heat at constant pressure [J/(mol·K)]")
        .def_readonly("cv", &ThermoState::cv, "Specific heat at constant volume [J/(mol·K)]")
        .def_readonly("h", &ThermoState::h, "Specific enthalpy [J/mol]")
        .def_readonly("s", &ThermoState::s, "Specific entropy [J/(mol·K)]")
        .def_readonly("u", &ThermoState::u, "Specific internal energy [J/mol]")
        .def_readonly("gamma", &ThermoState::gamma, "Isentropic expansion coefficient [-]")
        .def_readonly("a", &ThermoState::a, "Speed of sound [m/s]")
        .def_readonly("cp_mass", &ThermoState::cp_mass, "Mass-specific cp [J/(kg·K)]")
        .def_readonly("cv_mass", &ThermoState::cv_mass, "Mass-specific cv [J/(kg·K)]")
        .def_readonly("h_mass", &ThermoState::h_mass, "Mass-specific enthalpy [J/kg]")
        .def_readonly("s_mass", &ThermoState::s_mass, "Mass-specific entropy [J/(kg·K)]")
        .def_readonly("u_mass", &ThermoState::u_mass, "Mass-specific internal energy [J/kg]")
        .def_readonly("mw", &ThermoState::mw, "Molecular weight [g/mol]")
        .def("__repr__", [](const ThermoState& s) {
            return "<ThermoState: T=" + std::to_string(s.T) + " K, "
                   "P=" + std::to_string(s.P) + " Pa, "
                   "rho=" + std::to_string(s.rho) + " kg/m³>";
        });

    m.def(
        "thermo_state",
        &thermo_state,
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("P_ref") = 101325.0,
        "Compute all thermodynamic properties at once.\n\n"
        "Convenience function that computes all thermodynamic properties\n"
        "for a gas mixture in a single call. Returns ThermoState struct\n"
        "with read-only attributes for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  T     : temperature [K]\n"
        "  P     : pressure [Pa]\n"
        "  X     : mole fractions [-]\n"
        "  P_ref : reference pressure for entropy [Pa] (default: 101325.0)\n\n"
        "Returns: ThermoState object with attributes:\n"
        "  T, P, rho, cp, cv, h, s, u, gamma, a,\n"
        "  cp_mass, cv_mass, h_mass, s_mass, u_mass, mw\n\n"
        "Example:\n"
        "  >>> X = cb.standard_dry_air_composition()\n"
        "  >>> state = cb.thermo_state(T=300, P=101325, X=X)\n"
        "  >>> print(state.h)  # IDE autocomplete works!\n"
        "  8682.5\n"
        "  >>> print(state.gamma)\n"
        "  1.400"
    );

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
        "set_fuel_stream_for_phi",
        &set_fuel_stream_for_phi,
        py::arg("phi"),
        py::arg("fuel"),
        py::arg("oxidizer"),
        "Find fuel stream mass flow rate for target equivalence ratio phi.\n\n"
        "This is a direct algebraic helper (no iteration):\n"
        "  fuel.mdot = phi * fuel_stoich_mdot(oxidizer stream)\n\n"
        "Parameters:\n"
        "  phi      : target equivalence ratio [-]\n"
        "  fuel     : fuel Stream (T, P, X set; mdot ignored)\n"
        "  oxidizer : oxidizer Stream (with mdot set)\n\n"
        "Returns: fuel Stream with mdot set for requested phi."
    );

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
    // Friction factor correlations
    // -------------------------------------------------------------

    m.def(
        "friction_haaland",
        &friction_haaland,
        py::arg("Re"),
        py::arg("e_D"),
        "Haaland friction factor correlation (1983).\n\n"
        "Explicit approximation to Colebrook-White equation.\n"
        "Accuracy: ~2-3% vs Colebrook-White.\n\n"
        "f = 1 / [-1.8 * log10((e_D/3.7)^1.11 + 6.9/Re)]^2\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-] (must be > 0)\n"
        "  e_D : relative roughness ε/D [-] (must be >= 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for turbulent flow (Re > ~2300).\n"
        "For laminar flow, use f = 64/Re."
    );

    m.def(
        "friction_serghides",
        &friction_serghides,
        py::arg("Re"),
        py::arg("e_D"),
        "Serghides friction factor correlation (1984).\n\n"
        "Explicit approximation using Steffensen acceleration.\n"
        "Accuracy: <0.3% vs Colebrook-White (very accurate).\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-] (must be > 0)\n"
        "  e_D : relative roughness ε/D [-] (must be >= 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for turbulent flow (Re > ~2300)."
    );

    m.def(
        "friction_colebrook",
        &friction_colebrook,
        py::arg("Re"),
        py::arg("e_D"),
        py::arg("tol") = 1e-10,
        py::arg("max_iter") = 20,
        "Colebrook-White friction factor equation (1939).\n\n"
        "Implicit equation solved iteratively using Newton-Raphson.\n"
        "The reference standard for turbulent friction factor.\n\n"
        "1/√f = -2 * log10(ε/D/3.7 + 2.51/(Re*√f))\n\n"
        "Parameters:\n"
        "  Re       : Reynolds number [-] (must be > 0)\n"
        "  e_D      : relative roughness ε/D [-] (must be >= 0)\n"
        "  tol      : convergence tolerance [-] (default: 1e-10)\n"
        "  max_iter : maximum iterations (default: 20)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Uses Haaland correlation as initial guess."
    );

    m.def(
        "friction_petukhov",
        &friction_petukhov,
        py::arg("Re"),
        "Petukhov friction factor correlation (1970).\n\n"
        "For smooth pipes only (zero roughness).\n"
        "f = [0.790 * ln(Re) - 1.64]^(-2)\n\n"
        "Parameters:\n"
        "  Re : Reynolds number [-] (must be > 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for: 3000 < Re < 5x10^6\n"
        "Often used with Gnielinski/Petukhov heat transfer correlations."
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
        "Dittus-Boelter Nusselt number correlation (1930).\n\n"
        "Nu = 0.023 * Re^0.8 * Pr^n\n"
        "where n = 0.4 for heating, n = 0.3 for cooling.\n\n"
        "Parameters:\n"
        "  Re      : Reynolds number [-] (must be > 10,000)\n"
        "  Pr      : Prandtl number [-] (must be 0.6 < Pr < 160)\n"
        "  heating : True if fluid is being heated, False if cooled (default: True)\n\n"
        "Returns: Nusselt number Nu [-] (dimensionless)\n\n"
        "Valid for fully developed turbulent flow with L/D > 10.\n"
        "Heat transfer coefficient: h = Nu * k / L [W/(m²·K)]"
    );

    m.def(
        "nusselt_gnielinski",
        py::overload_cast<double, double, double>(&nusselt_gnielinski),
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("f"),
        "Gnielinski Nusselt number correlation (1976) with friction factor.\n\n"
        "Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - 1))\n\n"
        "More accurate than Dittus-Boelter, especially in transition region.\n\n"
        "Parameters:\n"
        "  Re : Reynolds number [-] (2300 < Re < 5x10^6)\n"
        "  Pr : Prandtl number [-] (0.5 < Pr < 2000)\n"
        "  f  : Darcy friction factor [-] (use friction_colebrook or similar)\n\n"
        "Returns: Nusselt number Nu [-] (dimensionless)"
    );

    m.def(
        "nusselt_gnielinski",
        py::overload_cast<double, double>(&nusselt_gnielinski),
        py::arg("Re"),
        py::arg("Pr"),
        "Gnielinski Nusselt number correlation (1976) with automatic friction.\n\n"
        "Uses Petukhov friction correlation for smooth pipes:\n"
        "f = [0.790 * ln(Re) - 1.64]^(-2)\n\n"
        "Parameters:\n"
        "  Re : Reynolds number [-] (2300 < Re < 5x10^6)\n"
        "  Pr : Prandtl number [-] (0.5 < Pr < 2000)\n\n"
        "Returns: Nusselt number Nu [-] (dimensionless)"
    );

    m.def(
        "nusselt_sieder_tate",
        &nusselt_sieder_tate,
        py::arg("Re"),
        py::arg("Pr"),
        py::arg("mu_ratio") = 1.0,
        "Sieder-Tate Nusselt number correlation (1936).\n\n"
        "Nu = 0.027 * Re^0.8 * Pr^(1/3) * (μ_bulk / μ_wall)^0.14\n\n"
        "Accounts for viscosity variation between bulk and wall temperatures.\n\n"
        "Parameters:\n"
        "  Re       : Reynolds number at bulk temperature [-] (Re > 10,000)\n"
        "  Pr       : Prandtl number at bulk temperature [-] (0.7 < Pr < 16,700)\n"
        "  mu_ratio : μ_bulk / μ_wall [-] (default: 1.0, typically 0.5-2.0)\n\n"
        "Returns: Nusselt number Nu [-] (dimensionless)\n\n"
        "Valid for L/D > 10."
    );

    m.def(
        "htc_from_nusselt",
        &htc_from_nusselt,
        py::arg("Nu"),
        py::arg("k"),
        py::arg("L"),
        "Heat transfer coefficient from Nusselt number.\n\n"
        "h = Nu * k / L\n\n"
        "Parameters:\n"
        "  Nu : Nusselt number [-] (dimensionless)\n"
        "  k  : thermal conductivity [W/(m·K)]\n"
        "  L  : characteristic length [m] (diameter for pipe flow)\n\n"
        "Returns: heat transfer coefficient h [W/(m²·K)]"
    );

    m.def(
        "htc_pipe",
        [](double T, double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double velocity, double diameter,
           const std::string& correlation,
           bool heating,
           double mu_ratio,
           double roughness) {
            auto X = to_vec(X_arr);
            return htc_pipe(T, P, X, velocity, diameter, correlation, heating, mu_ratio, roughness);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("velocity"),
        py::arg("diameter"),
        py::arg("correlation") = "gnielinski",
        py::arg("heating") = true,
        py::arg("mu_ratio") = 1.0,
        py::arg("roughness") = 0.0,
        "Composite heat transfer coefficient for pipe flow.\n\n"
        "Computes HTC, Nusselt number, and Reynolds number from thermodynamic state.\n"
        "Automatically calculates density, viscosity, thermal conductivity, and Prandtl number.\n\n"
        "Parameters:\n"
        "  T           : temperature [K]\n"
        "  P           : pressure [Pa]\n"
        "  X           : mole fractions [mol/mol]\n"
        "  velocity    : flow velocity [m/s]\n"
        "  diameter    : pipe diameter [m]\n"
        "  correlation : 'gnielinski' (default), 'dittus_boelter', 'sieder_tate', 'petukhov'\n"
        "  heating     : True for heating, False for cooling (affects Dittus-Boelter)\n"
        "  mu_ratio    : μ_bulk / μ_wall for Sieder-Tate viscosity correction (default: 1.0)\n"
        "  roughness   : absolute roughness [m] (default: 0.0 = smooth pipe)\n\n"
        "Returns: tuple (h, Nu, Re)\n"
        "  h  : heat transfer coefficient [W/(m²·K)]\n"
        "  Nu : Nusselt number [-]\n"
        "  Re : Reynolds number [-]\n\n"
        "Correlations:\n"
        "  - gnielinski    : 2300 < Re < 5e6, 0.5 < Pr < 2000 (best general-purpose)\n"
        "  - dittus_boelter: Re > 10000, 0.6 < Pr < 160 (simple, fast)\n"
        "  - sieder_tate   : Re > 10000, 0.7 < Pr < 16700 (viscosity correction)\n"
        "  - petukhov      : Re > 10000, 0.5 < Pr < 2000 (high accuracy)\n\n"
        "Automatically handles laminar flow (Re < 2300) with constant Nu."
    );

    m.def(
        "lmtd",
        &lmtd,
        py::arg("dT1"),
        py::arg("dT2"),
        "Log mean temperature difference (LMTD) for heat exchangers.\n\n"
        "LMTD = (ΔT1 - ΔT2) / ln(ΔT1 / ΔT2)\n\n"
        "For counter-flow:\n"
        "  ΔT1 = T_hot_in - T_cold_out\n"
        "  ΔT2 = T_hot_out - T_cold_in\n\n"
        "For parallel-flow:\n"
        "  ΔT1 = T_hot_in - T_cold_in\n"
        "  ΔT2 = T_hot_out - T_cold_out\n\n"
        "Parameters:\n"
        "  dT1 : temperature difference at one end [K]\n"
        "  dT2 : temperature difference at other end [K]\n\n"
        "Returns: log mean temperature difference [K]\n\n"
        "Note: If dT1 ≈ dT2, returns arithmetic mean to avoid 0/0."
    );

    // -------------------------------------------------------------
    // Geometry utilities
    // -------------------------------------------------------------

    m.def(
        "hydraulic_diameter",
        &hydraulic_diameter,
        py::arg("A"),
        py::arg("P_wetted"),
        "Hydraulic diameter for arbitrary cross-section.\n\n"
        "Dh = 4 * A / P_wetted\n\n"
        "Parameters:\n"
        "  A        : cross-sectional area [m²]\n"
        "  P_wetted : wetted perimeter [m]\n\n"
        "Returns: hydraulic diameter Dh [m]"
    );

    m.def(
        "hydraulic_diameter_rect",
        &hydraulic_diameter_rect,
        py::arg("a"),
        py::arg("b"),
        "Hydraulic diameter for rectangular duct.\n\n"
        "Dh = 2*a*b / (a + b)\n\n"
        "Parameters:\n"
        "  a : side length [m]\n"
        "  b : side length [m]\n\n"
        "Returns: hydraulic diameter Dh [m]"
    );

    m.def(
        "hydraulic_diameter_annulus",
        &hydraulic_diameter_annulus,
        py::arg("D_outer"),
        py::arg("D_inner"),
        "Hydraulic diameter for annular duct (concentric pipes).\n\n"
        "Dh = D_outer - D_inner\n\n"
        "Parameters:\n"
        "  D_outer : outer diameter [m]\n"
        "  D_inner : inner diameter [m]\n\n"
        "Returns: hydraulic diameter Dh [m]"
    );

    m.def(
        "residence_time",
        py::overload_cast<double, double>(&residence_time),
        py::arg("V"),
        py::arg("Q"),
        "Residence time from volumetric flow rate.\n\n"
        "τ = V / Q̇\n\n"
        "Time for fluid to pass through a volume.\n"
        "Used in reactor design (Damköhler number), combustor sizing.\n\n"
        "Parameters:\n"
        "  V : volume [m³]\n"
        "  Q : volumetric flow rate [m³/s]\n\n"
        "Returns: residence time τ [s]"
    );

    m.def(
        "residence_time_mdot",
        py::overload_cast<double, double, double>(&residence_time_mdot),
        py::arg("V"),
        py::arg("mdot"),
        py::arg("rho"),
        "Residence time from mass flow rate.\n\n"
        "τ = V·ρ / ṁ\n\n"
        "Parameters:\n"
        "  V    : volume [m³]\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  rho  : density [kg/m³]\n\n"
        "Returns: residence time τ [s]"
    );

    m.def(
        "pipe_mdot",
        &pipe_mdot,
        py::arg("v"),
        py::arg("D"),
        py::arg("rho"),
        "Pipe mass flow rate from velocity.\n\n"
        "mdot = rho * v * pi * D^2 / 4\n\n"
        "Parameters:\n"
        "  v   : velocity [m/s]\n"
        "  D   : diameter [m]\n"
        "  rho : density [kg/m^3]\n\n"
        "Returns: mass flow rate [kg/s]"
    );

    // Composite pipe flow function
    m.def(
        "pressure_drop_pipe",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double v,
           double D,
           double L,
           double roughness,
           const std::string& correlation) {
            auto X = to_vec(X_arr);
            return pressure_drop_pipe(T, P, X, v, D, L, roughness, correlation);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        py::arg("v"),
        py::arg("D"),
        py::arg("L"),
        py::arg("roughness") = 0.0,
        py::arg("correlation") = "haaland",
        "Composite pressure drop calculation for pipe flow.\n\n"
        "Combines thermodynamic properties, Reynolds number, friction factor,\n"
        "and Darcy-Weisbach equation in one call.\n\n"
        "Steps:\n"
        "  1. Calculate density and viscosity from (T, P, X)\n"
        "  2. Calculate Reynolds number: Re = rho * v * D / mu\n"
        "  3. Calculate friction factor using specified correlation\n"
        "  4. Calculate pressure drop: dP = f * (L/D) * (rho * v^2 / 2)\n\n"
        "Parameters:\n"
        "  T           : temperature [K]\n"
        "  P           : pressure [Pa]\n"
        "  X           : mole fractions [mol/mol]\n"
        "  v           : velocity [m/s]\n"
        "  D           : diameter [m]\n"
        "  L           : length [m]\n"
        "  roughness   : absolute roughness [m] (default: 0.0 = smooth)\n"
        "  correlation : friction correlation (default: 'haaland')\n"
        "                Options: 'haaland', 'serghides', 'colebrook', 'petukhov'\n\n"
        "Returns: tuple of (dP [Pa], Re [-], f [-])\n\n"
        "Example:\n"
        "  >>> dP, Re, f = pressure_drop_pipe(300, 101325, X_air, 10.0, 0.1, 100.0)\n"
        "  >>> print(f'Pressure drop: {dP:.1f} Pa, Re: {Re:.0f}, f: {f:.4f}')"
    );

    m.def(
        "space_velocity",
        &space_velocity,
        py::arg("Q"),
        py::arg("V"),
        "Space velocity (inverse of residence time).\n\n"
        "SV = Q̇ / V = 1/τ\n\n"
        "Parameters:\n"
        "  Q : volumetric flow rate [m³/s]\n"
        "  V : volume [m³]\n\n"
        "Returns: space velocity SV [1/s]"
    );

    m.def(
        "pipe_area",
        &pipe_area,
        py::arg("D"),
        "Circular pipe cross-sectional area.\n\n"
        "A = π * (D/2)²\n\n"
        "Parameters:\n"
        "  D : diameter [m]\n\n"
        "Returns: area [m²]"
    );

    m.def(
        "annular_area",
        &annular_area,
        py::arg("D_outer"),
        py::arg("D_inner"),
        "Annular cross-sectional area.\n\n"
        "A = π * ((D_outer/2)² - (D_inner/2)²)\n\n"
        "Parameters:\n"
        "  D_outer : outer diameter [m]\n"
        "  D_inner : inner diameter [m]\n\n"
        "Returns: area [m²]"
    );

    m.def(
        "pipe_volume",
        &pipe_volume,
        py::arg("D"),
        py::arg("L"),
        "Cylindrical pipe volume.\n\n"
        "V = π * (D/2)² * L\n\n"
        "Parameters:\n"
        "  D : diameter [m]\n"
        "  L : length [m]\n\n"
        "Returns: volume [m³]"
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

    // =========================================================================
    // Combustion State
    // =========================================================================

    py::class_<CombustionState>(m, "CombustionState",
        "Bundle of combustion properties with nested CompleteState objects.\n\n"
        "Contains reactant and product states with all thermodynamic and transport properties.\n\n"
        "Attributes:\n"
        "  phi                : equivalence ratio [-]\n"
        "  fuel_name          : optional fuel label (string)\n"
        "  reactants          : CompleteState at reactant conditions\n"
        "  products           : CompleteState at product conditions (T_ad)\n"
        "  mixture_fraction   : Bilger mixture fraction [-]\n"
        "  fuel_burn_fraction : fraction of fuel burned [0-1]")
        .def_readonly("phi", &CombustionState::phi, "Equivalence ratio [-]")
        .def_readonly("fuel_name", &CombustionState::fuel_name, "Fuel label (e.g., 'CH4', 'Natural Gas')")
        .def_readonly("reactants", &CombustionState::reactants, "Reactant CompleteState (all thermo + transport)")
        .def_readonly("products", &CombustionState::products, "Product CompleteState at T_ad (all thermo + transport)")
        .def_readonly("mixture_fraction", &CombustionState::mixture_fraction, "Bilger mixture fraction [-]")
        .def_readonly("fuel_burn_fraction", &CombustionState::fuel_burn_fraction, "Fraction of fuel burned [0-1]")
        .def("__repr__", [](const CombustionState& s) {
            return "<CombustionState: phi=" + std::to_string(s.phi) +
                   ", T_reactants=" + std::to_string(s.reactants.thermo.T) + " K" +
                   ", T_products=" + std::to_string(s.products.thermo.T) + " K" +
                   (s.fuel_name.empty() ? "" : ", fuel='" + s.fuel_name + "'") + ">";
        });

    m.def(
        "combustion_state",
        &combustion_state,
        py::arg("X_fuel"),
        py::arg("X_ox"),
        py::arg("phi"),
        py::arg("T_reactants"),
        py::arg("P"),
        py::arg("fuel_name") = "",
        "Compute combustion state from equivalence ratio.\n\n"
        "Typical use: calculations where phi is specified.\n\n"
        "Parameters:\n"
        "  X_fuel       : fuel composition (mole fractions) [-]\n"
        "  X_ox         : oxidizer composition (mole fractions) [-]\n"
        "  phi          : equivalence ratio [-] (INPUT)\n"
        "  T_reactants  : reactant temperature [K]\n"
        "  P            : pressure [Pa]\n"
        "  fuel_name    : optional fuel label (default: '')\n\n"
        "Returns: CombustionState with reactants, products, phi, mixture_fraction\n\n"
        "Example:\n"
        "  >>> X_CH4 = [0]*5 + [1.0] + [0]*8  # Pure methane\n"
        "  >>> X_air = cb.standard_dry_air_composition()\n"
        "  >>> state = cb.combustion_state(X_CH4, X_air, phi=1.0, T_reactants=300, P=101325)\n"
        "  >>> state.products.thermo.T  # Adiabatic flame temperature [K]"
    );

    m.def(
        "combustion_state_from_streams",
        &combustion_state_from_streams,
        py::arg("fuel_stream"),
        py::arg("ox_stream"),
        py::arg("fuel_name") = "",
        "Compute combustion state from measured streams.\n\n"
        "Typical use: lab measurements where mass flows are measured.\n"
        "Phi is COMPUTED from mass flow rates (output, not input).\n\n"
        "Parameters:\n"
        "  fuel_stream : fuel stream with mdot, T, X\n"
        "  ox_stream   : oxidizer stream with mdot, T, X\n"
        "  fuel_name   : optional fuel label (default: '')\n\n"
        "Returns: CombustionState with phi computed from flow rates\n\n"
        "Example:\n"
        "  >>> fuel = cb.Stream().set_T(300).set_X(X_CH4).set_mdot(0.01)\n"
        "  >>> air = cb.Stream().set_T(298).set_X(X_air).set_mdot(0.17)\n"
        "  >>> state = cb.combustion_state_from_streams(fuel, air, 'CH4')\n"
        "  >>> state.phi  # Computed from mass flows"
    );

    // =========================================================================
    // Acoustics
    // =========================================================================

    // Tube struct
    py::class_<Tube>(m, "Tube", "Cylindrical tube geometry for acoustic analysis")
        .def(py::init<double, double>(), py::arg("L"), py::arg("D"),
             "Create tube with length L [m] and diameter D [m]")
        .def_readwrite("L", &Tube::L, "Length [m]")
        .def_readwrite("D", &Tube::D, "Diameter [m]")
        .def("area", &Tube::area, "Cross-sectional area [m²]")
        .def("volume", &Tube::volume, "Volume [m³]")
        .def("perimeter", &Tube::perimeter, "Circumference [m]");

    // CdCorrelation enum (early binding for default arguments)
    py::enum_<CdCorrelation>(m, "CdCorrelation", "Orifice discharge coefficient correlations")
        .value("ReaderHarrisGallagher", CdCorrelation::ReaderHarrisGallagher, "ISO 5167-2 / ASME MFC-3M")
        .value("Stolz", CdCorrelation::Stolz, "ISO 5167:1980 (older)")
        .value("Miller", CdCorrelation::Miller, "Miller (1996) simplified");

    // Annulus struct
    py::class_<Annulus>(m, "Annulus", "Annular duct geometry for acoustic analysis")
        .def(py::init<double, double, double>(),
             py::arg("L"), py::arg("D_inner"), py::arg("D_outer"),
             "Create annulus with length L [m], inner diameter D_inner [m], outer diameter D_outer [m]")
        .def_readwrite("L", &Annulus::L, "Length [m]")
        .def_readwrite("D_inner", &Annulus::D_inner, "Inner diameter [m]")
        .def_readwrite("D_outer", &Annulus::D_outer, "Outer diameter [m]")
        .def_readwrite("n_azimuthal_max", &Annulus::n_azimuthal_max, "Maximum azimuthal mode number to search")
        .def("D_mean", &Annulus::D_mean, "Mean diameter [m]")
        .def("gap", &Annulus::gap, "Annular gap [m]")
        .def("area", &Annulus::area, "Cross-sectional area [m²]")
        .def("volume", &Annulus::volume, "Volume [m³]")
        .def("circumference", &Annulus::circumference, "Mean circumference [m]")
        .def("radius_inner", &Annulus::radius_inner, "Inner radius [m]")
        .def("radius_outer", &Annulus::radius_outer, "Outer radius [m]");

    py::class_<CanAnnularFlowGeometry>(
        m,
        "CanAnnularFlowGeometry",
        "Parameterized can-annular flow geometry with primary, transition, and annular sections"
    )
        .def(
            py::init<double, double, double, double, double, double>(),
            py::arg("L"),
            py::arg("D_inner"),
            py::arg("D_outer"),
            py::arg("L_primary"),
            py::arg("D_primary"),
            py::arg("L_transition"),
            "Create geometry with annular section and dedicated can transition parameters."
        )
        .def_readwrite("L", &CanAnnularFlowGeometry::L, "Annular section length [m]")
        .def_readwrite("D_inner", &CanAnnularFlowGeometry::D_inner, "Annular inner diameter [m]")
        .def_readwrite("D_outer", &CanAnnularFlowGeometry::D_outer, "Annular outer diameter [m]")
        .def_readwrite("L_primary", &CanAnnularFlowGeometry::L_primary, "Primary zone length [m]")
        .def_readwrite("D_primary", &CanAnnularFlowGeometry::D_primary, "Primary zone diameter [m]")
        .def_readwrite("L_transition", &CanAnnularFlowGeometry::L_transition, "Transition length [m]")
        .def("D_mean", &CanAnnularFlowGeometry::D_mean, "Annular mean diameter [m]")
        .def("gap", &CanAnnularFlowGeometry::gap, "Annular gap [m]")
        .def("area", &CanAnnularFlowGeometry::area, "Annular cross-sectional area [m²]")
        .def("volume", &CanAnnularFlowGeometry::volume, "Annular section volume [m³]")
        .def("circumference", &CanAnnularFlowGeometry::circumference, "Annular mean circumference [m]")
        .def("area_primary", &CanAnnularFlowGeometry::area_primary, "Primary circular area [m²]")
        .def("volume_primary", &CanAnnularFlowGeometry::volume_primary, "Primary zone volume [m³]")
        .def("volume_transition", &CanAnnularFlowGeometry::volume_transition, "Transition volume [m³]")
        .def("volume_total", &CanAnnularFlowGeometry::volume_total, "Total geometry volume [m³]")
        .def("length_total", &CanAnnularFlowGeometry::length_total, "Total geometry length [m]");

    // BoundaryCondition enum
    py::enum_<BoundaryCondition>(m, "BoundaryCondition", "Acoustic boundary condition")
        .value("Open", BoundaryCondition::Open, "Pressure release (p'=0)")
        .value("Closed", BoundaryCondition::Closed, "Rigid wall (u'=0)");

    // AcousticMode struct
    py::class_<AcousticMode>(m, "AcousticMode", "Acoustic mode with frequency and indices")
        .def_readonly("n_axial", &AcousticMode::n_axial, "Longitudinal mode number")
        .def_readonly("n_azimuthal", &AcousticMode::n_azimuthal, "Azimuthal mode number")
        .def_readonly("frequency", &AcousticMode::frequency, "Natural frequency [Hz]")
        .def("label", &AcousticMode::label, "Mode label (e.g., '1L', '2T', '1L1T')")
        .def("__repr__", [](const AcousticMode& m) {
            return "<AcousticMode " + m.label() + " @ " +
                   std::to_string(m.frequency) + " Hz>";
        });

    // Tube modes
    m.def(
        "tube_axial_modes",
        &tube_axial_modes,
        py::arg("tube"),
        py::arg("c"),
        py::arg("upstream"),
        py::arg("downstream"),
        py::arg("n_max") = 5,
        "Compute axial acoustic modes of a tube.\n\n"
        "tube      : Tube geometry\n"
        "c         : Speed of sound [m/s]\n"
        "upstream  : Boundary condition at x=0\n"
        "downstream: Boundary condition at x=L\n"
        "n_max     : Maximum mode number (default: 5)\n\n"
        "Returns list of AcousticMode sorted by frequency."
    );

    // Annulus modes
    m.def(
        "annulus_axial_modes",
        &annulus_axial_modes,
        py::arg("annulus"),
        py::arg("c"),
        py::arg("upstream"),
        py::arg("downstream"),
        py::arg("n_max") = 5,
        "Compute axial acoustic modes of an annulus."
    );

    m.def(
        "annulus_azimuthal_modes",
        &annulus_azimuthal_modes,
        py::arg("annulus"),
        py::arg("c"),
        py::arg("m_max") = 4,
        "Compute azimuthal (tangential) acoustic modes of an annulus.\n\n"
        "f_m = m * c / (π * D_mean)"
    );

    m.def(
        "annulus_modes",
        &annulus_modes,
        py::arg("annulus"),
        py::arg("c"),
        py::arg("upstream"),
        py::arg("downstream"),
        py::arg("n_max") = 5,
        py::arg("m_max") = 4,
        "Compute all acoustic modes (axial, azimuthal, combined) of an annulus.\n\n"
        "Returns list of AcousticMode sorted by frequency."
    );

    // Utility functions
    m.def(
        "modes_in_range",
        &modes_in_range,
        py::arg("modes"),
        py::arg("f_min"),
        py::arg("f_max"),
        "Filter modes within a frequency range [f_min, f_max] Hz."
    );

    m.def(
        "closest_mode",
        [](const std::vector<AcousticMode>& modes, double f_target) -> py::object {
            const auto* m = closest_mode(modes, f_target);
            if (m) return py::cast(*m);
            return py::none();
        },
        py::arg("modes"),
        py::arg("f_target"),
        "Find the mode closest to a target frequency.\n\n"
        "Returns AcousticMode or None if modes is empty."
    );

    m.def(
        "min_mode_separation",
        &min_mode_separation,
        py::arg("modes"),
        "Minimum frequency separation between adjacent modes [Hz]."
    );

    // Mean flow correction
    m.def(
        "axial_mode_upstream",
        &axial_mode_upstream,
        py::arg("f0"),
        py::arg("M"),
        "Upstream-propagating wave frequency with mean flow.\n\n"
        "f+ = f0 / (1 - M)\n\n"
        "f0 : quiescent frequency [Hz]\n"
        "M  : Mach number [-]"
    );

    m.def(
        "axial_mode_downstream",
        &axial_mode_downstream,
        py::arg("f0"),
        py::arg("M"),
        "Downstream-propagating wave frequency with mean flow.\n\n"
        "f- = f0 / (1 + M)"
    );

    m.def(
        "axial_mode_split",
        &axial_mode_split,
        py::arg("f0"),
        py::arg("M"),
        "Frequency split due to mean flow.\n\n"
        "Returns (f_upstream, f_downstream)."
    );

    // Helmholtz resonator
    m.def(
        "helmholtz_frequency",
        &helmholtz_frequency,
        py::arg("V"),
        py::arg("A_neck"),
        py::arg("L_neck"),
        py::arg("c"),
        py::arg("end_correction") = 0.85,
        "Helmholtz resonator frequency.\n\n"
        "f = (c / 2π) * √(A / (V * L_eff))\n\n"
        "V              : cavity volume [m³]\n"
        "A_neck         : neck area [m²]\n"
        "L_neck         : neck length [m]\n"
        "c              : speed of sound [m/s]\n"
        "end_correction : end correction factor (default: 0.85 flanged)"
    );

    // Strouhal
    m.def(
        "strouhal",
        &strouhal,
        py::arg("f"),
        py::arg("L"),
        py::arg("u"),
        "Strouhal number: St = f * L / u"
    );

    m.def(
        "frequency_from_strouhal",
        &frequency_from_strouhal,
        py::arg("St"),
        py::arg("L"),
        py::arg("u"),
        "Frequency from Strouhal number: f = St * u / L"
    );

    // Convenience
    m.def(
        "quarter_wave_frequency",
        &quarter_wave_frequency,
        py::arg("L"),
        py::arg("c"),
        "Quarter-wave resonator fundamental: f = c / (4L)"
    );

    m.def(
        "half_wave_frequency",
        &half_wave_frequency,
        py::arg("L"),
        py::arg("c"),
        "Half-wave resonator fundamental: f = c / (2L)"
    );

    // Viscothermal boundary layers
    m.def(
        "stokes_layer",
        &stokes_layer,
        py::arg("nu"),
        py::arg("f"),
        "Stokes (viscous) boundary layer thickness: δ_ν = √(2ν/ω)\n\n"
        "nu : kinematic viscosity [m²/s]\n"
        "f  : frequency [Hz]"
    );

    m.def(
        "thermal_layer",
        &thermal_layer,
        py::arg("alpha"),
        py::arg("f"),
        "Thermal penetration depth: δ_κ = √(2α/ω)\n\n"
        "alpha : thermal diffusivity [m²/s]\n"
        "f     : frequency [Hz]"
    );

    m.def(
        "effective_viscothermal_layer",
        &effective_viscothermal_layer,
        py::arg("delta_nu"),
        py::arg("delta_kappa"),
        py::arg("gamma"),
        "Effective viscothermal layer: δ_eff = δ_ν + (γ-1)·δ_κ"
    );

    // Quality factor screening
    m.def(
        "helmholtz_Q",
        &helmholtz_Q,
        py::arg("V"),
        py::arg("A_neck"),
        py::arg("L_neck"),
        py::arg("nu"),
        py::arg("alpha"),
        py::arg("gamma"),
        py::arg("f"),
        "Helmholtz resonator Q factor (viscothermal losses).\n\n"
        "Screening estimate - expect ±50% accuracy.\n"
        "Real designs require experimental validation."
    );

    m.def(
        "tube_Q",
        &tube_Q,
        py::arg("L"),
        py::arg("D"),
        py::arg("nu"),
        py::arg("alpha"),
        py::arg("gamma"),
        py::arg("f"),
        "Quarter/half-wave tube Q factor (viscothermal losses).\n\n"
        "Screening estimate - expect ±50% accuracy."
    );

    m.def(
        "damping_ratio",
        &damping_ratio,
        py::arg("Q"),
        "Damping ratio from quality factor: ζ = 1/(2Q)"
    );

    m.def(
        "bandwidth",
        &bandwidth,
        py::arg("f0"),
        py::arg("Q"),
        "Half-power bandwidth: Δf = f₀/Q [Hz]"
    );

    // AcousticProperties struct binding - Bundle of acoustic properties
    py::class_<AcousticProperties>(m, "AcousticProperties",
        "Bundle of acoustic properties computed from (f, rho, c, p_rms).\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  wavelength        : Acoustic wavelength [m]\n"
        "  frequency         : Frequency [Hz]\n"
        "  impedance         : Characteristic acoustic impedance [Pa·s/m]\n"
        "  particle_velocity : Particle velocity amplitude [m/s]\n"
        "  spl               : Sound pressure level [dB]")
        .def_readonly("wavelength", &AcousticProperties::wavelength, "Acoustic wavelength [m]")
        .def_readonly("frequency", &AcousticProperties::frequency, "Frequency [Hz]")
        .def_readonly("impedance", &AcousticProperties::impedance, "Characteristic acoustic impedance [Pa·s/m]")
        .def_readonly("particle_velocity", &AcousticProperties::particle_velocity, "Particle velocity amplitude [m/s]")
        .def_readonly("spl", &AcousticProperties::spl, "Sound pressure level [dB]")
        .def("__repr__", [](const AcousticProperties& p) {
            return "<AcousticProperties: f=" + std::to_string(p.frequency) + " Hz, "
                   "λ=" + std::to_string(p.wavelength) + " m, "
                   "SPL=" + std::to_string(p.spl) + " dB>";
        });

    m.def(
        "acoustic_properties",
        &acoustic_properties,
        py::arg("f"),
        py::arg("rho"),
        py::arg("c"),
        py::arg("p_rms") = 20e-6,
        py::arg("p_ref") = 20e-6,
        "Compute all acoustic properties at once.\n\n"
        "Convenience function that computes all acoustic properties in a\n"
        "single call. Returns AcousticProperties struct with read-only attributes\n"
        "for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  f     : frequency [Hz]\n"
        "  rho   : density [kg/m³]\n"
        "  c     : speed of sound [m/s]\n"
        "  p_rms : RMS pressure amplitude [Pa] (default: 20e-6 Pa = 0 dB reference)\n"
        "  p_ref : reference pressure for SPL [Pa] (default: 20e-6 Pa in air)\n\n"
        "Returns: AcousticProperties object with attributes:\n"
        "  wavelength, frequency, impedance, particle_velocity, spl\n\n"
        "Example:\n"
        "  >>> props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=1.0)\n"
        "  >>> print(f'Wavelength: {props.wavelength:.3f} m')\n"
        "  >>> print(f'Impedance: {props.impedance:.1f} Pa·s/m')\n"
        "  >>> print(f'SPL: {props.spl:.1f} dB')"
    );


    // =========================================================================
    // Transfer Matrix Method - Advanced Thermoacoustics
    // =========================================================================

    py::class_<TransferMatrix>(m, "TransferMatrix")
        .def(py::init<>())
        .def_readwrite("T11", &TransferMatrix::T11)
        .def_readwrite("T12", &TransferMatrix::T12)
        .def_readwrite("T21", &TransferMatrix::T21)
        .def_readwrite("T22", &TransferMatrix::T22)
        .def("__mul__", &TransferMatrix::operator*);

    py::class_<AcousticMedium>(m, "AcousticMedium")
        .def(py::init<>())
        .def_readwrite("rho", &AcousticMedium::rho)
        .def_readwrite("c", &AcousticMedium::c);

    py::class_<LinerOrificeGeometry>(m, "LinerOrificeGeometry")
        .def(py::init<>())
        .def_readwrite("d_orifice", &LinerOrificeGeometry::d_orifice)
        .def_readwrite("l_orifice", &LinerOrificeGeometry::l_orifice)
        .def_readwrite("porosity", &LinerOrificeGeometry::porosity)
        .def_readwrite("Cd", &LinerOrificeGeometry::Cd);

    py::class_<LinerFlowState>(m, "LinerFlowState")
        .def(py::init<>())
        .def_readwrite("u_bias", &LinerFlowState::u_bias)
        .def_readwrite("u_grazing", &LinerFlowState::u_grazing);

    py::class_<LinerCavity>(m, "LinerCavity")
        .def(py::init<>())
        .def_readwrite("depth", &LinerCavity::depth);

    m.def(
        "absorption_from_impedance_norm",
        &absorption_from_impedance_norm,
        py::arg("z_norm")
    );

    m.def(
        "liner_sdof_impedance_norm",
        &liner_sdof_impedance_norm,
        py::arg("freq"),
        py::arg("orifice"),
        py::arg("cavity"),
        py::arg("flow"),
        py::arg("medium")
    );

    m.def(
        "liner_sdof_absorption",
        &liner_sdof_absorption,
        py::arg("freq"),
        py::arg("orifice"),
        py::arg("cavity"),
        py::arg("flow"),
        py::arg("medium")
    );

    m.def(
        "sweep_liner_sdof_absorption",
        &sweep_liner_sdof_absorption,
        py::arg("freqs"),
        py::arg("orifice"),
        py::arg("cavity"),
        py::arg("flow"),
        py::arg("medium")
    );

    m.def(
        "liner_2dof_serial_impedance_norm",
        &liner_2dof_serial_impedance_norm,
        py::arg("freq"),
        py::arg("face_orifice"),
        py::arg("septum_orifice"),
        py::arg("depth_1"),
        py::arg("depth_2"),
        py::arg("face_flow"),
        py::arg("septum_flow"),
        py::arg("medium")
    );

    m.def(
        "liner_2dof_serial_absorption",
        &liner_2dof_serial_absorption,
        py::arg("freq"),
        py::arg("face_orifice"),
        py::arg("septum_orifice"),
        py::arg("depth_1"),
        py::arg("depth_2"),
        py::arg("face_flow"),
        py::arg("septum_flow"),
        py::arg("medium")
    );

    m.def(
        "sweep_liner_2dof_serial_absorption",
        &sweep_liner_2dof_serial_absorption,
        py::arg("freqs"),
        py::arg("face_orifice"),
        py::arg("septum_orifice"),
        py::arg("depth_1"),
        py::arg("depth_2"),
        py::arg("face_flow"),
        py::arg("septum_flow"),
        py::arg("medium")
    );

    m.def("orifice_impedance_with_flow", &orifice_impedance_with_flow,
        py::arg("freq"), py::arg("u_bias"), py::arg("u_grazing"),
        py::arg("d_orifice"), py::arg("l_orifice"), py::arg("porosity"),
        py::arg("Cd"), py::arg("rho"), py::arg("c"));

    m.def("quarter_wave_resonator_tmm", &quarter_wave_resonator_tmm,
        py::arg("freq"), py::arg("L_tube"), py::arg("A_duct"), py::arg("A_tube"),
        py::arg("d_orifice"), py::arg("l_orifice"), py::arg("porosity"), py::arg("Cd"),
        py::arg("u_bias"), py::arg("u_grazing"), py::arg("rho"), py::arg("c"));

    m.def("is_whistling_risk", &is_whistling_risk,
        py::arg("freq"), py::arg("u_bias"), py::arg("d_orifice"));

    // =========================================================================
    // Can-Annular Combustor Acoustics (Bloch-Floquet Theory)
    // =========================================================================

    py::class_<CanAnnularGeometry>(m, "CanAnnularGeometry")
        .def(py::init<>())
        .def_readwrite("n_cans", &CanAnnularGeometry::n_cans)
        .def_readwrite("length_can", &CanAnnularGeometry::length_can)
        .def_readwrite("area_can", &CanAnnularGeometry::area_can)
        .def_readwrite("radius_plenum", &CanAnnularGeometry::radius_plenum)
        .def_readwrite("area_plenum", &CanAnnularGeometry::area_plenum)
        .def("__repr__", [](const CanAnnularGeometry& g) {
            return "CanAnnularGeometry(n_cans=" + std::to_string(g.n_cans) +
                   ", length_can=" + std::to_string(g.length_can) +
                   ", area_can=" + std::to_string(g.area_can) +
                   ", radius_plenum=" + std::to_string(g.radius_plenum) +
                   ", area_plenum=" + std::to_string(g.area_plenum) + ")";
        });

    py::class_<BlochMode>(m, "BlochMode")
        .def(py::init<>())
        .def_readwrite("m_azimuthal", &BlochMode::m_azimuthal)
        .def_readwrite("frequency", &BlochMode::frequency)
        .def_readwrite("n_cans", &BlochMode::n_cans)
        .def("symmetry_type", &BlochMode::symmetry_type)
        .def("__repr__", [](const BlochMode& m) {
            return "BlochMode(m=" + std::to_string(m.m_azimuthal) +
                   ", f=" + std::to_string(m.frequency) + " Hz, " +
                   m.symmetry_type() + ")";
        });

    m.def("can_annular_eigenmodes", &can_annular_eigenmodes,
        py::arg("geom"),
        py::arg("c_can"),
        py::arg("c_plenum"),
        py::arg("rho_can"),
        py::arg("rho_plenum"),
        py::arg("f_max") = 2000.0,
        py::arg("bc_can_top") = BoundaryCondition::Closed,
        "Find all can-annular combustor eigenmodes using Bloch-Floquet theory.\n\n"
        "Solves the dispersion relation D(omega) = Y_can + Y_annulus = 0\n"
        "for a system of N cans coupled by an annular plenum.\n\n"
        "Uses magnitude minimization + Muller refinement (fast, reliable).\n"
        "For research-grade Argument Principle, use can_annular_eigenmodes_with_method.\n\n"
        "Parameters:\n"
        "  geom       : CanAnnularGeometry with N cans, lengths, areas\n"
        "  c_can      : speed of sound in can [m/s]\n"
        "  c_plenum   : speed of sound in plenum [m/s] (often different!)\n"
        "  rho_can    : density in can [kg/m^3]\n"
        "  rho_plenum : density in plenum [kg/m^3]\n"
        "  f_max      : maximum frequency to search [Hz] (default: 2000)\n"
        "  bc_can_top : boundary condition at can top (default: Closed)\n\n"
        "Returns: List of BlochMode sorted by frequency\n\n"
        "Method: Continuous sliding window scan + magnitude minimization\n"
        "Accuracy: Finds modes wherever they occur (no harmonic assumptions)\n\n"
        "References:\n"
        "  - Noiray et al. (2011): Unified framework for combustion instability\n"
        "  - Evesque & Polifke (2005): Low-order acoustic modelling\n\n"
        "Example:\n"
        "  >>> geom = cb.CanAnnularGeometry()\n"
        "  >>> geom.n_cans = 24\n"
        "  >>> geom.length_can = 0.5  # m\n"
        "  >>> geom.area_can = 0.01   # m^2\n"
        "  >>> geom.radius_plenum = 0.5  # m\n"
        "  >>> geom.area_plenum = 0.02   # m^2\n"
        "  >>> modes = cb.can_annular_eigenmodes(geom, c_can=500, c_plenum=600,\n"
        "  ...                                    rho_can=1.0, rho_plenum=1.2)\n"
        "  >>> for mode in modes[:5]:\n"
        "  ...     print(f'm={mode.m_azimuthal}, f={mode.frequency:.1f} Hz')"
    );

    // Adapter function to convert CanAnnularFlowGeometry to CanAnnularGeometry
    m.def("to_acoustic_geometry", &to_acoustic_geometry,
        py::arg("flow_geom"),
        py::arg("n_cans"),
        py::arg("L_can"),
        py::arg("D_can"),
        "Convert CanAnnularFlowGeometry to CanAnnularGeometry for acoustic solvers.\n\n"
        "This adapter bridges flow geometry (with primary/transition/annular sections)\n"
        "to acoustic geometry (cans + annular plenum).\n\n"
        "Parameters:\n"
        "  flow_geom : CanAnnularFlowGeometry from geometry.h\n"
        "  n_cans    : Number of cans (required for acoustic model)\n"
        "  L_can     : Can length [m] (acoustic model uses single section)\n"
        "  D_can     : Can diameter [m] (for circular cross-section)\n\n"
        "Returns: CanAnnularGeometry for can_annular_eigenmodes()\n\n"
        "Notes:\n"
        "  - Plenum radius derived from flow_geom.D_mean() / 2\n"
        "  - Plenum area from flow_geom.area()\n"
        "  - Can area computed as π*(D_can/2)²\n\n"
        "Example:\n"
        "  >>> flow_geom = cb.CanAnnularFlowGeometry(L=0.3, D_inner=0.4, D_outer=0.6,\n"
        "  ...                                        L_primary=0.2, D_primary=0.1,\n"
        "  ...                                        L_transition=0.05)\n"
        "  >>> geom = cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=0.5, D_can=0.1)\n"
        "  >>> modes = cb.can_annular_eigenmodes(geom, c_can=500, c_plenum=550,\n"
        "  ...                                    rho_can=1.2, rho_plenum=1.0)"
    );

    // Annular Duct Geometry and Modes
    // AnnularDuctGeometry is now an alias for Annulus - keep for backward compatibility
    using combaero::Annulus;
    using combaero::AnnularMode;
    using combaero::annular_duct_eigenmodes;
    using combaero::annular_duct_modes_analytical;

    // AnnularDuctGeometry is deprecated, use Annulus instead
    // This binding provides backward compatibility with the old API
    py::class_<Annulus>(m, "AnnularDuctGeometry",
        "Annular duct geometry for acoustic solvers.\n\n"
        "Deprecated: Use Annulus from geometry.h instead.\n\n"
        "Note: AnnularDuctGeometry is now an alias for Annulus with acoustic solver compatibility.")
        .def(py::init<>())
        .def_property("length",
            [](const Annulus& g) { return g.L; },
            [](Annulus& g, double val) { g.L = val; },
            "Axial length [m]")
        .def_property("radius_inner",
            [](const Annulus& g) { return g.radius_inner(); },
            [](Annulus& g, double val) { g.D_inner = val * 2.0; },
            "Inner radius [m]")
        .def_property("radius_outer",
            [](const Annulus& g) { return g.radius_outer(); },
            [](Annulus& g, double val) { g.D_outer = val * 2.0; },
            "Outer radius [m]")
        .def_readwrite("n_azimuthal_max", &Annulus::n_azimuthal_max,
            "Maximum azimuthal mode number to search")
        .def("area", &Annulus::area, "Cross-sectional area [m²]")
        .def("__repr__", [](const Annulus& g) {
            return "AnnularDuctGeometry(length=" + std::to_string(g.L) +
                   ", r_inner=" + std::to_string(g.radius_inner()) +
                   ", r_outer=" + std::to_string(g.radius_outer()) +
                   ", n_azimuthal_max=" + std::to_string(g.n_azimuthal_max) + ")";
        });

    py::class_<AnnularMode>(m, "AnnularMode")
        .def(py::init<>())
        .def_readwrite("m_azimuthal", &AnnularMode::m_azimuthal)
        .def_readwrite("n_axial", &AnnularMode::n_axial)
        .def_readwrite("frequency", &AnnularMode::frequency)
        .def("mode_type", &AnnularMode::mode_type)
        .def("__repr__", [](const AnnularMode& mode) {
            return "AnnularMode(m=" + std::to_string(mode.m_azimuthal) +
                   ", n=" + std::to_string(mode.n_axial) +
                   ", f=" + std::to_string(mode.frequency) + " Hz, " +
                   mode.mode_type() + ")";
        });

    m.def("annular_duct_eigenmodes", &annular_duct_eigenmodes,
        py::arg("geom"),
        py::arg("c"),
        py::arg("rho"),
        py::arg("f_max"),
        py::arg("bc_ends") = BoundaryCondition::Closed,
        "Find all annular duct eigenmodes using Argument Principle.\n\n"
        "Pure annular waveguide (no cans) with Bloch-Floquet periodic BCs.\n"
        "Uses Nyquist contour integration for robust root finding.\n\n"
        "Parameters:\n"
        "  geom    : Annulus geometry (from geometry.h)\n"
        "  c       : speed of sound [m/s]\n"
        "  rho     : density [kg/m^3]\n"
        "  f_max   : maximum frequency [Hz]\n"
        "  bc_ends : boundary condition at duct ends (default: Closed)\n\n"
        "Returns: List of AnnularMode sorted by frequency\n\n"
        "Method: Argument Principle (gold standard for finding ALL roots)\n"
        "Guarantees: No missed modes, handles complex frequencies\n\n"
        "Example:\n"
        "  >>> geom = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)\n"
        "  >>> geom.n_azimuthal_max = 5\n"
        "  >>> modes = cb.annular_duct_eigenmodes(geom, c=500, rho=1.2, f_max=1000)\n"
        "  >>> for mode in modes[:5]:\n"
        "  ...     print(f'm={mode.m_azimuthal}, n={mode.n_axial}, f={mode.frequency:.1f} Hz')"
    );

    m.def("annular_duct_modes_analytical", &annular_duct_modes_analytical,
        py::arg("geom"),
        py::arg("c"),
        py::arg("f_max"),
        py::arg("bc_ends") = BoundaryCondition::Closed
    );

    // Utility functions
    m.def(
        "wavelength",
        &wavelength,
        py::arg("f"),
        py::arg("c"),
        "Acoustic wavelength from frequency.\n\n"
        "λ = c / f\n\n"
        "Parameters:\n"
        "  f : frequency [Hz]\n"
        "  c : speed of sound [m/s]\n\n"
        "Returns: wavelength [m]"
    );

    m.def(
        "frequency_from_wavelength",
        &frequency_from_wavelength,
        py::arg("lambda"),
        py::arg("c"),
        "Frequency from wavelength.\n\n"
        "f = c / λ\n\n"
        "Parameters:\n"
        "  lambda : wavelength [m]\n"
        "  c      : speed of sound [m/s]\n\n"
        "Returns: frequency [Hz]"
    );

    m.def(
        "acoustic_impedance",
        &acoustic_impedance,
        py::arg("rho"),
        py::arg("c"),
        "Characteristic acoustic impedance.\n\n"
        "Z = ρ · c\n\n"
        "Parameters:\n"
        "  rho : density [kg/m³]\n"
        "  c   : speed of sound [m/s]\n\n"
        "Returns: impedance [Pa·s/m]"
    );

    m.def(
        "sound_pressure_level",
        &sound_pressure_level,
        py::arg("p_rms"),
        py::arg("p_ref") = 20e-6,
        "Sound pressure level in dB.\n\n"
        "SPL = 20 · log₁₀(p_rms / p_ref)\n\n"
        "Parameters:\n"
        "  p_rms : RMS pressure amplitude [Pa]\n"
        "  p_ref : reference pressure [Pa] (default: 20e-6 Pa in air)\n\n"
        "Returns: SPL [dB]"
    );

    m.def(
        "particle_velocity",
        &particle_velocity,
        py::arg("p"),
        py::arg("rho"),
        py::arg("c"),
        "Particle velocity from pressure amplitude.\n\n"
        "u = p / (ρ · c)\n\n"
        "Parameters:\n"
        "  p   : pressure amplitude [Pa]\n"
        "  rho : density [kg/m³]\n"
        "  c   : speed of sound [m/s]\n\n"
        "Returns: particle velocity [m/s]"
    );

    // Residence time
    m.def(
        "residence_time",
        py::overload_cast<double, double>(&residence_time),
        py::arg("V"),
        py::arg("Q"),
        "Residence time τ = V / Q̇ [s].\n\n"
        "V : volume [m³]\n"
        "Q : volumetric flow rate [m³/s]"
    );

    m.def(
        "residence_time_tube",
        py::overload_cast<const Tube&, double>(&residence_time),
        py::arg("tube"),
        py::arg("Q"),
        "Residence time for a tube [s]."
    );

    m.def(
        "residence_time_annulus",
        py::overload_cast<const Annulus&, double>(&residence_time),
        py::arg("annulus"),
        py::arg("Q"),
        "Residence time for an annulus [s]."
    );

    m.def(
        "residence_time_can_annular",
        py::overload_cast<const CanAnnularFlowGeometry&, double>(&residence_time),
        py::arg("geom"),
        py::arg("Q"),
        "Residence time for CanAnnularFlowGeometry based on total volume [s]."
    );

    m.def(
        "residence_time_mdot",
        py::overload_cast<double, double, double>(&residence_time_mdot),
        py::arg("V"),
        py::arg("mdot"),
        py::arg("rho"),
        "Residence time from mass flow: τ = V·ρ / ṁ [s].\n\n"
        "V    : volume [m³]\n"
        "mdot : mass flow rate [kg/s]\n"
        "rho  : density [kg/m³]"
    );

    m.def(
        "residence_time_mdot_can_annular",
        py::overload_cast<const CanAnnularFlowGeometry&, double, double>(&residence_time_mdot),
        py::arg("geom"),
        py::arg("mdot"),
        py::arg("rho"),
        "Residence time for CanAnnularFlowGeometry from mass flow [s]."
    );

    m.def(
        "space_velocity",
        &space_velocity,
        py::arg("Q"),
        py::arg("V"),
        "Space velocity SV = Q̇ / V = 1/τ [1/s]."
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

    m.def(
        "solve_orifice_mdot",
        &solve_orifice_mdot,
        py::arg("geom"),
        py::arg("dP"),
        py::arg("rho"),
        py::arg("mu"),
        py::arg("P_upstream") = 101325.0,
        py::arg("kappa") = 0.0,
        py::arg("correlation") = CdCorrelation::ReaderHarrisGallagher,
        py::arg("tol") = 1e-6,
        py::arg("max_iter") = 20,
        "Iterative solver for orifice mass flow with Cd-Re coupling.\n\n"
        "Solves the coupled system:\n"
        "  mdot = ε · Cd(Re_D) · A · √(2 · ρ · ΔP)\n"
        "  Re_D = 4 · mdot / (π · D · μ)\n\n"
        "Iterates until mdot converges (typically 3-5 iterations).\n\n"
        "Parameters:\n"
        "  geom        : OrificeGeometry\n"
        "  dP          : differential pressure [Pa]\n"
        "  rho         : upstream density [kg/m³]\n"
        "  mu          : dynamic viscosity [Pa·s]\n"
        "  P_upstream  : absolute upstream pressure [Pa] (default: 101325)\n"
        "  kappa       : isentropic exponent cp/cv (0 = incompressible)\n"
        "  correlation : Cd correlation (default: ReaderHarrisGallagher)\n"
        "  tol         : relative convergence tolerance (default: 1e-6)\n"
        "  max_iter    : maximum iterations (default: 20)\n\n"
        "Returns: converged mass flow rate [kg/s]\n\n"
        "Raises: RuntimeError if fails to converge\n\n"
        "Example:\n"
        "  >>> geom = OrificeGeometry(d=0.05, D=0.1)\n"
        "  >>> mdot = solve_orifice_mdot(geom, 10000, 1.2, 1.8e-5)\n"
        "  >>> print(f'Mass flow: {mdot:.4f} kg/s')"
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
        "Solve for discharge coefficient from measurement.\n\n"
        "Cd = mdot / (A * sqrt(2 * rho * dP))\n\n"
        "Parameters:\n"
        "  geom : OrificeGeometry\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  dP   : pressure drop [Pa]\n"
        "  rho  : density [kg/m^3]\n\n"
        "Returns: discharge coefficient Cd [-]"
    );

    m.def(
        "expansibility_factor",
        &expansibility_factor,
        py::arg("beta"),
        py::arg("dP"),
        py::arg("P_upstream"),
        py::arg("kappa"),
        "Expansibility factor for compressible gas flow through orifices.\n\n"
        "Accounts for gas expansion as it accelerates through the orifice.\n"
        "Formula from ISO 5167-2:2003, Section 5.3.2.2.\n\n"
        "ε = 1 - (0.351 + 0.256·β⁴ + 0.93·β⁸) · [1 - (1 - τ)^(1/κ)]\n\n"
        "where τ = ΔP/P_upstream, β = d/D, κ = cp/cv\n\n"
        "Usage: mdot = ε · Cd · A · √(2 · ρ_upstream · ΔP)\n\n"
        "Parameters:\n"
        "  beta        : diameter ratio d/D [-]\n"
        "  dP          : differential pressure [Pa]\n"
        "  P_upstream  : absolute upstream pressure [Pa]\n"
        "  kappa       : isentropic exponent cp/cv [-]\n\n"
        "Returns: expansibility factor ε [-] (0 < ε ≤ 1)\n\n"
        "Valid for: 0.1 ≤ β ≤ 0.75, τ ≤ 0.25, ideal gas"
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
        py::arg("Re_d"),
        "Thickness correction factor for thick-plate orifices.\n\n"
        "Idelchik model: Cd rises (reattachment) then falls (friction).\n\n"
        "Multiplies thin-plate Cd to account for flow reattachment.\n"
        "Returns 1.0 for t/d <= 0.02 (thin plate)."
    );

    // OrificeFlowResult struct binding - Bundle of orifice flow properties
    py::class_<OrificeFlowResult>(m, "OrificeFlowResult",
        "Bundle of orifice flow properties computed from (geom, dP, T, P, mu, Z).\n\n"
        "All properties are read-only attributes computed in a single call.\n"
        "Provides IDE autocomplete and type safety.\n\n"
        "Attributes:\n"
        "  mdot          : Mass flow rate [kg/s]\n"
        "  v             : Velocity through orifice throat [m/s]\n"
        "  Re_D          : Pipe Reynolds number (based on D) [-]\n"
        "  Re_d          : Orifice Reynolds number (based on d) [-]\n"
        "  Cd            : Discharge coefficient [-]\n"
        "  epsilon       : Expansibility factor [-] (1.0 for incompressible)\n"
        "  rho_corrected : Density corrected for compressibility [kg/m³]")
        .def_readonly("mdot", &OrificeFlowResult::mdot, "Mass flow rate [kg/s]")
        .def_readonly("v", &OrificeFlowResult::v, "Velocity through orifice throat [m/s]")
        .def_readonly("Re_D", &OrificeFlowResult::Re_D, "Pipe Reynolds number (based on D) [-]")
        .def_readonly("Re_d", &OrificeFlowResult::Re_d, "Orifice Reynolds number (based on d) [-]")
        .def_readonly("Cd", &OrificeFlowResult::Cd, "Discharge coefficient [-]")
        .def_readonly("epsilon", &OrificeFlowResult::epsilon, "Expansibility factor [-]")
        .def_readonly("rho_corrected", &OrificeFlowResult::rho_corrected, "Density corrected for compressibility [kg/m³]")
        .def("__repr__", [](const OrificeFlowResult& r) {
            return "<OrificeFlowResult: mdot=" + std::to_string(r.mdot) + " kg/s, "
                   "v=" + std::to_string(r.v) + " m/s, "
                   "Cd=" + std::to_string(r.Cd) + ", "
                   "Re_d=" + std::to_string(r.Re_d) + ">";
        });

    m.def(
        "orifice_flow",
        &orifice_flow,
        py::arg("geom"),
        py::arg("dP"),
        py::arg("T"),
        py::arg("P"),
        py::arg("mu"),
        py::arg("Z") = 1.0,
        py::arg("X") = std::vector<double>(),
        py::arg("kappa") = 0.0,
        py::arg("correlation") = CdCorrelation::ReaderHarrisGallagher,
        "Compute all orifice flow properties at once with real gas correction.\n\n"
        "Convenience function that computes all orifice flow properties in a\n"
        "single call. Returns OrificeFlowResult struct with read-only attributes\n"
        "for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  geom        : OrificeGeometry\n"
        "  dP          : differential pressure [Pa]\n"
        "  T           : temperature [K]\n"
        "  P           : absolute upstream pressure [Pa]\n"
        "  mu          : dynamic viscosity [Pa·s]\n"
        "  Z           : compressibility factor [-] (default: 1.0 = ideal gas)\n"
        "  X           : mole fractions [-] (optional, uses air if empty)\n"
        "  kappa       : isentropic exponent cp/cv [-] (default: 0.0 = incompressible)\n"
        "  correlation : Cd correlation (default: ReaderHarrisGallagher)\n\n"
        "Returns: OrificeFlowResult object with attributes:\n"
        "  mdot, v, Re_D, Re_d, Cd, epsilon, rho_corrected\n\n"
        "Real Gas Support:\n"
        "  The compressibility factor Z corrects density for non-ideal behavior:\n"
        "    rho_corrected = rho_ideal / Z\n"
        "  For ideal gas: Z = 1.0 (default)\n"
        "  For real gas: Z typically ranges 0.7-1.0 depending on P, T, composition\n"
        "  Users can obtain Z from REFPROP, Peng-Robinson, or other EOS libraries\n\n"
        "Example:\n"
        "  >>> geom = cb.OrificeGeometry()\n"
        "  >>> geom.d = 0.05  # 50 mm orifice\n"
        "  >>> geom.D = 0.1   # 100 mm pipe\n"
        "  >>> result = cb.orifice_flow(geom, dP=50000, T=300, P=10e6, mu=1.1e-5, Z=0.85)\n"
        "  >>> print(f'Mass flow: {result.mdot:.3f} kg/s')\n"
        "  >>> print(f'Velocity: {result.v:.2f} m/s')\n"
        "  >>> print(f'Cd: {result.Cd:.3f}')"
    );

    // Utility functions
    m.def(
        "orifice_velocity_from_mdot",
        &orifice_velocity_from_mdot,
        py::arg("mdot"),
        py::arg("rho"),
        py::arg("d"),
        py::arg("Z") = 1.0,
        "Velocity through orifice from mass flow rate with real gas correction.\n\n"
        "v = mdot / (rho_corrected * A)\n"
        "where rho_corrected = rho / Z, A = π·d²/4\n\n"
        "Parameters:\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  rho  : density [kg/m³]\n"
        "  d    : orifice diameter [m]\n"
        "  Z    : compressibility factor [-] (default: 1.0 = ideal gas)\n\n"
        "Returns: velocity [m/s]"
    );

    m.def(
        "orifice_area_from_beta",
        &orifice_area_from_beta,
        py::arg("D"),
        py::arg("beta"),
        "Orifice area from beta ratio.\n\n"
        "A = π·(D·beta/2)²\n\n"
        "Parameters:\n"
        "  D    : pipe diameter [m]\n"
        "  beta : diameter ratio d/D [-]\n\n"
        "Returns: orifice area [m²]"
    );

    m.def(
        "beta_from_diameters",
        &beta_from_diameters,
        py::arg("d"),
        py::arg("D"),
        "Beta ratio from diameters.\n\n"
        "beta = d / D\n\n"
        "Parameters:\n"
        "  d : orifice diameter [m]\n"
        "  D : pipe diameter [m]\n\n"
        "Returns: beta ratio [-]"
    );

    m.def(
        "orifice_Re_d_from_mdot",
        &orifice_Re_d_from_mdot,
        py::arg("mdot"),
        py::arg("d"),
        py::arg("mu"),
        "Orifice Reynolds number from mass flow rate.\n\n"
        "Re_d = 4·mdot / (π·d·μ)\n\n"
        "Parameters:\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  d    : orifice diameter [m]\n"
        "  mu   : dynamic viscosity [Pa·s]\n\n"
        "Returns: Reynolds number [-]"
    );

    // =========================================================================
    // Pipe Roughness Database
    // =========================================================================

    m.def(
        "pipe_roughness",
        &pipe_roughness,
        py::arg("material"),
        "Get absolute roughness for a standard pipe material.\n\n"
        "Returns the absolute roughness ε [m] for common pipe materials.\n"
        "Material names are case-insensitive.\n\n"
        "Parameters:\n"
        "  material : pipe material name (str)\n\n"
        "Returns: absolute roughness ε [m]\n\n"
        "Available materials:\n"
        "  - 'smooth', 'drawn_tubing', 'pvc', 'plastic'\n"
        "  - 'commercial_steel', 'new_steel', 'wrought_iron'\n"
        "  - 'galvanized_iron', 'galvanized_steel', 'rusted_steel'\n"
        "  - 'cast_iron', 'asphalted_cast_iron'\n"
        "  - 'concrete', 'rough_concrete'\n"
        "  - 'riveted_steel', 'wood_stave', 'corrugated_metal'\n\n"
        "Raises: ValueError if material not found\n\n"
        "References:\n"
        "  - Moody (1944), Colebrook (1939), White (2011),\n"
        "    Munson (2013), Crane (2009)\n\n"
        "Example:\n"
        "  >>> eps = pipe_roughness('commercial_steel')\n"
        "  >>> print(f'Roughness: {eps*1e6:.1f} μm')  # 45.0 μm"
    );

    m.def(
        "standard_pipe_roughness",
        &standard_pipe_roughness,
        "Get all standard pipe roughness values.\n\n"
        "Returns a dictionary mapping material names to absolute\n"
        "roughness values [m].\n\n"
        "Returns: dict[str, float] - material -> roughness [m]\n\n"
        "Example:\n"
        "  >>> roughness_db = standard_pipe_roughness()\n"
        "  >>> for material, eps in sorted(roughness_db.items()):\n"
        "  ...     print(f'{material:20s}: {eps*1e6:8.2f} μm')"
    );

    // =========================================================================
    // Materials Database
    // =========================================================================

    // Material thermal conductivity functions
    m.def(
        "k_inconel718",
        &combaero::materials::k_inconel718,
        py::arg("T"),
        "Thermal conductivity of Inconel 718 [W/(m·K)].\n\n"
        "Ni-based superalloy commonly used in turbine applications.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 300-1200 K\n"
        "Source: Haynes International, Special Metals Corporation"
    );

    m.def(
        "k_haynes230",
        &combaero::materials::k_haynes230,
        py::arg("T"),
        "Thermal conductivity of Haynes 230 [W/(m·K)].\n\n"
        "Ni-Cr-W-Mo alloy with high-temperature oxidation resistance.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 300-1400 K\n"
        "Source: Haynes International Technical Data"
    );

    m.def(
        "k_stainless_steel_316",
        &combaero::materials::k_stainless_steel_316,
        py::arg("T"),
        "Thermal conductivity of Stainless Steel 316 [W/(m·K)].\n\n"
        "Austenitic stainless steel with corrosion resistance.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 300-1200 K\n"
        "Source: NIST, ASM Handbook"
    );

    m.def(
        "k_aluminum_6061",
        &combaero::materials::k_aluminum_6061,
        py::arg("T"),
        "Thermal conductivity of Aluminum 6061 [W/(m·K)].\n\n"
        "Aluminum alloy, limited to <300°C due to T6 temper stability.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 200-600 K\n"
        "Source: ASM Handbook, Aluminum Association"
    );

    m.def(
        "k_tbc_ysz",
        &combaero::materials::k_tbc_ysz,
        py::arg("T"),
        py::arg("hours") = 0.0,
        py::arg("is_ebpvd") = false,
        "Thermal conductivity of YSZ thermal barrier coating [W/(m·K)].\n\n"
        "7-8 wt% Y2O3 stabilized zirconia with sintering model.\n"
        "Conductivity increases over time at high temperature (>1073 K) due to pore closure.\n\n"
        "Parameters:\n"
        "  T        : temperature [K]\n"
        "  hours    : operating hours at temperature (default: 0 = as-sprayed)\n"
        "  is_ebpvd : True for EB-PVD coating, False for APS (default: False)\n\n"
        "Valid range: 300-1700 K\n"
        "Source: NASA TM-2010-216765 (Zhu/Miller sintering model)\n\n"
        "Coating types:\n"
        "  - APS (Atmospheric Plasma Spray): Lower initial k, splat boundaries\n"
        "  - EB-PVD (Electron Beam PVD): Higher initial k, columnar structure\n\n"
        "Example:\n"
        "  >>> k_aps_fresh = cb.k_tbc_ysz(T=1500)  # APS as-sprayed\n"
        "  >>> k_aps_aged = cb.k_tbc_ysz(T=1500, hours=1000)  # APS after 1000h\n"
        "  >>> k_ebpvd = cb.k_tbc_ysz(T=1500, is_ebpvd=True)  # EB-PVD as-sprayed"
    );

    m.def(
        "list_materials",
        &combaero::materials::list_materials,
        "List all available materials in the database.\n\n"
        "Returns: List of material names (strings)\n\n"
        "Example:\n"
        "  >>> materials = cb.list_materials()\n"
        "  >>> print(materials)\n"
        "  ['inconel718', 'haynes230', 'ss316', 'al6061', 'ysz']"
    );

    // =========================================================================
    // Advanced Cooling Correlations
    // =========================================================================

    // Rib-enhanced cooling
    m.def(
        "rib_enhancement_factor",
        &combaero::cooling::rib_enhancement_factor,
        py::arg("e_D"),
        py::arg("P_e"),
        py::arg("alpha"),
        "Rib enhancement factor for heat transfer [-].\n\n"
        "Accounts for increased surface area and turbulence from ribs.\n\n"
        "Parameters:\n"
        "  e_D   : rib height to hydraulic diameter ratio [-]\n"
        "  P_e   : rib pitch to rib height ratio [-]\n"
        "  alpha : rib angle [degrees]\n\n"
        "Valid range: e_D = 0.02-0.1, P_e = 5-20, alpha = 30-90 deg\n"
        "Source: Han et al. (1988)\n\n"
        "Returns: Enhancement factor [-] (typically 1.5-4.0)"
    );

    m.def(
        "rib_friction_multiplier",
        &combaero::cooling::rib_friction_multiplier,
        py::arg("e_D"),
        py::arg("P_e"),
        "Rib friction factor multiplier [-].\n\n"
        "Accounts for increased pressure drop from ribs.\n\n"
        "Parameters:\n"
        "  e_D : rib height to hydraulic diameter ratio [-]\n"
        "  P_e : rib pitch to rib height ratio [-]\n\n"
        "Valid range: e_D = 0.02-0.1, P_e = 5-20\n"
        "Source: Han et al. (1988)\n\n"
        "Returns: Friction multiplier [-] (typically 2-10)"
    );

    // Impingement cooling
    m.def(
        "impingement_nusselt",
        &combaero::cooling::impingement_nusselt,
        py::arg("Re_jet"),
        py::arg("Pr"),
        py::arg("z_D"),
        py::arg("x_D") = 0.0,
        py::arg("y_D") = 0.0,
        "Impingement jet Nusselt number correlation [-].\n\n"
        "For single jet or jet array impinging on flat plate.\n\n"
        "Parameters:\n"
        "  Re_jet : Reynolds number based on jet diameter [-]\n"
        "  Pr     : Prandtl number [-]\n"
        "  z_D    : jet-to-plate distance / jet diameter [-]\n"
        "  x_D    : streamwise spacing / jet diameter [-] (default: 0 for single jet)\n"
        "  y_D    : spanwise spacing / jet diameter [-] (default: 0 for single jet)\n\n"
        "Valid range: Re_jet = 5000-80000, z_D = 1-12, x_D/y_D = 4-16\n"
        "Source: Florschuetz et al. (1981), Martin (1977)\n\n"
        "Returns: Average Nusselt number [-]"
    );

    // Film cooling
    m.def(
        "film_cooling_effectiveness",
        &combaero::cooling::film_cooling_effectiveness,
        py::arg("x_D"),
        py::arg("M"),
        py::arg("DR"),
        py::arg("alpha_deg"),
        "Adiabatic film cooling effectiveness [-].\n\n"
        "Measures how well coolant film protects surface from hot gas.\n\n"
        "Parameters:\n"
        "  x_D       : downstream distance / hole diameter [-]\n"
        "  M         : blowing ratio (rho_c*v_c)/(rho_inf*v_inf) [-]\n"
        "  DR        : density ratio rho_c/rho_inf [-]\n"
        "  alpha_deg : injection angle [degrees]\n\n"
        "Valid range: M = 0.3-2.5, DR = 1.2-2.0, alpha = 20-90 deg, x_D >= 0\n"
        "Source: Baldauf et al. (2002)\n\n"
        "Returns: Adiabatic effectiveness eta [-] (0 = no cooling, 1 = perfect)"
    );

    m.def(
        "film_cooling_effectiveness_avg",
        &combaero::cooling::film_cooling_effectiveness_avg,
        py::arg("x_D"),
        py::arg("M"),
        py::arg("DR"),
        py::arg("alpha_deg"),
        py::arg("s_D") = 3.0,
        "Laterally averaged film cooling effectiveness [-].\n\n"
        "Averaged across span for design calculations.\n\n"
        "Parameters:\n"
        "  x_D       : downstream distance / hole diameter [-]\n"
        "  M         : blowing ratio [-]\n"
        "  DR        : density ratio [-]\n"
        "  alpha_deg : injection angle [degrees]\n"
        "  s_D       : hole spacing / hole diameter [-] (default: 3.0)\n\n"
        "Valid range: M = 0.3-2.5, DR = 1.2-2.0, s_D = 2-6\n"
        "Source: Baldauf et al. (2002)\n\n"
        "Returns: Laterally averaged effectiveness eta_avg [-]"
    );

    m.def(
        "film_cooling_multirow_sellers",
        &combaero::cooling::film_cooling_multirow_sellers,
        py::arg("row_positions_xD"),
        py::arg("eval_xD"),
        py::arg("M"),
        py::arg("DR"),
        py::arg("alpha_deg"),
        "Multi-row film cooling using Sellers (1963) superposition.\n\n"
        "Combines effectiveness from multiple upstream rows using the\n"
        "Sellers superposition principle: eta_total = 1 - product(1 - eta_i).\n\n"
        "Parameters:\n"
        "  row_positions_xD : streamwise positions of hole rows [x/D]\n"
        "  eval_xD          : evaluation location [x/D]\n"
        "  M                : blowing ratio [-]\n"
        "  DR               : density ratio [-]\n"
        "  alpha_deg        : injection angle [deg]\n\n"
        "Returns: total adiabatic effectiveness [-]\n\n"
        "References:\n"
        "  - Sellers (1963): Superposition principle for multiple rows\n"
        "  - Baldauf et al. (2002): Single-row effectiveness correlation\n\n"
        "Accuracy: +/-15-20% (flat plate; curvature requires Ito correction)\n\n"
        "Example:\n"
        "  >>> rows = [0, 10, 20]  # Three rows at x/D = 0, 10, 20\n"
        "  >>> eta = cb.film_cooling_multirow_sellers(rows, eval_xD=30, M=0.5, DR=1.8, alpha_deg=30)\n"
        "  >>> print(eta)  # Combined effectiveness at x/D = 30"
    );

    m.def(
        "effusion_effectiveness",
        &combaero::cooling::effusion_effectiveness,
        py::arg("x_D"),
        py::arg("M"),
        py::arg("DR"),
        py::arg("porosity"),
        py::arg("s_D"),
        py::arg("alpha_deg"),
        "Effusion cooling effectiveness for high-density hole arrays.\n\n"
        "Used in combustor liners with continuous coolant injection.\n"
        "Based on Lefebvre (1984) momentum flux ratio approach with\n"
        "crossflow accumulation effects (L'Ecuyer & Matsuura 1985).\n\n"
        "Parameters:\n"
        "  x_D      : streamwise distance from first row / hole diameter [-]\n"
        "  M        : blowing ratio (rho_c*v_c)/(rho_inf*v_inf) [-]\n"
        "  DR       : density ratio rho_c/rho_inf [-]\n"
        "  porosity : hole area / total surface area [-]\n"
        "  s_D      : hole pitch / hole diameter [-]\n"
        "  alpha_deg: injection angle [degrees]\n\n"
        "Returns: adiabatic effectiveness eta [-] (0 = no cooling, 1 = perfect)\n\n"
        "Valid ranges:\n"
        "  M = 1-4, DR = 1.2-2.0, porosity = 0.02-0.10\n"
        "  s_D = 4-8, alpha = 20-45 deg\n\n"
        "Accuracy: +/-20-25% (high-density arrays, combustor liner geometry)\n\n"
        "References:\n"
        "  - Lefebvre (1984): Gas Turbine Combustion\n"
        "  - L'Ecuyer & Matsuura (1985): Crossflow effects\n\n"
        "Example:\n"
        "  >>> eta = cb.effusion_effectiveness(x_D=10, M=2.0, DR=1.8, porosity=0.05, s_D=6.0, alpha_deg=30)\n"
        "  >>> print(eta)  # Effectiveness at x/D = 10"
    );

    m.def(
        "effusion_discharge_coefficient",
        &combaero::cooling::effusion_discharge_coefficient,
        py::arg("Re_d"),
        py::arg("P_ratio"),
        py::arg("alpha_deg"),
        py::arg("L_D") = 4.0,
        "Effusion hole discharge coefficient.\n\n"
        "Accounts for compressibility, inclination, and hole geometry.\n"
        "Uses isentropic flow relations for subsonic flow.\n\n"
        "Parameters:\n"
        "  Re_d      : hole Reynolds number based on diameter [-]\n"
        "  P_ratio   : plenum/mainstream pressure ratio [-]\n"
        "  alpha_deg : hole inclination angle [degrees]\n"
        "  L_D       : hole length / diameter [-] (default: 4.0)\n\n"
        "Returns: discharge coefficient Cd [-]\n\n"
        "Valid ranges:\n"
        "  Re_d > 3000, 1.02 < P_ratio < 1.15\n"
        "  20 < alpha < 45 deg, 2 < L_D < 8\n\n"
        "Accuracy: +/-10-15%\n\n"
        "References:\n"
        "  - Hay & Spencer (1992): Compressible flow through inclined holes\n"
        "  - Gritsch et al. (1998): Shaped hole discharge\n\n"
        "Example:\n"
        "  >>> Cd = cb.effusion_discharge_coefficient(Re_d=10000, P_ratio=1.05, alpha_deg=30, L_D=4.0)\n"
        "  >>> print(Cd)  # Discharge coefficient"
    );

    m.def(
        "pin_fin_nusselt",
        &combaero::cooling::pin_fin_nusselt,
        py::arg("Re_d"),
        py::arg("Pr"),
        py::arg("L_D"),
        py::arg("S_D"),
        py::arg("X_D"),
        py::arg("is_staggered") = true,
        "Pin fin array Nusselt number for internal cooling.\n\n"
        "Staggered or inline arrangements of cylindrical pins.\n"
        "Based on Metzger et al. (1982) correlations.\n\n"
        "Parameters:\n"
        "  Re_d         : Reynolds number based on pin diameter [-]\n"
        "  Pr           : Prandtl number [-]\n"
        "  L_D          : pin length / diameter [-]\n"
        "  S_D          : spanwise spacing / diameter [-]\n"
        "  X_D          : streamwise spacing / diameter [-]\n"
        "  is_staggered : True for staggered, False for inline (default: True)\n\n"
        "Returns: Average Nusselt number Nu_d [-]\n\n"
        "Valid ranges:\n"
        "  Re_d = 3000-90000, Pr = 0.5-1.0, L_D = 0.5-4.0\n"
        "  S_D = 1.5-4.0, X_D = 1.5-4.0\n\n"
        "Accuracy: +/-15-20%\n\n"
        "References:\n"
        "  - Metzger et al. (1982): Effects of Pin Shape and Array Orientation\n\n"
        "Example:\n"
        "  >>> Nu = cb.pin_fin_nusselt(Re_d=20000, Pr=0.7, L_D=2.0, S_D=2.5, X_D=2.5)\n"
        "  >>> print(Nu)  # Nusselt number for staggered array"
    );

    m.def(
        "dimple_nusselt_enhancement",
        &combaero::cooling::dimple_nusselt_enhancement,
        py::arg("Re_Dh"),
        py::arg("d_Dh"),
        py::arg("h_d"),
        py::arg("S_d"),
        "Dimpled surface heat transfer enhancement factor.\n\n"
        "Semi-spherical indentations for internal cooling.\n"
        "Returns enhancement ratio Nu_dimple/Nu_smooth.\n\n"
        "Parameters:\n"
        "  Re_Dh : Reynolds number based on channel hydraulic diameter [-]\n"
        "  d_Dh  : dimple diameter / channel height [-]\n"
        "  h_d   : dimple depth / diameter [-]\n"
        "  S_d   : dimple spacing / diameter [-]\n\n"
        "Returns: Enhancement factor Nu_dimple/Nu_smooth [-]\n\n"
        "Valid ranges:\n"
        "  Re_Dh = 10000-80000, d_Dh = 0.1-0.3\n"
        "  h_d = 0.1-0.3, S_d = 1.5-3.0\n\n"
        "Accuracy: +/-20%\n\n"
        "References:\n"
        "  - Chyu et al. (1997): Heat Transfer of Arrays of Semi-Spherical Indentations\n\n"
        "Key advantage: Lower friction penalty than ribs (1.5-2x vs 6-10x)\n\n"
        "Example:\n"
        "  >>> enhancement = cb.dimple_nusselt_enhancement(Re_Dh=30000, d_Dh=0.2, h_d=0.2, S_d=2.0)\n"
        "  >>> print(enhancement)  # Typically 1.8-2.5x"
    );

    m.def(
        "dimple_friction_multiplier",
        &combaero::cooling::dimple_friction_multiplier,
        py::arg("Re_Dh"),
        py::arg("d_Dh"),
        py::arg("h_d"),
        "Dimpled surface friction multiplier.\n\n"
        "Friction penalty for dimpled surfaces (f_dimple/f_smooth).\n\n"
        "Parameters:\n"
        "  Re_Dh : Reynolds number based on channel hydraulic diameter [-]\n"
        "  d_Dh  : dimple diameter / channel height [-]\n"
        "  h_d   : dimple depth / diameter [-]\n\n"
        "Returns: Friction multiplier f_dimple/f_smooth [-]\n\n"
        "Valid ranges:\n"
        "  Re_Dh = 10000-80000, d_Dh = 0.1-0.3, h_d = 0.1-0.3\n\n"
        "Accuracy: +/-15%\n\n"
        "References:\n"
        "  - Chyu et al. (1997)\n\n"
        "Typical range: 1.3-2.2x (much lower than ribs)\n\n"
        "Example:\n"
        "  >>> f_ratio = cb.dimple_friction_multiplier(Re_Dh=30000, d_Dh=0.2, h_d=0.2)\n"
        "  >>> print(f_ratio)  # Typically 1.5-2.0x"
    );

    m.def(
        "adiabatic_wall_temperature",
        &combaero::cooling::adiabatic_wall_temperature,
        py::arg("T_hot"),
        py::arg("T_coolant"),
        py::arg("eta"),
        "Adiabatic wall temperature from cooling effectiveness.\n\n"
        "Uses: eta = (T_hot - T_aw) / (T_hot - T_coolant).\n\n"
        "Parameters:\n"
        "  T_hot     : mainstream hot-gas temperature [K]\n"
        "  T_coolant : coolant supply temperature [K]\n"
        "  eta       : adiabatic effectiveness [-]\n\n"
        "Returns: adiabatic wall temperature T_aw [K]"
    );

    m.def(
        "cooled_wall_heat_flux",
        &combaero::cooling::cooled_wall_heat_flux,
        py::arg("T_hot"),
        py::arg("T_coolant"),
        py::arg("h_hot"),
        py::arg("h_coolant"),
        py::arg("eta"),
        py::arg("t_wall"),
        py::arg("k_wall"),
        "Film/effusion-cooled wall heat flux with single wall layer.\n\n"
        "Computes q = U * (T_aw - T_coolant), where T_aw is from cooling\n"
        "effectiveness and U includes both convection sides plus wall conduction.\n\n"
        "Parameters:\n"
        "  T_hot      : mainstream hot-gas temperature [K]\n"
        "  T_coolant  : coolant supply temperature [K]\n"
        "  h_hot      : hot-side HTC [W/(m^2*K)]\n"
        "  h_coolant  : coolant-side HTC [W/(m^2*K)]\n"
        "  eta        : adiabatic effectiveness [-]\n"
        "  t_wall     : wall thickness [m]\n"
        "  k_wall     : wall thermal conductivity [W/(m*K)]\n\n"
        "Returns: wall heat flux q [W/m^2]"
    );

    // =========================================================================
    // Units System
    // =========================================================================

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
