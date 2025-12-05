#include "../include/state.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include <algorithm>
#include <stdexcept>
#include <limits>

// -------------------------------------------------------------
// State property getters
// -------------------------------------------------------------

double State::mw() const { return mwmix(X); }
double State::cp() const { return ::cp(T, X); }
double State::h() const { return ::h(T, X); }
double State::s() const { return ::s(T, X, P); }
double State::cv() const { return ::cv(T, X); }
double State::u() const { return ::u(T, X); }
double State::rho() const { return ::density(T, P, X); }
double State::R() const { return ::specific_gas_constant(X); }
double State::gamma() const { return ::isentropic_expansion_coefficient(T, X); }
double State::a() const { return ::speed_of_sound(T, X); }

// -------------------------------------------------------------
// Transport properties
// -------------------------------------------------------------

double State::mu() const { return ::viscosity(T, P, X); }
double State::k() const { return ::thermal_conductivity(T, P, X); }
double State::nu() const { return ::kinematic_viscosity(T, P, X); }
double State::Pr() const { return ::prandtl(T, P, X); }
double State::alpha() const { return ::thermal_diffusivity(T, P, X); }

// -------------------------------------------------------------
// Stream mixing
// -------------------------------------------------------------

Stream mix(const std::vector<Stream>& streams, double P_out)
{
    if (streams.empty()) {
        throw std::runtime_error("mix: at least one stream required");
    }

    if (streams.size() == 1) {
        Stream out = streams[0];
        if (P_out > 0.0) {
            out.state.P = P_out;
        }
        return out;
    }

    const std::size_t n_species = streams[0].X().size();

    // Validate all streams have same species count
    for (const auto& s : streams) {
        if (s.X().size() != n_species) {
            throw std::runtime_error("mix: all streams must have same number of species");
        }
    }

    // Determine output pressure (minimum if not specified)
    if (P_out < 0.0) {
        P_out = std::numeric_limits<double>::max();
        for (const auto& s : streams) {
            P_out = std::min(P_out, s.P());
        }
    }

    // Mass balance: total mass flow
    double mdot_total = 0.0;
    for (const auto& s : streams) {
        mdot_total += s.mdot;
    }

    if (mdot_total <= 0.0) {
        throw std::runtime_error("mix: total mass flow must be positive");
    }

    // Species balance: accumulate molar flows
    // n_dot_k = sum_i (mdot_i / MW_i) * X_i[k]
    std::vector<double> n_dot(n_species, 0.0);
    double n_dot_total = 0.0;

    for (const auto& s : streams) {
        double MW_i = s.mw();  // g/mol
        double n_dot_i = s.mdot / (MW_i * 1e-3);  // mol/s (MW in kg/mol)
        for (std::size_t k = 0; k < n_species; ++k) {
            n_dot[k] += n_dot_i * s.X()[k];
        }
        n_dot_total += n_dot_i;
    }

    // Convert to mole fractions
    std::vector<double> X_out(n_species);
    for (std::size_t k = 0; k < n_species; ++k) {
        X_out[k] = n_dot[k] / n_dot_total;
    }

    // Enthalpy balance: H_out * mdot_out = sum_i (H_i * mdot_i)
    // where H is specific enthalpy [J/kg]
    double H_dot_total = 0.0;  // [W] = [J/s]
    for (const auto& s : streams) {
        double h_molar = s.h();  // J/mol
        double MW_i = s.mw();    // g/mol
        double h_mass = h_molar / (MW_i * 1e-3);  // J/kg
        H_dot_total += h_mass * s.mdot;
    }

    double h_out_mass = H_dot_total / mdot_total;  // J/kg

    // Convert to molar enthalpy for T solver
    double MW_out = mwmix(X_out);  // g/mol
    double h_out_molar = h_out_mass * (MW_out * 1e-3);  // J/mol

    // Solve for T_out: h(T_out, X_out) = h_out_molar
    // Use mass-weighted average T as initial guess
    double T_guess = 0.0;
    for (const auto& s : streams) {
        T_guess += s.T() * s.mdot;
    }
    T_guess /= mdot_total;

    double T_out = calc_T_from_h(h_out_molar, X_out, T_guess);

    // Build output stream
    Stream out;
    out.state.T = T_out;
    out.state.P = P_out;
    out.state.X = X_out;
    out.mdot = mdot_total;

    return out;
}
