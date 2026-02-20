#ifndef STATE_H
#define STATE_H

#include <string>
#include <vector>

// NOTE: For better performance, consider making T, P, X private with
// getter/setters that invalidate a transport cache. Currently transport
// properties are recomputed on every call.

// Thermodynamic state of a gas mixture at STATIC conditions (T, P).
// All property getters (h, cp, s, rho, a, ...) are evaluated at the
// given static temperature T and static pressure P.
//
// Exception: CompressibleFlowSolution::stagnation holds a State at
// stagnation conditions (T0, P0) — the only deliberate use of State
// to carry total rather than static values.
//
// For the network solver, use MixtureState (NETWORK_ROADMAP.md) which
// carries both static (T, P) and total (T_total, P_total) explicitly.
struct State {
    double T = 298.15;               // Temperature [K]
    double P = 101325.0;             // Pressure [Pa]
    std::vector<double> X;           // Mole fractions [-]

    // Property getters (implemented in state.cpp)
    double mw() const;
    double cp() const;
    double h() const;
    double s() const;
    double cv() const;
    double u() const;
    double rho() const;
    double R() const;
    double gamma() const;
    double a() const;

    // Transport properties (recomputed on each call)
    double mu() const;     // Dynamic viscosity [Pa·s]
    double k() const;      // Thermal conductivity [W/(m·K)]
    double nu() const;     // Kinematic viscosity [m²/s]
    double Pr() const;     // Prandtl number [-]
    double alpha() const;  // Thermal diffusivity [m²/s]

    // Convenience setters (return *this for chaining)
    State& set_T(double T_new) { T = T_new; return *this; }
    State& set_P(double P_new) { P = P_new; return *this; }
    State& set_X(const std::vector<double>& X_new) { X = X_new; return *this; }
};

// A stream is a state with a mass flow rate
struct Stream {
    State state;
    double mdot = 0.0;               // Mass flow rate [kg/s]

    // Convenience accessors
    double T() const { return state.T; }
    double P() const { return state.P; }
    const std::vector<double>& X() const { return state.X; }

    // Property getters (delegate to state)
    double mw() const { return state.mw(); }
    double cp() const { return state.cp(); }
    double h() const { return state.h(); }
    double s() const { return state.s(); }
    double rho() const { return state.rho(); }

    // Setters
    Stream& set_T(double T_new) { state.T = T_new; return *this; }
    Stream& set_P(double P_new) { state.P = P_new; return *this; }
    Stream& set_X(const std::vector<double>& X_new) { state.X = X_new; return *this; }
    Stream& set_mdot(double mdot_new) { mdot = mdot_new; return *this; }
};

// -------------------------------------------------------------
// Stream mixing
// -------------------------------------------------------------

// Mix multiple streams, solving mass and enthalpy balance
// P_out: output pressure (default -1 means use minimum inlet pressure)
Stream mix(const std::vector<Stream>& streams, double P_out = -1.0);

// -------------------------------------------------------------
// Transport State Bundle
// -------------------------------------------------------------

// Bundle of transport properties for a gas mixture
// All properties computed from (T, P, X) in a single call
struct TransportState {
    double T;      // Temperature [K] (input, echoed back)
    double P;      // Pressure [Pa] (input, echoed back)
    double rho;    // Density [kg/m³]
    double mu;     // Dynamic viscosity [Pa·s]
    double k;      // Thermal conductivity [W/(m·K)]
    double nu;     // Kinematic viscosity [m²/s]
    double alpha;  // Thermal diffusivity [m²/s]
    double Pr;     // Prandtl number [-]
    double cp;     // Specific heat at constant pressure [J/(kg·K)]
    double cv;     // Specific heat at constant volume [J/(kg·K)]
    double gamma;  // Heat capacity ratio cp/cv [-]
    double a;      // Speed of sound [m/s]
};

// -------------------------------------------------------------
// Thermo State Bundles
// -------------------------------------------------------------

// Bundle of thermodynamic properties for a gas mixture
// All properties computed from (T, P, X) in a single call
struct ThermoState {
    double T;         // Temperature [K] (input, echoed back)
    double P;         // Pressure [Pa] (input, echoed back)
    double rho;       // Density [kg/m³]
    double cp;        // Specific heat at constant pressure [J/(mol·K)]
    double cv;        // Specific heat at constant volume [J/(mol·K)]
    double h;         // Specific enthalpy [J/mol]
    double s;         // Specific entropy [J/(mol·K)]
    double u;         // Specific internal energy [J/mol]
    double gamma;     // Isentropic expansion coefficient [-]
    double a;         // Speed of sound [m/s]
    double cp_mass;   // Mass-specific cp [J/(kg·K)]
    double cv_mass;   // Mass-specific cv [J/(kg·K)]
    double h_mass;    // Mass-specific enthalpy [J/kg]
    double s_mass;    // Mass-specific entropy [J/(kg·K)]
    double u_mass;    // Mass-specific internal energy [J/kg]
    double mw;        // Molecular weight [g/mol]
};

// AirProperties is an alias for TransportState.
// air_properties(T, P, humidity) returns a TransportState with all fields populated.
using AirProperties = TransportState;

// -------------------------------------------------------------
// Equilibrium Result Bundle
// -------------------------------------------------------------

// Result of any equilibrium solve — carries the final state plus diagnostics.
// Designed to be network-node-friendly: converged flag propagates through a
// network without silently passing bad values.
struct EquilibriumResult {
    State state;       // Equilibrium T, P, X
    double T_in;       // Input temperature [K]
    double delta_T;    // state.T - T_in [K] (negative = endothermic)
    bool converged;    // Solver convergence flag
};

// -------------------------------------------------------------
// Combustion Method
// -------------------------------------------------------------

// Selects the product model used by combustion_state / combustion_state_from_streams.
//   Complete    : complete combustion to CO2 + H2O (fast, default)
//   Equilibrium : complete combustion followed by reforming + WGS equilibrium
//                 (slower, physically correct for rich or high-T conditions)
// Designed for extension: a future WGS-only or full-Gibbs option can be added here.
enum class CombustionMethod {
    Complete,
    Equilibrium,
};

// Bundle of thermodynamic and transport properties for a gas mixture
// Combines ThermoState + TransportState in a single call
struct CompleteState {
    ThermoState thermo;       // All 16 thermodynamic properties
    TransportState transport; // All 9 transport properties
};

// -------------------------------------------------------------
// Combustion State Bundle
// -------------------------------------------------------------

// Bundle of combustion properties with nested CompleteState objects
struct CombustionState {
    double phi;                    // Equivalence ratio [-]
    std::string fuel_name;         // Optional fuel label (e.g., "CH4", "Natural Gas")
    CombustionMethod method;       // Product model used (Complete or Equilibrium)
    CompleteState reactants;       // Reactant state (all thermo + transport properties)
    CompleteState products;        // Product state (all thermo + transport properties)
    double mixture_fraction;       // Bilger mixture fraction [-]
    double fuel_burn_fraction;     // Fraction of fuel burned [0-1]
};

#endif // STATE_H
