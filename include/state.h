#ifndef STATE_H
#define STATE_H

#include <vector>

// Thermodynamic state of a gas mixture
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

#endif // STATE_H
