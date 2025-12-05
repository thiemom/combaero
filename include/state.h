#ifndef STATE_H
#define STATE_H

#include <vector>

// Thermodynamic state of a gas mixture
struct State {
    double T;                        // Temperature [K]
    double P = 101325.0;             // Pressure [Pa]
    std::vector<double> X;           // Mole fractions [-]
};

// A stream is a state with a mass flow rate
struct Stream {
    State state;
    double mdot = 0.0;               // Mass flow rate [kg/s]
};

#endif // STATE_H
