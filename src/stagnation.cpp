// Isentropic stagnation/static conversion utilities.
// Uses the variable-cp thermo backend (thermo.h) throughout.

#include "../include/stagnation.h"
#include "../include/thermo.h"
#include "../include/transport.h"
#include <cmath>
#include <stdexcept>

// -------------------------------------------------------------
// Kinetic energy / stagnation enthalpy
// -------------------------------------------------------------

double h0_from_static(double h_static_J_per_kg, double v) {
    return h_static_J_per_kg + 0.5 * v * v;
}

double v_from_h0(double h0, double h_static_J_per_kg) {
    double dh = h0 - h_static_J_per_kg;
    return (dh > 0.0) ? std::sqrt(2.0 * dh) : 0.0;
}

// -------------------------------------------------------------
// Mach number
// -------------------------------------------------------------

double mach_number(double v, double T, const std::vector<double>& X) {
    double a = speed_of_sound(T, X);
    if (a <= 0.0) {
        throw std::invalid_argument("mach_number: speed of sound is non-positive");
    }
    return v / a;
}

// -------------------------------------------------------------
// Stagnation temperature
// -------------------------------------------------------------

// Helper: h_mass [J/kg] at temperature T for composition X
static double h_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return h(T, X) / mw_g * 1000.0;  // J/mol / (g/mol) * 1000 = J/kg
}

// Helper: cp_mass [J/(kg*K)] at temperature T for composition X
static double cp_mass_at(double T, const std::vector<double>& X) {
    double mw_g = mwmix(X);  // g/mol
    return cp(T, X) / mw_g * 1000.0;  // J/(mol*K) / (g/mol) * 1000 = J/(kg*K)
}

double T0_from_static_v(double T, double v, const std::vector<double>& X,
                        double tol, std::size_t max_iter) {
    if (T <= 0.0) {
        throw std::invalid_argument("T0_from_static_v: T must be positive");
    }
    if (v < 0.0) {
        throw std::invalid_argument("T0_from_static_v: v must be non-negative");
    }

    double h_static = h_mass_at(T, X);
    double h0 = h_static + 0.5 * v * v;

    // Newton iteration: find T0 such that h(T0) = h0
    // F(T0) = h(T0) - h0 = 0,  dF/dT0 = cp(T0)
    double T0 = T + 0.5 * v * v / cp_mass_at(T, X);  // ideal-gas initial guess

    constexpr double T_MIN = 200.0;
    constexpr double T_MAX = 6000.0;

    for (std::size_t it = 0; it < max_iter; ++it) {
        double F = h_mass_at(T0, X) - h0;
        double dF = cp_mass_at(T0, X);

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T0;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T0 += dT;
        T0 = std::max(T_MIN, std::min(T0, T_MAX));

        if (std::abs(dT) < tol * (1.0 + std::abs(T0))) {
            return T0;
        }
    }
    return T0;
}

double T0_from_static(double T, double M, const std::vector<double>& X,
                      double tol, std::size_t max_iter) {
    if (M < 0.0) {
        throw std::invalid_argument("T0_from_static: M must be non-negative");
    }
    double a = speed_of_sound(T, X);
    double v = M * a;
    return T0_from_static_v(T, v, X, tol, max_iter);
}

double T_from_stagnation(double T0, double M, const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (T0 <= 0.0) {
        throw std::invalid_argument("T_from_stagnation: T0 must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("T_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return T0;

    double h0 = h_mass_at(T0, X);

    // Newton iteration: find T such that h(T) + M^2*a(T)^2/2 = h0
    // F(T) = h(T) + M^2*a(T)^2/2 - h0
    // dF/dT = cp(T) + M^2 * a * da/dT
    // da/dT = (1/(2a)) * d(gamma*R*T)/dT  -- approximate as a/(2T) for ideal gas
    double T = T0 / (1.0 + 0.2 * M * M);  // calorically perfect gas initial guess

    constexpr double T_MIN = 200.0;
    constexpr double T_MAX = 6000.0;

    for (std::size_t it = 0; it < max_iter; ++it) {
        double a = speed_of_sound(T, X);
        double v = M * a;
        double F = h_mass_at(T, X) + 0.5 * v * v - h0;
        // dF/dT ~ cp + M^2 * a * (a/(2T)) = cp + M^2 * a^2 / (2T)
        double dF = cp_mass_at(T, X) + M * M * a * a / (2.0 * T);

        if (std::abs(dF) < 1e-30) break;

        double dT = -F / dF;
        double max_step = 0.5 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0.0) ? max_step : -max_step;
        }
        T += dT;
        T = std::max(T_MIN, std::min(T, T_MAX));

        if (std::abs(dT) < tol * (1.0 + std::abs(T))) {
            return T;
        }
    }
    return T;
}

// -------------------------------------------------------------
// Stagnation pressure (isentropic)
// -------------------------------------------------------------

double P0_from_static(double P, double T, double M,
                      const std::vector<double>& X,
                      double tol, std::size_t max_iter) {
    if (P <= 0.0) {
        throw std::invalid_argument("P0_from_static: P must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("P0_from_static: M must be non-negative");
    }
    if (M == 0.0) return P;

    // Isentropic: s(T0, P0) = s(T, P)
    // s0 = s(T, P)  (entropy at static conditions)
    // T0 = T0_from_static(T, M, X)
    // Solve for P0: s(T0, P0) = s0
    // s(T, P) = s(T, P_ref) - R*ln(P/P_ref)  =>  s(T0, P0) = s0
    // => s(T0, P_ref) - R*ln(P0/P_ref) = s0
    // => P0 = P_ref * exp((s(T0, P_ref) - s0) / R_specific)

    constexpr double P_REF = 101325.0;
    double s0 = s(T, X, P, P_REF);
    double T0 = T0_from_static(T, M, X, tol, max_iter);

    // specific_gas_constant(X) already returns R_mix [J/(kg*K)]
    double R_specific = specific_gas_constant(X);  // J/(kg*K)

    // s() returns J/(mol*K); convert to mass basis [J/(kg*K)]
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    // s(T0, P0) = s(T0, P_ref) - R_specific * ln(P0/P_ref) = s0
    // => ln(P0/P_ref) = (s(T0, P_ref) - s0) / R_specific
    double s_T0_Pref = s(T0, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass   = s0 / mw_kg;                        // J/(kg*K)
    double ln_ratio = (s_T0_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

double P_from_stagnation(double P0, double T0, double M,
                         const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (P0 <= 0.0) {
        throw std::invalid_argument("P_from_stagnation: P0 must be positive");
    }
    if (M < 0.0) {
        throw std::invalid_argument("P_from_stagnation: M must be non-negative");
    }
    if (M == 0.0) return P0;

    constexpr double P_REF = 101325.0;
    double s0 = s(T0, X, P0, P_REF);
    double T_static = T_from_stagnation(T0, M, X, tol, max_iter);

    double R_specific = specific_gas_constant(X);  // J/(kg*K)
    double mw_kg = mwmix(X) / 1000.0;  // kg/mol

    double s_T_Pref = s(T_static, X, P_REF, P_REF) / mw_kg;  // J/(kg*K)
    double s0_mass  = s0 / mw_kg;                              // J/(kg*K)
    double ln_ratio = (s_T_Pref - s0_mass) / R_specific;
    return P_REF * std::exp(ln_ratio);
}

// -------------------------------------------------------------
// Adiabatic wall temperature
// -------------------------------------------------------------

double recovery_factor(double Pr, bool turbulent) {
    if (Pr <= 0.0) {
        throw std::invalid_argument("recovery_factor: Pr must be positive");
    }
    return turbulent ? std::cbrt(Pr) : std::sqrt(Pr);
}

double T_adiabatic_wall(double T_static, double v,
                        double T, double P, const std::vector<double>& X,
                        bool turbulent) {
    if (T_static <= 0.0) {
        throw std::invalid_argument("T_adiabatic_wall: T_static must be positive");
    }
    if (v < 0.0) {
        throw std::invalid_argument("T_adiabatic_wall: v must be non-negative");
    }
    if (v == 0.0) return T_static;

    double Pr = prandtl(T, P, X);
    double r = recovery_factor(Pr, turbulent);
    double cp_val = cp_mass_at(T, X);

    return T_static + r * 0.5 * v * v / cp_val;
}

double T_adiabatic_wall_mach(double T_static, double M,
                              double T, double P, const std::vector<double>& X,
                              bool turbulent) {
    if (M < 0.0) {
        throw std::invalid_argument("T_adiabatic_wall_mach: M must be non-negative");
    }
    if (M == 0.0) return T_static;

    double a = speed_of_sound(T, X);
    double v = M * a;
    return T_adiabatic_wall(T_static, v, T, P, X, turbulent);
}
