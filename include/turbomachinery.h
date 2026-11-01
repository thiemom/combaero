#pragma once

#include <vector>

// Compression and expansion along a path of given efficiency (#893): the
// outlet total temperature of a compressor or turbine at a total-to-total
// pressure ratio, for combaero's ideal-gas mixtures (cp a function of T).
//
// Polytropic efficiency eta_p: v dP = eta_p dh (compressor), dh = eta_p v dP
// (turbine), along the path. For an ideal gas, cp dT / T = (R / eta_p) dP / P
// (compressor), so with S(T) = sum_k X_k s_k(T) at the reference pressure
//   S(T2) - S(T1) = (R / eta_p) ln PR       (compression, PR = P2 / P1 >= 1)
//   S(T2) - S(T1) = -eta_p R ln PR          (expansion,   PR = P1 / P2 >= 1)
// which is exact for a variable cp: no path integration is needed.
//
// Isentropic efficiency eta_s: T2s from S(T2s) - S(T1) = +-R ln PR, then
//   h2 = h1 + (h2s - h1) / eta_s            (compression)
//   h2 = h1 - eta_s (h1 - h2s)              (expansion).
//
// The composition is frozen along the path (no reaction). At a fixed PR the
// outlet temperature does not depend on the absolute pressure.

namespace combaero {

struct TurbomachinePath {
    double T2 = 0.0;           // outlet total temperature [K]
    double dh = 0.0;           // h_mass(T2) - h_mass(T1) [J/kg]; > 0 into the flow
    double eta_poly = 0.0;     // polytropic efficiency of the path [-]
    double eta_is = 0.0;       // isentropic efficiency of the path [-]
    double T2s = 0.0;          // isentropic outlet temperature [K]
    // Jacobian of T2
    double dT2_dT1 = 0.0;
    double dT2_dPR = 0.0;
    double dT2_deta = 0.0;     // in the efficiency given
    std::vector<double> dT2_dY;  // in the mass fractions
    // Jacobian of dh
    double ddh_dT1 = 0.0;
    double ddh_dPR = 0.0;
    double ddh_deta = 0.0;
    std::vector<double> ddh_dY;
};

// Compression (expansion = false, PR = P2/P1) or expansion (expansion = true,
// PR = P1/P2), PR >= 1. Exactly one of the efficiencies applies:
// polytropic = true takes eta as eta_poly, false as eta_is; 0 < eta <= 1.
// X are mole fractions. Throws std::invalid_argument on PR < 1 or eta out of
// (0, 1].
TurbomachinePath turbomachine_path(double T1, const std::vector<double>& X, double PR,
                                   double eta, bool polytropic, bool expansion);

}  // namespace combaero
