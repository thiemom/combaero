#include "turbomachinery.h"

#include "composition.h"
#include "thermo.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace combaero {

namespace {

using thermo::R_GAS;

// S(T) = sum_k X_k s_k(T) at the reference pressure [J/(mol K)]: the mixing
// term and the pressure level cancel in every difference taken here.
double s_ref(double T, const std::vector<double>& X) { return s(T, X, 101325.0); }

// T with S(T) = S_target, starting from T_guess.
double T_from_s_ref(double S_target, const std::vector<double>& X, double T_guess) {
    return calc_T_from_s(S_target, 101325.0, X, T_guess, 1.0e-10, 100);
}

// dT/dX_k of the root of S(T) = S(T1) + c at fixed c: the species' own
// entropy change over [T1, T] against the mixture's slope cp/T there.
std::vector<double> dT_dX_isentrope(double T1, double T, double cp_T,
                                    const std::size_t n) {
    std::vector<double> g(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        g[k] = -(s_species(k, T) - s_species(k, T1)) / (cp_T / T);
    }
    return g;
}

}  // namespace

TurbomachinePath turbomachine_path(double T1, const std::vector<double>& X, double PR,
                                   double eta, bool polytropic, bool expansion) {
    if (!(PR >= 1.0)) {
        throw std::invalid_argument(
            "turbomachine_path: pressure ratio " + std::to_string(PR) +
            " < 1 (PR is P2/P1 for a compression, P1/P2 for an expansion)");
    }
    if (!(eta > 0.0 && eta <= 1.0)) {
        throw std::invalid_argument("turbomachine_path: efficiency " + std::to_string(eta) +
                                    " is outside (0, 1]: above 1 the path would destroy "
                                    "entropy");
    }
    if (!(T1 > 0.0)) {
        throw std::invalid_argument("turbomachine_path: T1 must be > 0");
    }
    const std::size_t n = X.size();
    const double xi = std::log(PR);
    const double sgn = expansion ? -1.0 : 1.0;  // the sign of the isentropic S shift
    const double S1 = s_ref(T1, X);
    const double cp1 = cp(T1, X);  // [J/(mol K)]
    const double h1 = h(T1, X);    // [J/mol]

    TurbomachinePath out;
    out.dT2_dY.assign(n, 0.0);
    out.ddh_dY.assign(n, 0.0);

    // The isentrope, reported for both definitions.
    const double T2s = T_from_s_ref(S1 + sgn * R_GAS * xi, X, T1);
    const double cp2s = cp(T2s, X);
    const double h2s = h(T2s, X);
    out.T2s = T2s;

    double T2 = T1;
    double dT2_dT1 = 1.0;
    double dT2_dxi = 0.0;
    double dT2_deta = 0.0;
    std::vector<double> dT2_dX(n, 0.0);

    if (polytropic) {
        // S(T2) - S(T1) = c, c = R xi / eta (compression), -eta R xi (expansion).
        const double c = expansion ? -eta * R_GAS * xi : R_GAS * xi / eta;
        const double dc_dxi = expansion ? -eta * R_GAS : R_GAS / eta;
        const double dc_deta = expansion ? -R_GAS * xi : -R_GAS * xi / (eta * eta);
        T2 = T_from_s_ref(S1 + c, X, T1);
        const double slope2 = cp(T2, X) / T2;
        dT2_dT1 = (cp1 / T1) / slope2;
        dT2_dxi = dc_dxi / slope2;
        dT2_deta = dc_deta / slope2;
        dT2_dX = dT_dX_isentrope(T1, T2, cp(T2, X), n);
    } else {
        // h2 = h1 + (h2s - h1) / eta (compression), (1 - eta) h1 + eta h2s
        // (expansion). The isentrope's own derivatives first.
        const double slope2s = cp2s / T2s;
        const double dT2s_dT1 = (cp1 / T1) / slope2s;
        const double dT2s_dxi = sgn * R_GAS / slope2s;
        const std::vector<double> dT2s_dX = dT_dX_isentrope(T1, T2s, cp2s, n);
        double h2 = 0.0;
        double dh2_dT1 = 0.0;
        double dh2_dxi = 0.0;
        double dh2_deta = 0.0;
        std::vector<double> dh2_dX(n, 0.0);
        if (expansion) {
            h2 = (1.0 - eta) * h1 + eta * h2s;
            dh2_dT1 = (1.0 - eta) * cp1 + eta * cp2s * dT2s_dT1;
            dh2_dxi = eta * cp2s * dT2s_dxi;
            dh2_deta = h2s - h1;
            for (std::size_t k = 0; k < n; ++k) {
                dh2_dX[k] = (1.0 - eta) * h_species(k, T1) +
                            eta * (h_species(k, T2s) + cp2s * dT2s_dX[k]);
            }
        } else {
            h2 = h1 + (h2s - h1) / eta;
            dh2_dT1 = cp1 + (cp2s * dT2s_dT1 - cp1) / eta;
            dh2_dxi = cp2s * dT2s_dxi / eta;
            dh2_deta = -(h2s - h1) / (eta * eta);
            for (std::size_t k = 0; k < n; ++k) {
                dh2_dX[k] = h_species(k, T1) +
                            (h_species(k, T2s) + cp2s * dT2s_dX[k] - h_species(k, T1)) / eta;
            }
        }
        T2 = calc_T_from_h(h2, X, T2s, 1.0e-10, 100);
        const double cp2 = cp(T2, X);
        dT2_dT1 = dh2_dT1 / cp2;
        dT2_dxi = dh2_dxi / cp2;
        dT2_deta = dh2_deta / cp2;
        for (std::size_t k = 0; k < n; ++k) {
            // h(T2, X) = h2: h_k(T2) + cp2 dT2/dX_k = dh2/dX_k.
            dT2_dX[k] = (dh2_dX[k] - h_species(k, T2)) / cp2;
        }
    }

    out.T2 = T2;
    out.dT2_dT1 = dT2_dT1;
    out.dT2_dPR = dT2_dxi / PR;
    out.dT2_deta = dT2_deta;
    out.dT2_dY = mole_gradient_to_mass(X, dT2_dX);

    // Both efficiencies of the path; at PR = 1 the path is a point and both
    // take the value given.
    const double h2 = h(T2, X);
    const double S2 = s_ref(T2, X);
    if (xi == 0.0 || T2 == T1) {
        out.eta_poly = eta;
        out.eta_is = eta;
    } else if (expansion) {
        out.eta_is = (h1 - h2) / (h1 - h2s);
        out.eta_poly = (S1 - S2) / (R_GAS * xi);
    } else {
        out.eta_is = (h2s - h1) / (h2 - h1);
        out.eta_poly = R_GAS * xi / (S2 - S1);
    }
    if (polytropic) {
        out.eta_poly = eta;
    } else {
        out.eta_is = eta;
    }

    // dh per unit mass, and its Jacobian: (h(T2) - h(T1)) / MW.
    const double MW = mwmix(X) / 1000.0;  // [kg/mol]
    const double cp2 = cp(T2, X);
    out.dh = (h2 - h1) / MW;
    out.ddh_dT1 = (cp2 * dT2_dT1 - cp1) / MW;
    out.ddh_dPR = cp2 * out.dT2_dPR / MW;
    out.ddh_deta = cp2 * dT2_deta / MW;
    std::vector<double> ddh_dX(n, 0.0);
    for (std::size_t k = 0; k < n; ++k) {
        const double dN = h_species(k, T2) - h_species(k, T1) + cp2 * dT2_dX[k];
        ddh_dX[k] = dN / MW - (h2 - h1) * (species_molar_mass(k) / 1000.0) / (MW * MW);
    }
    out.ddh_dY = mole_gradient_to_mass(X, ddh_dX);
    return out;
}

}  // namespace combaero
