#include "incompressible.h"

#include <cmath>
#include <stdexcept>

// Mathematical constant
constexpr double PI = 3.14159265358979323846;

// -------------------------------------------------------------
// Bernoulli equation
// -------------------------------------------------------------

double bernoulli_P2(double P1, double v1, double v2, double rho,
                    double dz, double g)
{
    // P1 + ½ρv1² + ρgh1 = P2 + ½ρv2² + ρgh2
    // P2 = P1 + ½ρ(v1² - v2²) - ρg·dz
    return P1 + 0.5 * rho * (v1 * v1 - v2 * v2) - rho * g * dz;
}

double bernoulli_v2(double P1, double P2, double v1, double rho,
                    double dz, double g)
{
    // P1 + ½ρv1² = P2 + ½ρv2² + ρg·dz
    // v2² = v1² + 2(P1 - P2)/ρ - 2g·dz
    double v2_sq = v1 * v1 + 2.0 * (P1 - P2) / rho - 2.0 * g * dz;
    if (v2_sq < 0.0) {
        throw std::invalid_argument(
            "bernoulli_v2: insufficient energy for flow (v2² < 0)");
    }
    return std::sqrt(v2_sq);
}

// -------------------------------------------------------------
// Orifice / restriction flow
// -------------------------------------------------------------

double orifice_mdot(double P1, double P2, double A, double Cd, double rho)
{
    if (P1 < P2) {
        throw std::invalid_argument("orifice_mdot: P1 must be >= P2");
    }
    if (A <= 0.0 || Cd <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("orifice_mdot: A, Cd, rho must be positive");
    }
    // ṁ = Cd · A · √(2 · ρ · ΔP)
    double dP = P1 - P2;
    return Cd * A * std::sqrt(2.0 * rho * dP);
}

double orifice_Q(double P1, double P2, double A, double Cd, double rho)
{
    // Q = ṁ / ρ
    return orifice_mdot(P1, P2, A, Cd, rho) / rho;
}

double orifice_velocity(double P1, double P2, double rho)
{
    if (P1 < P2) {
        throw std::invalid_argument("orifice_velocity: P1 must be >= P2");
    }
    if (rho <= 0.0) {
        throw std::invalid_argument("orifice_velocity: rho must be positive");
    }
    // v = √(2 · ΔP / ρ)
    return std::sqrt(2.0 * (P1 - P2) / rho);
}

double orifice_area(double mdot, double P1, double P2, double Cd, double rho)
{
    if (P1 <= P2) {
        throw std::invalid_argument("orifice_area: P1 must be > P2");
    }
    if (mdot <= 0.0 || Cd <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("orifice_area: mdot, Cd, rho must be positive");
    }
    // A = ṁ / (Cd · √(2 · ρ · ΔP))
    double dP = P1 - P2;
    return mdot / (Cd * std::sqrt(2.0 * rho * dP));
}

double orifice_dP(double mdot, double A, double Cd, double rho)
{
    if (mdot < 0.0 || A <= 0.0 || Cd <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("orifice_dP: invalid parameters");
    }
    // ΔP = (ṁ / (Cd · A))² / (2 · ρ)
    double term = mdot / (Cd * A);
    return term * term / (2.0 * rho);
}

// -------------------------------------------------------------
// Pipe flow (Darcy-Weisbach)
// -------------------------------------------------------------

double pipe_dP(double v, double L, double D, double f, double rho)
{
    if (D <= 0.0 || L < 0.0 || f < 0.0 || rho <= 0.0) {
        throw std::invalid_argument("pipe_dP: invalid parameters");
    }
    // ΔP = f · (L/D) · (ρ · v² / 2)
    return f * (L / D) * 0.5 * rho * v * v;
}

double pipe_dP_mdot(double mdot, double L, double D, double f, double rho)
{
    double v = pipe_velocity(mdot, D, rho);
    return pipe_dP(v, L, D, f, rho);
}

double pipe_velocity(double mdot, double D, double rho)
{
    if (D <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("pipe_velocity: D and rho must be positive");
    }
    // v = ṁ / (ρ · A) where A = π·D²/4
    double A = PI * D * D / 4.0;
    return mdot / (rho * A);
}

double pipe_mdot(double v, double D, double rho)
{
    if (D <= 0.0 || rho <= 0.0) {
        throw std::invalid_argument("pipe_mdot: D and rho must be positive");
    }
    // ṁ = ρ · v · A where A = π·D²/4
    double A = PI * D * D / 4.0;
    return rho * v * A;
}

// -------------------------------------------------------------
// Hydraulic utilities
// -------------------------------------------------------------

double dynamic_pressure(double v, double rho)
{
    // q = ½ · ρ · v²
    return 0.5 * rho * v * v;
}

double velocity_from_q(double q, double rho)
{
    if (q < 0.0 || rho <= 0.0) {
        throw std::invalid_argument("velocity_from_q: q >= 0 and rho > 0 required");
    }
    // v = √(2 · q / ρ)
    return std::sqrt(2.0 * q / rho);
}
