#include "incompressible.h"
#include "friction.h"
#include "transport.h"
#include "thermo.h"

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

// -------------------------------------------------------------
// Pressure loss coefficient (zeta / K)
// -------------------------------------------------------------

double pressure_loss(double v, double rho, double zeta)
{
    if (rho <= 0.0 || zeta < 0.0) {
        throw std::invalid_argument("pressure_loss: rho > 0 and zeta >= 0 required");
    }
    // ΔP = ζ · ½ρv²
    return zeta * 0.5 * rho * v * v;
}

double velocity_from_pressure_loss(double dP, double rho, double zeta)
{
    if (dP < 0.0 || rho <= 0.0 || zeta <= 0.0) {
        throw std::invalid_argument(
            "velocity_from_pressure_loss: dP >= 0, rho > 0, zeta > 0 required");
    }
    // v = √(2·ΔP / (ζ·ρ))
    return std::sqrt(2.0 * dP / (zeta * rho));
}

double zeta_from_Cd(double Cd)
{
    if (Cd <= 0.0) {
        throw std::invalid_argument("zeta_from_Cd: Cd must be positive");
    }
    // ζ = 1 / Cd²
    return 1.0 / (Cd * Cd);
}

double Cd_from_zeta(double zeta)
{
    if (zeta <= 0.0) {
        throw std::invalid_argument("Cd_from_zeta: zeta must be positive");
    }
    // Cd = 1 / √ζ
    return 1.0 / std::sqrt(zeta);
}

// -------------------------------------------------------------
// High-level thermo-aware API
// -------------------------------------------------------------

// Helper: dispatch friction factor from correlation name
static double friction_from_correlation(double Re, double e_D,
                                        const std::string& correlation)
{
    if (correlation == "haaland")  return friction_haaland(Re, e_D);
    if (correlation == "serghides") return friction_serghides(Re, e_D);
    if (correlation == "colebrook") return friction_colebrook(Re, e_D);
    if (correlation == "petukhov")  return friction_petukhov(Re);
    throw std::invalid_argument(
        "incompressible: unknown friction correlation '" + correlation + "'. "
        "Valid options: 'haaland', 'serghides', 'colebrook', 'petukhov'");
}

IncompressibleFlowSolution orifice_flow_thermo(
    double T, double P, const std::vector<double>& X,
    double P_back, double A, double Cd)
{
    if (T <= 0.0)
        throw std::invalid_argument("orifice_flow: T must be positive");
    if (P <= 0.0)
        throw std::invalid_argument("orifice_flow: P must be positive");
    if (P_back < 0.0)
        throw std::invalid_argument("orifice_flow: P_back must be non-negative");
    if (P < P_back)
        throw std::invalid_argument("orifice_flow: P must be >= P_back");
    if (A <= 0.0)
        throw std::invalid_argument("orifice_flow: A must be positive");
    if (Cd <= 0.0 || Cd > 1.0)
        throw std::invalid_argument("orifice_flow: Cd must be in (0, 1]");

    const double rho = density(T, P, X);
    const double dP  = P - P_back;
    const double v   = (dP > 0.0) ? std::sqrt(2.0 * dP / rho) : 0.0;
    const double mdot = Cd * A * rho * v;

    IncompressibleFlowSolution sol;
    sol.mdot = mdot;
    sol.v    = Cd * v;   // effective throat velocity (Cd * ideal velocity)
    sol.dP   = dP;
    sol.Re   = 0.0;      // not meaningful for an orifice (no length scale)
    sol.rho  = rho;
    sol.f    = Cd;       // store Cd in f field for orifice results
    return sol;
}

IncompressibleFlowSolution pipe_flow(
    double T, double P, const std::vector<double>& X,
    double u, double L, double D, double f)
{
    if (T <= 0.0)
        throw std::invalid_argument("pipe_flow: T must be positive");
    if (P <= 0.0)
        throw std::invalid_argument("pipe_flow: P must be positive");
    if (u < 0.0)
        throw std::invalid_argument("pipe_flow: u must be non-negative");
    if (L < 0.0)
        throw std::invalid_argument("pipe_flow: L must be non-negative");
    if (D <= 0.0)
        throw std::invalid_argument("pipe_flow: D must be positive");
    if (f < 0.0)
        throw std::invalid_argument("pipe_flow: f must be non-negative");

    const double rho  = density(T, P, X);
    const double mu   = viscosity(T, P, X);
    const double Re   = (mu > 0.0) ? rho * u * D / mu : 0.0;
    const double dP   = pipe_dP(u, L, D, f, rho);
    const double A    = PI * D * D / 4.0;
    const double mdot = rho * u * A;

    IncompressibleFlowSolution sol;
    sol.mdot = mdot;
    sol.v    = u;
    sol.dP   = dP;
    sol.Re   = Re;
    sol.rho  = rho;
    sol.f    = f;
    return sol;
}

IncompressibleFlowSolution pipe_flow_rough(
    double T, double P, const std::vector<double>& X,
    double u, double L, double D,
    double roughness,
    const std::string& correlation)
{
    if (T <= 0.0)
        throw std::invalid_argument("pipe_flow_rough: T must be positive");
    if (P <= 0.0)
        throw std::invalid_argument("pipe_flow_rough: P must be positive");
    if (u < 0.0)
        throw std::invalid_argument("pipe_flow_rough: u must be non-negative");
    if (L < 0.0)
        throw std::invalid_argument("pipe_flow_rough: L must be non-negative");
    if (D <= 0.0)
        throw std::invalid_argument("pipe_flow_rough: D must be positive");
    if (roughness < 0.0)
        throw std::invalid_argument("pipe_flow_rough: roughness must be non-negative");

    const double rho  = density(T, P, X);
    const double mu   = viscosity(T, P, X);
    const double Re   = (mu > 0.0) ? rho * u * D / mu : 0.0;
    const double e_D  = roughness / D;
    const double f    = friction_from_correlation(Re, e_D, correlation);
    const double dP   = pipe_dP(u, L, D, f, rho);
    const double A    = PI * D * D / 4.0;
    const double mdot = rho * u * A;

    IncompressibleFlowSolution sol;
    sol.mdot = mdot;
    sol.v    = u;
    sol.dP   = dP;
    sol.Re   = Re;
    sol.rho  = rho;
    sol.f    = f;
    return sol;
}

std::tuple<double, double, double> pressure_drop_pipe(
    double T, double P, const std::vector<double>& X,
    double v, double D, double L,
    double roughness,
    const std::string& correlation)
{
    if (T <= 0.0)
        throw std::invalid_argument("pressure_drop_pipe: T must be positive");
    if (P <= 0.0)
        throw std::invalid_argument("pressure_drop_pipe: P must be positive");
    if (D <= 0.0)
        throw std::invalid_argument("pressure_drop_pipe: D must be positive");
    if (L < 0.0)
        throw std::invalid_argument("pressure_drop_pipe: L must be non-negative");
    if (v < 0.0)
        throw std::invalid_argument("pressure_drop_pipe: v must be non-negative");
    if (roughness < 0.0)
        throw std::invalid_argument("pressure_drop_pipe: roughness must be non-negative");

    const double rho = density(T, P, X);
    const double Re  = reynolds(T, P, X, v, D);
    const double e_D = roughness / D;
    const double f   = friction_from_correlation(Re, e_D, correlation);
    const double dP  = pipe_dP(v, L, D, f, rho);

    return std::make_tuple(dP, Re, f);
}
