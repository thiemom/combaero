// Compressible isentropic flow for ideal gas with variable cp(T).
// Solves isentropic relations numerically using Newton iteration.

#include "compressible.h"
#include "thermo.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

using combaero::thermo::R_GAS;

namespace {

// Reference pressure for entropy calculations
constexpr double P_REF = 101325.0;

// Solve for temperature at pressure P such that s(T, P) = s0 (isentropic)
// Uses Newton iteration with ds/dT = cp/T
double solve_T_isentropic(double P, double s0, double T_init,
                          const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    double T = T_init;
    
    // Clamp temperature bounds to valid thermo data range
    constexpr double T_MIN = 200.0;
    constexpr double T_MAX = 6000.0;
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        double s_curr = s(T, X, P, P_REF);
        double F = s_curr - s0;
        
        // ds/dT at constant P = cp/T
        double dF = cp(T, X) / T;
        
        if (std::abs(dF) < 1e-30) break;
        
        double dT = -F / dF;
        
        // Limit step size to avoid overshooting
        double max_step = 0.5 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0) ? max_step : -max_step;
        }
        
        T += dT;
        
        // Clamp temperature to reasonable bounds
        T = std::max(T_MIN, std::min(T, T_MAX));
        
        if (std::abs(dT) < tol * (1.0 + std::abs(T))) {
            return T;
        }
    }
    
    return T;
}

// Compute mass flux G = rho * v for isentropic expansion from (T0, P0) to P
// Returns G in kg/(m²·s)
double compute_mass_flux(double T0, double P0, double P, double s0, double h0,
                         const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (P >= P0) return 0.0;
    
    double T = solve_T_isentropic(P, s0, T0, X, tol, max_iter);
    
    double h_curr = h(T, X);
    double dh = h0 - h_curr;
    
    if (dh <= 0.0) return 0.0;
    
    double v = std::sqrt(2.0 * dh);
    double rho_curr = density(T, P, X);
    
    return rho_curr * v;
}

// Golden section search to find pressure that maximizes mass flux (critical point)
double find_critical_pressure(double T0, double P0, double s0, double h0,
                              const std::vector<double>& X,
                              double tol, std::size_t max_iter) {
    const double phi = 0.5 * (3.0 - std::sqrt(5.0));  // Golden ratio conjugate
    
    // For ideal gas, P*/P0 ~ 0.528 for gamma=1.4
    // Search in a reasonable range around this
    double a = 0.3 * P0;   // Lower bound (below typical critical ratio)
    double b = 0.99 * P0;  // Upper bound (just below stagnation)
    
    double x1 = b - phi * (b - a);
    double x2 = a + phi * (b - a);
    
    double f1 = compute_mass_flux(T0, P0, x1, s0, h0, X, tol, max_iter);
    double f2 = compute_mass_flux(T0, P0, x2, s0, h0, X, tol, max_iter);
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        if (std::abs(b - a) < tol * (1.0 + std::abs(a) + std::abs(b))) break;
        
        if (f1 < f2) {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = compute_mass_flux(T0, P0, x2, s0, h0, X, tol, max_iter);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = compute_mass_flux(T0, P0, x1, s0, h0, X, tol, max_iter);
        }
    }
    
    return (f1 > f2) ? x1 : x2;
}

}  // anonymous namespace

// -------------------------------------------------------------
// Forward problem: nozzle flow
// -------------------------------------------------------------

CompressibleFlowSolution nozzle_flow(
    double T0, double P0, double P_back, double A_eff,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (T0 <= 0.0) {
        throw std::invalid_argument("nozzle_flow: T0 must be positive");
    }
    if (P0 <= 0.0) {
        throw std::invalid_argument("nozzle_flow: P0 must be positive");
    }
    if (P_back <= 0.0) {
        throw std::invalid_argument("nozzle_flow: P_back must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("nozzle_flow: A_eff must be positive");
    }
    
    CompressibleFlowSolution sol;
    
    // Set up stagnation state
    sol.stagnation.T = T0;
    sol.stagnation.P = P0;
    sol.stagnation.X = X;
    
    double h0 = sol.stagnation.h();
    double s0 = sol.stagnation.s();
    
    // Find critical pressure (where mass flux is maximum)
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    
    double P_outlet;
    if (P_back <= P_crit) {
        // Choked flow: outlet is at critical conditions
        P_outlet = P_crit;
        sol.choked = true;
    } else {
        // Subsonic flow: outlet is at back pressure
        P_outlet = P_back;
        sol.choked = false;
    }
    
    // Compute outlet state
    double T_outlet = solve_T_isentropic(P_outlet, s0, T0, X, tol, max_iter);
    sol.outlet.T = T_outlet;
    sol.outlet.P = P_outlet;
    sol.outlet.X = X;
    
    double h_outlet = sol.outlet.h();
    double dh = h0 - h_outlet;
    
    sol.v = (dh > 0.0) ? std::sqrt(2.0 * dh) : 0.0;
    
    double a_outlet = sol.outlet.a();
    sol.M = (a_outlet > 0.0) ? sol.v / a_outlet : 0.0;
    
    double rho_outlet = sol.outlet.rho();
    double G_outlet = rho_outlet * sol.v;
    
    sol.mdot = A_eff * G_outlet;
    
    return sol;
}

// -------------------------------------------------------------
// Inverse problems
// -------------------------------------------------------------

double solve_A_eff_from_mdot(
    double T0, double P0, double P_back, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_A_eff_from_mdot: mdot_target must be positive");
    }
    
    // Compute mass flux at operating conditions
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    // Find critical conditions
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    
    double P_outlet = (P_back <= P_crit) ? P_crit : P_back;
    double G_outlet = compute_mass_flux(T0, P0, P_outlet, s0, h0, X, tol, max_iter);
    
    if (G_outlet <= 0.0) {
        throw std::runtime_error("solve_A_eff_from_mdot: zero mass flux at operating conditions");
    }
    
    return mdot_target / G_outlet;
}

double solve_P_back_from_mdot(
    double T0, double P0, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_P_back_from_mdot: mdot_target must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("solve_P_back_from_mdot: A_eff must be positive");
    }
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    // Find critical conditions
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    double mdot_choked = A_eff * G_crit;
    
    if (mdot_target > mdot_choked * (1.0 + tol)) {
        throw std::runtime_error("solve_P_back_from_mdot: mdot_target exceeds choked mass flow");
    }
    
    // If essentially choked, return critical pressure
    if (mdot_target >= mdot_choked * (1.0 - tol)) {
        return P_crit;
    }
    
    // Bisection search for P_back in [P_crit, P0]
    double G_target = mdot_target / A_eff;
    double P_lo = P_crit;
    double P_hi = P0;
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        double P_mid = 0.5 * (P_lo + P_hi);
        double G_mid = compute_mass_flux(T0, P0, P_mid, s0, h0, X, tol, max_iter);
        
        if (std::abs(G_mid - G_target) < tol * G_target) {
            return P_mid;
        }
        
        // Mass flux decreases as P increases toward P0
        if (G_mid > G_target) {
            P_lo = P_mid;
        } else {
            P_hi = P_mid;
        }
        
        if ((P_hi - P_lo) < tol * P_lo) {
            return P_mid;
        }
    }
    
    return 0.5 * (P_lo + P_hi);
}

double solve_P0_from_mdot(
    double T0, double P_back, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: mdot_target must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: A_eff must be positive");
    }
    if (P_back <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: P_back must be positive");
    }
    
    // Bisection search for P0
    // P0 must be > P_back for any flow
    double P0_lo = P_back * 1.001;
    double P0_hi = P_back * 1000.0;  // Start with large upper bound
    
    // Expand upper bound if needed
    for (int expand = 0; expand < 10; ++expand) {
        auto sol = nozzle_flow(T0, P0_hi, P_back, A_eff, X, tol, max_iter);
        if (sol.mdot >= mdot_target) break;
        P0_hi *= 10.0;
    }
    
    // Bisection
    for (std::size_t it = 0; it < max_iter; ++it) {
        double P0_mid = 0.5 * (P0_lo + P0_hi);
        auto sol = nozzle_flow(T0, P0_mid, P_back, A_eff, X, tol, max_iter);
        
        if (std::abs(sol.mdot - mdot_target) < tol * mdot_target) {
            return P0_mid;
        }
        
        if (sol.mdot < mdot_target) {
            P0_lo = P0_mid;
        } else {
            P0_hi = P0_mid;
        }
        
        if ((P0_hi - P0_lo) < tol * P0_lo) {
            return P0_mid;
        }
    }
    
    return 0.5 * (P0_lo + P0_hi);
}

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

double critical_pressure_ratio(
    double T0, double P0, const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    return P_crit / P0;
}

double mach_from_pressure_ratio(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (P >= P0) return 0.0;
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    double T = solve_T_isentropic(P, s0, T0, X, tol, max_iter);
    
    double h_curr = h(T, X);
    double dh = h0 - h_curr;
    
    if (dh <= 0.0) return 0.0;
    
    double v = std::sqrt(2.0 * dh);
    double a_curr = speed_of_sound(T, X);
    
    return (a_curr > 0.0) ? v / a_curr : 0.0;
}

double mass_flux_isentropic(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    return compute_mass_flux(T0, P0, P, s0, h0, X, tol, max_iter);
}

// -------------------------------------------------------------
// Fanno flow implementation
// -------------------------------------------------------------

namespace {

// Solve for T given P such that h(T) + u²/2 = h0, with u from mass conservation.
// Returns T, and sets u_out and rho_out.
double solve_T_from_energy(
    double P, double h0_mass, double mdot, double A,
    const std::vector<double>& X, double mw_kg,  // mw in kg/mol
    double T_guess,
    double& u_out, double& rho_out,
    double tol = 1e-8, std::size_t max_iter = 50)
{
    double T = T_guess;
    constexpr double T_MIN = 200.0;
    constexpr double T_MAX = 6000.0;
    double mw_g = mw_kg * 1000.0;  // g/mol for enthalpy conversion
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        // Density from ideal gas: rho = P * mw / (R * T), mw in kg/mol
        double rho = P * mw_kg / (R_GAS * T);
        
        // Velocity from mass conservation: u = mdot / (rho * A)
        double u = mdot / (rho * A);
        
        // Specific enthalpy [J/kg]
        // h() returns J/mol, mw_g in g/mol
        // h_mass = h_mol / mw_g * 1000 = J/mol / (g/mol) * 1000 = J/g * 1000 = J/kg
        double h_mol = h(T, X);  // J/mol
        double h_mass = h_mol / mw_g * 1000.0;  // J/kg
        
        // Energy equation residual: h + u²/2 - h0 = 0
        double F = h_mass + 0.5 * u * u - h0_mass;
        
        // Derivative: dF/dT = dh/dT + u * du/dT
        // dh/dT = cp (mass basis)
        // du/dT = d(mdot/(rho*A))/dT = mdot/A * d(1/rho)/dT
        //       = mdot/A * d(R*T/(P*mw_kg))/dT = mdot*R/(A*P*mw_kg) = u/T
        double cp_mol = cp(T, X);  // J/(mol·K)
        double cp_mass = cp_mol / mw_g * 1000.0;  // J/(kg·K)
        double du_dT = u / T;
        double dF = cp_mass + u * du_dT;
        
        if (std::abs(dF) < 1e-30) break;
        
        double dT = -F / dF;
        
        // Limit step size
        double max_step = 0.3 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0) ? max_step : -max_step;
        }
        
        T += dT;
        T = std::clamp(T, T_MIN, T_MAX);
        
        if (std::abs(dT) < tol * T) {
            break;
        }
    }
    
    // Final values
    rho_out = P * mw_kg / (R_GAS * T);
    u_out = mdot / (rho_out * A);
    
    return T;
}

// Compute dp/dx from momentum equation: dp/dx = -f/(2D) * rho * u²
double dpdx_fanno(double rho, double u, double f, double D) {
    return -f / (2.0 * D) * rho * u * u;
}

}  // namespace

FannoSolution fanno_pipe(
    double T_in, double P_in, double u_in,
    double L, double D, double f,
    const std::vector<double>& X,
    std::size_t n_steps,
    bool store_profile)
{
    // Validate inputs
    if (T_in <= 0.0) {
        throw std::invalid_argument("fanno_pipe: T_in must be positive");
    }
    if (P_in <= 0.0) {
        throw std::invalid_argument("fanno_pipe: P_in must be positive");
    }
    if (L <= 0.0) {
        throw std::invalid_argument("fanno_pipe: L must be positive");
    }
    if (D <= 0.0) {
        throw std::invalid_argument("fanno_pipe: D must be positive");
    }
    if (f < 0.0) {
        throw std::invalid_argument("fanno_pipe: f must be non-negative");
    }
    if (n_steps == 0) {
        throw std::invalid_argument("fanno_pipe: n_steps must be positive");
    }
    
    FannoSolution sol;
    sol.L = L;
    sol.D = D;
    sol.f = f;
    
    // Inlet state
    sol.inlet.T = T_in;
    sol.inlet.P = P_in;
    sol.inlet.X = X;
    double mw_g = sol.inlet.mw();  // g/mol
    double mw_kg = mw_g / 1000.0;  // kg/mol
    double A = M_PI * D * D / 4.0;  // m²
    
    // Inlet density and mass flow
    double rho_in = sol.inlet.rho();
    sol.mdot = rho_in * u_in * A;
    
    // Stagnation enthalpy (conserved): h0 = h + u²/2
    double h_in_mol = h(T_in, X);  // J/mol
    double h_in_mass = h_in_mol / mw_g * 1000.0;  // J/kg
    sol.h0 = h_in_mass + 0.5 * u_in * u_in;
    
    // Check inlet Mach number
    double a_in = sol.inlet.a();
    double M_in = u_in / a_in;
    if (M_in >= 1.0) {
        throw std::invalid_argument("fanno_pipe: inlet flow is supersonic (M >= 1), not supported");
    }
    
    // Store inlet profile point
    if (store_profile) {
        FannoStation st;
        st.x = 0.0;
        st.P = P_in;
        st.T = T_in;
        st.rho = rho_in;
        st.u = u_in;
        st.M = M_in;
        st.h = h_in_mass;
        st.s = sol.inlet.s() / mw_g * 1000.0;  // J/(kg·K)
        sol.profile.push_back(st);
    }
    
    // Integration using RK4
    double dx = L / static_cast<double>(n_steps);
    double x = 0.0;
    double P = P_in;
    double T = T_in;
    double rho = rho_in;
    double u = u_in;
    
    for (std::size_t step = 0; step < n_steps; ++step) {
        // RK4 integration of dp/dx
        // k1
        double k1 = dpdx_fanno(rho, u, f, D);
        
        // k2: evaluate at x + dx/2, P + k1*dx/2
        double P2 = P + 0.5 * k1 * dx;
        if (P2 <= 0.0) { sol.choked = true; sol.L_choke = x; break; }
        double T2, u2, rho2;
        T2 = solve_T_from_energy(P2, sol.h0, sol.mdot, A, X, mw_kg, T, u2, rho2);
        double k2 = dpdx_fanno(rho2, u2, f, D);
        
        // k3: evaluate at x + dx/2, P + k2*dx/2
        double P3 = P + 0.5 * k2 * dx;
        if (P3 <= 0.0) { sol.choked = true; sol.L_choke = x; break; }
        double T3, u3, rho3;
        T3 = solve_T_from_energy(P3, sol.h0, sol.mdot, A, X, mw_kg, T, u3, rho3);
        double k3 = dpdx_fanno(rho3, u3, f, D);
        
        // k4: evaluate at x + dx, P + k3*dx
        double P4 = P + k3 * dx;
        if (P4 <= 0.0) { sol.choked = true; sol.L_choke = x; break; }
        double T4, u4, rho4;
        T4 = solve_T_from_energy(P4, sol.h0, sol.mdot, A, X, mw_kg, T, u4, rho4);
        double k4 = dpdx_fanno(rho4, u4, f, D);
        
        // Update P
        double P_new = P + dx * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        
        if (P_new <= 0.0) {
            sol.choked = true;
            sol.L_choke = x;
            break;
        }
        
        // Update state
        x += dx;
        P = P_new;
        T = solve_T_from_energy(P, sol.h0, sol.mdot, A, X, mw_kg, T, u, rho);
        
        // Check for choking (M approaching 1)
        State current;
        current.T = T;
        current.P = P;
        current.X = X;
        double a = current.a();
        double M = u / a;
        
        if (M >= 0.999) {
            sol.choked = true;
            sol.L_choke = x;
            break;
        }
        
        // Store profile point
        if (store_profile) {
            FannoStation st;
            st.x = x;
            st.P = P;
            st.T = T;
            st.rho = rho;
            st.u = u;
            st.M = M;
            st.h = h(T, X) / mw_g * 1000.0;
            st.s = current.s() / mw_g * 1000.0;
            sol.profile.push_back(st);
        }
    }
    
    // Set outlet state
    sol.outlet.T = T;
    sol.outlet.P = P;
    sol.outlet.X = X;
    
    return sol;
}

FannoSolution fanno_pipe(
    const State& inlet, double u_in,
    double L, double D, double f,
    std::size_t n_steps,
    bool store_profile)
{
    return fanno_pipe(inlet.T, inlet.P, u_in, L, D, f, inlet.X, n_steps, store_profile);
}

double fanno_max_length(
    double T_in, double P_in, double u_in,
    double D, double f,
    const std::vector<double>& X,
    double tol, std::size_t max_iter)
{
    // Binary search for length that causes choking
    // Start with an estimate based on ideal gas relations
    
    State inlet;
    inlet.T = T_in;
    inlet.P = P_in;
    inlet.X = X;
    double M_in = u_in / inlet.a();
    
    if (M_in >= 1.0) {
        return 0.0;  // Already choked
    }
    
    // Initial estimate: use ideal gas Fanno relation
    // 4fL*/D = (1-M²)/(γM²) + (γ+1)/(2γ) * ln((γ+1)M² / (2 + (γ-1)M²))
    double gamma = inlet.gamma();
    double M2 = M_in * M_in;
    double term1 = (1.0 - M2) / (gamma * M2);
    double term2 = (gamma + 1.0) / (2.0 * gamma) * 
                   std::log((gamma + 1.0) * M2 / (2.0 + (gamma - 1.0) * M2));
    double L_star_estimate = D * (term1 + term2) / (4.0 * f);
    
    // Binary search
    double L_low = 0.0;
    double L_high = 2.0 * L_star_estimate;
    
    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double L_mid = 0.5 * (L_low + L_high);
        
        try {
            auto sol = fanno_pipe(T_in, P_in, u_in, L_mid, D, f, X, 200, false);
            
            if (sol.choked) {
                L_high = L_mid;
            } else {
                // Check outlet Mach
                double A = M_PI * D * D / 4.0;
                double u_out = sol.mdot / (sol.outlet.rho() * A);
                double M_out = u_out / sol.outlet.a();
                
                if (M_out > 0.99) {
                    L_high = L_mid;
                } else {
                    L_low = L_mid;
                }
            }
        } catch (...) {
            L_high = L_mid;
        }
        
        if ((L_high - L_low) / L_high < tol) {
            break;
        }
    }
    
    return 0.5 * (L_low + L_high);
}
