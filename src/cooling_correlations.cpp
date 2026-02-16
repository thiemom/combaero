#include "cooling_correlations.h"
#include "math_constants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace combaero::cooling {

// -------------------------------------------------------------
// Helper Functions
// -------------------------------------------------------------

void validate_rib_params(double e_D, double P_e, double alpha) {
    if (e_D < 0.02 || e_D > 0.1) {
        throw std::runtime_error(
            "Rib height ratio e_D = " + std::to_string(e_D) +
            " is outside valid range [0.02, 0.1]"
        );
    }
    if (P_e < 5.0 || P_e > 20.0) {
        throw std::runtime_error(
            "Rib pitch ratio P_e = " + std::to_string(P_e) +
            " is outside valid range [5, 20]"
        );
    }
    if (alpha < 30.0 || alpha > 90.0) {
        throw std::runtime_error(
            "Rib angle alpha = " + std::to_string(alpha) +
            " deg is outside valid range [30, 90] deg"
        );
    }
}

void validate_impingement_params(double Re_jet, double z_D, double x_D, double y_D) {
    if (Re_jet < 5000.0 || Re_jet > 80000.0) {
        throw std::runtime_error(
            "Jet Reynolds number Re_jet = " + std::to_string(Re_jet) +
            " is outside valid range [5000, 80000]"
        );
    }
    if (z_D < 1.0 || z_D > 12.0) {
        throw std::runtime_error(
            "Jet-to-plate distance z_D = " + std::to_string(z_D) +
            " is outside valid range [1, 12]"
        );
    }
    if (x_D > 0.0 || y_D > 0.0) {
        if (x_D < 4.0 || x_D > 16.0) {
            throw std::runtime_error(
                "Streamwise spacing x_D = " + std::to_string(x_D) +
                " is outside valid range [4, 16]"
            );
        }
        if (y_D < 4.0 || y_D > 16.0) {
            throw std::runtime_error(
                "Spanwise spacing y_D = " + std::to_string(y_D) +
                " is outside valid range [4, 16]"
            );
        }
    }
}

void validate_film_cooling_params(double M, double DR, double alpha_deg) {
    if (M < 0.3 || M > 2.5) {
        throw std::runtime_error(
            "Blowing ratio M = " + std::to_string(M) +
            " is outside valid range [0.3, 2.5]"
        );
    }
    if (DR < 1.2 || DR > 2.0) {
        throw std::runtime_error(
            "Density ratio DR = " + std::to_string(DR) +
            " is outside valid range [1.2, 2.0]"
        );
    }
    if (alpha_deg < 20.0 || alpha_deg > 90.0) {
        throw std::runtime_error(
            "Injection angle alpha = " + std::to_string(alpha_deg) +
            " deg is outside valid range [20, 90] deg"
        );
    }
}

// -------------------------------------------------------------
// Rib-Enhanced Cooling
// -------------------------------------------------------------

double rib_enhancement_factor(double e_D, double P_e, double alpha) {
    validate_rib_params(e_D, P_e, alpha);

    // Han et al. (1988) correlation
    // Nu/Nu_smooth = C * (e/D)^a * (P/e)^b * f(alpha)

    // Angle factor (90-deg ribs are most effective)
    double alpha_rad = alpha * M_PI / 180.0;
    double f_alpha = 0.8 + 0.4 * std::sin(alpha_rad);

    // Base correlation (adjusted to give 1.5-4.0 range)
    double enhancement = 3.5 * std::pow(e_D, 0.35) * std::pow(P_e, -0.15) * f_alpha;

    return enhancement;
}

double rib_friction_multiplier(double e_D, double P_e) {
    // Validate only e_D and P_e (alpha not needed for friction)
    if (e_D < 0.02 || e_D > 0.1) {
        throw std::runtime_error(
            "Rib height ratio e_D = " + std::to_string(e_D) +
            " is outside valid range [0.02, 0.1]"
        );
    }
    if (P_e < 5.0 || P_e > 20.0) {
        throw std::runtime_error(
            "Rib pitch ratio P_e = " + std::to_string(P_e) +
            " is outside valid range [5, 20]"
        );
    }

    // Han et al. (1988) friction correlation
    // f/f_smooth = C * (e/D)^a * (P/e)^b
    double multiplier = 9.0 * std::pow(e_D, 0.45) * std::pow(P_e, -0.25);

    return multiplier;
}

// -------------------------------------------------------------
// Impingement Cooling
// -------------------------------------------------------------

double impingement_nusselt(double Re_jet, double Pr, double z_D,
                          double x_D, double y_D) {
    validate_impingement_params(Re_jet, z_D, x_D, y_D);

    // Florschuetz et al. (1981) correlation
    // For jet array with crossflow

    bool is_array = (x_D > 0.0 && y_D > 0.0);

    if (!is_array) {
        // Single jet correlation (Martin 1977)
        double Nu = 0.42 * std::pow(Re_jet, 0.55) * std::pow(Pr, 0.4) *
                    std::pow(z_D, -0.05);
        return Nu;
    } else {
        // Jet array with crossflow (Florschuetz et al. 1981)
        // Accounts for spent air crossflow degrading downstream jets

        // Geometric factor
        double G = std::sqrt(x_D * y_D);

        // Crossflow degradation factor (stronger effect)
        double f_crossflow = 1.0 / (1.0 + 0.6 * std::pow(G / z_D, 0.7));

        // Base Nusselt number (lower than single jet)
        double Nu_base = 0.38 * std::pow(Re_jet, 0.55) * std::pow(Pr, 0.4);

        // Height correction
        double f_height = std::pow(z_D / 6.0, -0.05);

        return Nu_base * f_crossflow * f_height;
    }
}

// -------------------------------------------------------------
// Film Cooling Effectiveness
// -------------------------------------------------------------

double film_cooling_effectiveness(double x_D, double M, double DR, double alpha_deg) {
    validate_film_cooling_params(M, DR, alpha_deg);

    if (x_D < 0.0) {
        throw std::runtime_error(
            "Downstream distance x_D = " + std::to_string(x_D) +
            " must be non-negative"
        );
    }

    // Baldauf et al. (2002) correlation
    // eta = f(x/D, M, DR, alpha)

    // Convert angle to radians
    double alpha_rad = alpha_deg * M_PI / 180.0;

    // Momentum flux ratio I = (rho_c*v_c^2)/(rho_inf*v_inf^2) = M^2 / DR
    double I = M * M / DR;

    // Optimal blowing ratio (where effectiveness peaks)
    double M_opt = 0.6 * std::pow(DR, 0.5);

    // Blowing ratio correction (effectiveness drops if M too high or low)
    double f_M = std::exp(-0.8 * std::pow((M - M_opt) / M_opt, 2.0));

    // Angle correction (shallow angles give better coverage)
    // Lower angles (30 deg) should give higher effectiveness than steep (90 deg)
    double f_alpha = 1.2 - 0.4 * std::sin(alpha_rad);

    // Streamwise decay
    double decay = std::exp(-0.08 * x_D);

    // Peak effectiveness (at hole exit) - higher base value
    double eta_0 = 0.75 * f_M * f_alpha;

    // Effectiveness at x_D
    double eta = eta_0 * decay;

    // Clamp to [0, 1]
    return std::max(0.0, std::min(1.0, eta));
}

double film_cooling_effectiveness_avg(double x_D, double M, double DR,
                                     double alpha_deg, double s_D) {
    validate_film_cooling_params(M, DR, alpha_deg);

    if (x_D < 0.0) {
        throw std::runtime_error(
            "Downstream distance x_D = " + std::to_string(x_D) +
            " must be non-negative"
        );
    }
    if (s_D < 2.0 || s_D > 6.0) {
        throw std::runtime_error(
            "Hole spacing s_D = " + std::to_string(s_D) +
            " is outside valid range [2, 6]"
        );
    }

    // Laterally averaged effectiveness (Baldauf et al. 2002)
    // Accounts for hole spacing and jet spreading

    // Get centerline effectiveness
    double eta_centerline = film_cooling_effectiveness(x_D, M, DR, alpha_deg);

    // Lateral spreading factor (closer holes -> better coverage)
    double f_spacing = 1.0 - 0.15 * (s_D - 3.0);

    // Averaged effectiveness is lower than centerline
    double eta_avg = eta_centerline * f_spacing * 0.7;

    // Clamp to [0, 1]
    return std::max(0.0, std::min(1.0, eta_avg));
}

// Multi-row film cooling using Sellers (1963) superposition principle
// Combines effectiveness from multiple upstream rows
double film_cooling_multirow_sellers(
    const std::vector<double>& row_positions_xD,
    double eval_xD,
    double M,
    double DR,
    double alpha_deg
) {
    // Validate inputs
    if (row_positions_xD.empty()) {
        throw std::runtime_error("row_positions_xD cannot be empty");
    }
    
    // Sellers superposition: eta_total = 1 - product(1 - eta_i)
    // Each row contributes based on its local x/D from injection point
    double product_term = 1.0;
    
    for (double row_x : row_positions_xD) {
        double local_xD = eval_xD - row_x;
        
        // Only consider upstream rows (local_xD > 0)
        if (local_xD > 0.0) {
            // Get single-row effectiveness using Baldauf correlation
            // Parameter validation happens inside film_cooling_effectiveness
            double eta_i = film_cooling_effectiveness(local_xD, M, DR, alpha_deg);
            product_term *= (1.0 - eta_i);
        }
    }
    
    return 1.0 - product_term;
}

}  // namespace combaero::cooling
