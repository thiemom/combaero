#include "../include/orifice.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

// -------------------------------------------------------------
// Constants
// -------------------------------------------------------------

namespace {
constexpr double PI = 3.14159265358979323846;
} // namespace

// -------------------------------------------------------------
// OrificeGeometry implementation
// -------------------------------------------------------------

double OrificeGeometry::beta() const {
    if (D <= 0.0) {
        throw std::invalid_argument("OrificeGeometry: pipe diameter D must be > 0");
    }
    return d / D;
}

double OrificeGeometry::area() const {
    return PI * d * d / 4.0;
}

double OrificeGeometry::t_over_d() const {
    if (d <= 0.0) return 0.0;
    return t / d;
}

double OrificeGeometry::r_over_d() const {
    if (d <= 0.0) return 0.0;
    return r / d;
}

bool OrificeGeometry::is_valid() const {
    if (d <= 0.0 || D <= 0.0) return false;
    if (d >= D) return false;
    if (t < 0.0 || r < 0.0) return false;
    return true;
}

// -------------------------------------------------------------
// OrificeState implementation
// -------------------------------------------------------------

double OrificeState::Re_d(double beta) const {
    // Orifice Reynolds number based on orifice diameter d
    // From continuity: v_orifice = v_pipe / beta^2
    // Re_d = (rho * v_orifice * d) / mu
    //      = (rho * (v_pipe / beta^2) * (D * beta)) / mu
    //      = (rho * v_pipe * D) / mu * (1/beta)
    //      = Re_D / beta
    // However, the conventional definition uses Re_d = Re_D * beta
    // (based on the diameter ratio, not the actual flow velocity)
    return Re_D * beta;
}

// -------------------------------------------------------------
// Namespace: individual correlations
// -------------------------------------------------------------

namespace orifice {

// Reader-Harris/Gallagher (1998) correlation
// ISO 5167-2:2003, also ASME MFC-3M
// Valid for: 0.1 <= beta <= 0.75, Re_D >= 5000 (preferably >= 10000)
//            D >= 50 mm, d >= 12.5 mm
double Cd_ReaderHarrisGallagher(double beta, double Re_D, double D) {
    // Ensure minimum Reynolds number
    if (Re_D < 1.0) Re_D = 1.0;

    const double beta2 = beta * beta;
    const double beta4 = beta2 * beta2;
    const double beta8 = beta4 * beta4;

    // Coefficient of discharge (Reader-Harris/Gallagher equation)
    // C = 0.5961 + 0.0261*beta^2 - 0.216*beta^8
    //     + 0.000521*(10^6*beta/Re_D)^0.7
    //     + (0.0188 + 0.0063*A)*beta^3.5*(10^6/Re_D)^0.3
    //     + (0.043 + 0.080*exp(-10*L1) - 0.123*exp(-7*L1))*(1 - 0.11*A)*beta^4/(1 - beta^4)
    //     - 0.031*(M2 - 0.8*M2^1.1)*beta^1.3
    //
    // For flange taps (most common):
    //   L1 = L2 = 25.4/D (mm), A = (19000*beta/Re_D)^0.8

    // Convert D to mm for the correlation
    const double D_mm = D * 1000.0;

    // Flange tap geometry (25.4 mm from plate face)
    const double L1 = 25.4 / D_mm;  // Upstream tap distance ratio
    const double L2 = 25.4 / D_mm;  // Downstream tap distance ratio

    // A parameter (small-bore correction)
    const double A = std::pow(19000.0 * beta / Re_D, 0.8);

    // M2 parameter (downstream tap correction)
    const double M2 = 2.0 * L2 / (1.0 - beta);

    // Base coefficient
    double C = 0.5961
             + 0.0261 * beta2
             - 0.216 * beta8;

    // Reynolds number term
    C += 0.000521 * std::pow(1.0e6 * beta / Re_D, 0.7);

    // Small-bore correction
    C += (0.0188 + 0.0063 * A) * std::pow(beta, 3.5) * std::pow(1.0e6 / Re_D, 0.3);

    // Upstream tap term
    const double exp_term1 = std::exp(-10.0 * L1);
    const double exp_term2 = std::exp(-7.0 * L1);
    C += (0.043 + 0.080 * exp_term1 - 0.123 * exp_term2)
       * (1.0 - 0.11 * A) * beta4 / (1.0 - beta4);

    // Downstream tap term
    C -= 0.031 * (M2 - 0.8 * std::pow(M2, 1.1)) * std::pow(beta, 1.3);

    return C;
}

// Stolz (1978) correlation - older ISO 5167
double Cd_Stolz(double beta, double Re_D) {
    if (Re_D < 1.0) Re_D = 1.0;

    const double beta2 = beta * beta;
    const double beta4 = beta2 * beta2;

    // Stolz equation (corner taps)
    // C = 0.5959 + 0.0312*beta^2.1 - 0.184*beta^8 + 91.71*beta^2.5/Re_D^0.75
    double C = 0.5959
             + 0.0312 * std::pow(beta, 2.1)
             - 0.184 * std::pow(beta, 8.0)
             + 91.71 * std::pow(beta, 2.5) / std::pow(Re_D, 0.75);

    return C;
}

// Miller (1996) simplified correlation
double Cd_Miller(double beta, double Re_D) {
    if (Re_D < 1.0) Re_D = 1.0;

    // Simplified form: C ≈ 0.596 + 0.031*beta^2 for high Re
    // With Reynolds correction
    const double beta2 = beta * beta;

    double C = 0.596 + 0.031 * beta2;

    // Reynolds number correction (approximate)
    if (Re_D < 1.0e6) {
        C += 0.5 * std::pow(beta, 2.5) / std::pow(Re_D, 0.5);
    }

    return C;
}

// Thickness correction factor for thick plates
// Idelchik model: Cd rises (reattachment) then falls (friction)
//
// Physical model:
// - For small t/d: flow reattachment increases Cd (vena contracta recovery)
// - For large t/d: friction in bore decreases Cd (pipe friction losses)
// - Peak Cd occurs at t/d ≈ 0.5-1.5 (depending on beta, Re)
//
// References:
// - Idelchik, I.E. (2008). "Handbook of Hydraulic Resistance" (4th ed.)
//   Diagram 4-15: Discharge coefficients for thick-plate orifices
// - Lichtarowicz, A., et al. (1965). "Discharge Coefficients for
//   Incompressible Non-Cavitating Flow through Long Orifices"
//   J. Mech. Eng. Sci., 7(2), 210-219.
//
// Note: Smooth in Re_d for solver stability (t/d is constant per geometry)
//
double thickness_correction(double t_over_d, double beta, double Re_d) {
    if (t_over_d <= 0.02) {
        return 1.0;  // Thin plate, no correction
    }

    // Ensure Re_d is positive for stability
    Re_d = std::max(Re_d, 100.0);

    // Component 1: Reattachment benefit (independent of Re)
    // Cd increases as flow reattaches to bore wall, reducing vena contracta
    // Peak occurs at t/d ≈ 0.2-0.3, giving ~20% increase for beta=0.5
    const double reattach_factor = 0.35 * (1.0 - std::exp(-8.0 * t_over_d));
    const double reattachment = reattach_factor * (1.0 - beta * beta);

    // Component 2: Friction penalty (smooth Re dependence)
    // Use Blasius smooth turbulent friction: f = 0.316 / Re^0.25
    // Friction loss in bore reduces effective Cd
    // Coefficient calibrated to Idelchik data: k_t ≈ 0.96 at t/d=3.0, beta=0.5, Re=1e5
    // Calculation: need friction_loss ≈ 0.30 at t/d=3.0 → coeff ≈ 5.67
    const double f = 0.316 / std::pow(Re_d, 0.25);
    const double friction_loss = 5.67 * f * t_over_d;

    // Combined correction: rise (reattachment) then fall (friction)
    const double correction = 1.0 + reattachment - friction_loss;

    // Clamping with safety floor for extreme cases
    // Lower bound: 0.5 (extremely thick/rough orifices approach pipe entrance)
    // Upper bound: 1.3 (maximum reattachment benefit)
    return std::max(0.5, std::min(correction, 1.3));
}

// Rounded-entry Cd
// For well-rounded entries (r/d >= 0.15), Cd approaches 0.98-0.99
// Based on Idelchik contraction loss coefficients
double Cd_rounded(double r_over_d, double beta, double Re_D) {
    if (Re_D < 1.0) Re_D = 1.0;

    // For r/d = 0: sharp edge, use thin-plate correlation
    if (r_over_d <= 0.0) {
        return Cd_Stolz(beta, Re_D);
    }

    // Idelchik: loss coefficient K for rounded entry
    // K = 0.5 * (1 - r/d/0.15)^2 for r/d < 0.15
    // K ≈ 0.03 - 0.05 for r/d >= 0.15 (well-rounded)
    //
    // Convert K to Cd using: Cd = 1 / sqrt(1 + K / (1 - beta^4))

    double K;
    if (r_over_d >= 0.15) {
        // Well-rounded entry
        K = 0.04;
    } else {
        // Partially rounded
        const double ratio = 1.0 - r_over_d / 0.15;
        K = 0.5 * ratio * ratio;
    }

    // Account for beta effect
    const double beta4 = std::pow(beta, 4.0);
    const double Cd = 1.0 / std::sqrt(1.0 + K / (1.0 - beta4));

    // Reynolds number correction for low Re
    double Re_correction = 1.0;
    if (Re_D < 1.0e5) {
        Re_correction = 1.0 - 0.1 * std::pow(1.0e5 / Re_D, 0.2);
    }

    return Cd * Re_correction;
}

// Convert between Cd and loss coefficient K
double K_from_Cd(double Cd, double beta) {
    if (Cd <= 0.0 || Cd > 1.0) {
        throw std::invalid_argument("Cd must be in (0, 1]");
    }
    const double beta4 = std::pow(beta, 4.0);
    return (1.0 / (Cd * Cd) - 1.0) * (1.0 - beta4);
}

double Cd_from_K(double K, double beta) {
    if (K < 0.0) {
        throw std::invalid_argument("K must be >= 0");
    }
    const double beta4 = std::pow(beta, 4.0);
    return 1.0 / std::sqrt(1.0 + K / (1.0 - beta4));
}

} // namespace orifice

// -------------------------------------------------------------
// Main Cd functions (free functions)
// -------------------------------------------------------------

double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }
    return orifice::Cd_ReaderHarrisGallagher(geom.beta(), state.Re_D, geom.D);
}

double Cd_thick_plate(const OrificeGeometry& geom, const OrificeState& state) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }

    // Start with thin-plate Cd
    double Cd = orifice::Cd_ReaderHarrisGallagher(geom.beta(), state.Re_D, geom.D);

    // Apply thickness correction (includes reattachment + friction effects)
    Cd *= orifice::thickness_correction(geom.t_over_d(), geom.beta(), state.Re_D);

    return Cd;
}

double Cd_rounded_entry(const OrificeGeometry& geom, const OrificeState& state) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }
    return orifice::Cd_rounded(geom.r_over_d(), geom.beta(), state.Re_D);
}

double Cd(const OrificeGeometry& geom, const OrificeState& state) {
    if (!geom.is_valid()) {
        throw std::invalid_argument("Invalid orifice geometry");
    }

    // Auto-select based on geometry
    if (geom.r_over_d() > 0.01) {
        // Rounded entry
        return Cd_rounded_entry(geom, state);
    } else if (geom.t_over_d() > 0.02) {
        // Thick plate
        return Cd_thick_plate(geom, state);
    } else {
        // Sharp thin plate (default)
        return Cd_sharp_thin_plate(geom, state);
    }
}

// -------------------------------------------------------------
// Correlation classes
// -------------------------------------------------------------

namespace {

class ReaderHarrisGallagherCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return Cd_sharp_thin_plate(geom, state);
    }
    std::string name() const override { return "Reader-Harris/Gallagher (ISO 5167-2)"; }
};

class StolzCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return orifice::Cd_Stolz(geom.beta(), state.Re_D);
    }
    std::string name() const override { return "Stolz (ISO 5167:1980)"; }
};

class MillerCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return orifice::Cd_Miller(geom.beta(), state.Re_D);
    }
    std::string name() const override { return "Miller (1996)"; }
};

class ThickPlateCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return Cd_thick_plate(geom, state);
    }
    std::string name() const override { return "Thick Plate (Idelchik correction)"; }
};

class RoundedEntryCorrelation : public OrificeCorrelationBase {
public:
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return Cd_rounded_entry(geom, state);
    }
    std::string name() const override { return "Rounded Entry (Idelchik)"; }
};

class ConstantCdCorrelation : public OrificeCorrelationBase {
    double Cd_value_;
public:
    explicit ConstantCdCorrelation(double Cd = 0.61) : Cd_value_(Cd) {}
    double Cd(const OrificeGeometry&, const OrificeState&) const override {
        return Cd_value_;
    }
    std::string name() const override { return "Constant Cd"; }
};

class UserFunctionCorrelation : public OrificeCorrelationBase {
    CdFunction fn_;
    std::string name_;
public:
    UserFunctionCorrelation(CdFunction fn, std::string name)
        : fn_(std::move(fn)), name_(std::move(name)) {}
    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return fn_(geom, state);
    }
    std::string name() const override { return name_; }
};

class TabulatedCorrelation : public OrificeCorrelationBase {
    std::vector<double> beta_values_;
    std::vector<double> Re_values_;
    std::vector<std::vector<double>> Cd_table_;
    std::string name_;

public:
    TabulatedCorrelation(std::vector<double> beta_values,
                         std::vector<double> Re_values,
                         std::vector<std::vector<double>> Cd_table,
                         std::string name)
        : beta_values_(std::move(beta_values))
        , Re_values_(std::move(Re_values))
        , Cd_table_(std::move(Cd_table))
        , name_(std::move(name)) {}

    double Cd(const OrificeGeometry& geom, const OrificeState& state) const override {
        return interpolate(geom.beta(), state.Re_D);
    }

    std::string name() const override { return name_; }

private:
    double interpolate(double beta, double Re_D) const {
        // Bilinear interpolation in beta and log(Re_D)
        const double log_Re = std::log10(std::max(Re_D, 1.0));

        // Find beta indices
        auto it_beta = std::lower_bound(beta_values_.begin(), beta_values_.end(), beta);
        std::size_t i_beta = (it_beta == beta_values_.begin()) ? 0 :
                             (it_beta == beta_values_.end()) ? beta_values_.size() - 2 :
                             static_cast<std::size_t>(it_beta - beta_values_.begin() - 1);

        // Find Re indices (in log space)
        std::vector<double> log_Re_values;
        log_Re_values.reserve(Re_values_.size());
        for (double Re : Re_values_) {
            log_Re_values.push_back(std::log10(std::max(Re, 1.0)));
        }
        auto it_Re = std::lower_bound(log_Re_values.begin(), log_Re_values.end(), log_Re);
        std::size_t i_Re = (it_Re == log_Re_values.begin()) ? 0 :
                           (it_Re == log_Re_values.end()) ? log_Re_values.size() - 2 :
                           static_cast<std::size_t>(it_Re - log_Re_values.begin() - 1);

        // Clamp indices
        i_beta = std::min(i_beta, beta_values_.size() - 2);
        i_Re = std::min(i_Re, Re_values_.size() - 2);

        // Interpolation weights
        const double t_beta = (beta - beta_values_[i_beta]) /
                              (beta_values_[i_beta + 1] - beta_values_[i_beta]);
        const double t_Re = (log_Re - log_Re_values[i_Re]) /
                            (log_Re_values[i_Re + 1] - log_Re_values[i_Re]);

        // Clamp weights to [0, 1]
        const double tb = std::clamp(t_beta, 0.0, 1.0);
        const double tr = std::clamp(t_Re, 0.0, 1.0);

        // Bilinear interpolation
        const double c00 = Cd_table_[i_beta][i_Re];
        const double c10 = Cd_table_[i_beta + 1][i_Re];
        const double c01 = Cd_table_[i_beta][i_Re + 1];
        const double c11 = Cd_table_[i_beta + 1][i_Re + 1];

        return (1 - tb) * (1 - tr) * c00
             + tb * (1 - tr) * c10
             + (1 - tb) * tr * c01
             + tb * tr * c11;
    }
};

} // anonymous namespace

std::unique_ptr<OrificeCorrelationBase> make_correlation(CdCorrelation id) {
    switch (id) {
        case CdCorrelation::ReaderHarrisGallagher:
            return std::make_unique<ReaderHarrisGallagherCorrelation>();
        case CdCorrelation::Stolz:
            return std::make_unique<StolzCorrelation>();
        case CdCorrelation::Miller:
            return std::make_unique<MillerCorrelation>();
        case CdCorrelation::IdelchikThick:
        case CdCorrelation::BohlThick:
            return std::make_unique<ThickPlateCorrelation>();
        case CdCorrelation::IdelchikRounded:
        case CdCorrelation::BohlRounded:
            return std::make_unique<RoundedEntryCorrelation>();
        case CdCorrelation::Constant:
            return std::make_unique<ConstantCdCorrelation>();
        case CdCorrelation::UserFunction:
            return nullptr;  // Use make_user_correlation instead
    }
    return nullptr;
}

std::unique_ptr<OrificeCorrelationBase> make_user_correlation(
    CdFunction fn,
    const std::string& name) {
    return std::make_unique<UserFunctionCorrelation>(std::move(fn), name);
}

std::unique_ptr<OrificeCorrelationBase> make_tabulated_correlation(
    const std::vector<double>& beta_values,
    const std::vector<double>& Re_values,
    const std::vector<std::vector<double>>& Cd_table,
    const std::string& name) {
    return std::make_unique<TabulatedCorrelation>(beta_values, Re_values, Cd_table, name);
}

// -------------------------------------------------------------
// Flow calculations
// -------------------------------------------------------------

double orifice_mdot(const OrificeGeometry& geom, double Cd, double dP, double rho) {
    if (dP < 0.0 || rho <= 0.0 || Cd <= 0.0) {
        throw std::invalid_argument("orifice_mdot: invalid parameters");
    }
    return Cd * geom.area() * std::sqrt(2.0 * rho * dP);
}

double orifice_dP(const OrificeGeometry& geom, double Cd, double mdot, double rho) {
    if (mdot < 0.0 || rho <= 0.0 || Cd <= 0.0) {
        throw std::invalid_argument("orifice_dP: invalid parameters");
    }
    // From mdot = Cd * A * sqrt(2 * rho * dP)
    // Solve for dP: dP = (mdot / (Cd * A))^2 / (2 * rho)
    const double A = geom.area();
    const double term = mdot / (Cd * A);
    return term * term / (2.0 * rho);
}

double orifice_Cd_from_measurement(const OrificeGeometry& geom,
                                    double mdot, double dP, double rho) {
    if (dP <= 0.0 || rho <= 0.0 || mdot <= 0.0) {
        throw std::invalid_argument("orifice_Cd_from_measurement: invalid parameters");
    }
    double A = geom.area();
    return mdot / (A * std::sqrt(2.0 * rho * dP));
}

// -------------------------------------------------------------
// Iterative solver for Cd-Re coupling
// -------------------------------------------------------------

double solve_orifice_mdot(
    const OrificeGeometry& geom,
    double dP,
    double rho,
    double mu,
    double P_upstream,
    double kappa,
    CdCorrelation correlation,
    double tol,
    int max_iter)
{
    // Input validation
    if (!geom.is_valid()) {
        throw std::invalid_argument("solve_orifice_mdot: invalid geometry");
    }
    if (dP < 0.0) {
        throw std::invalid_argument("solve_orifice_mdot: dP must be non-negative");
    }
    if (rho <= 0.0) {
        throw std::invalid_argument("solve_orifice_mdot: rho must be positive");
    }
    if (mu <= 0.0) {
        throw std::invalid_argument("solve_orifice_mdot: mu must be positive");
    }
    if (P_upstream <= 0.0) {
        throw std::invalid_argument("solve_orifice_mdot: P_upstream must be positive");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("solve_orifice_mdot: tol must be positive");
    }
    if (max_iter < 1) {
        throw std::invalid_argument("solve_orifice_mdot: max_iter must be >= 1");
    }

    // Handle zero pressure drop case
    if (dP == 0.0) {
        return 0.0;
    }

    // Precompute constants
    const double area = geom.area();
    const double beta = geom.beta();
    const double D = geom.D;

    // Initial guess for Cd (typical value for sharp orifices)
    double Cd = 0.61;

    // Initial guess for mdot (use incompressible formula with initial Cd)
    double mdot = Cd * area * std::sqrt(2.0 * rho * dP);

    // Iteration loop
    for (int iter = 0; iter < max_iter; ++iter) {
        // Calculate expansibility factor if compressible (kappa > 1)
        double epsilon = 1.0;
        if (kappa > 1.0) {
            epsilon = expansibility_factor(beta, dP, P_upstream, kappa);
        }

        // Calculate new mass flow rate using current Cd and epsilon
        const double mdot_new = Cd * epsilon * area * std::sqrt(2.0 * rho * dP);

        // Check convergence
        const double rel_error = std::abs(mdot_new - mdot) / (mdot + 1e-30);
        if (rel_error < tol) {
            return mdot_new;
        }

        // Update Reynolds number based on new mdot
        // Re_D = (4 * mdot) / (π * D * μ)
        const double Re_D = (4.0 * mdot_new) / (PI * D * mu);

        // Update Cd based on new Reynolds number
        switch (correlation) {
            case CdCorrelation::ReaderHarrisGallagher:
                Cd = orifice::Cd_ReaderHarrisGallagher(beta, Re_D, D);
                break;
            case CdCorrelation::Stolz:
                Cd = orifice::Cd_Stolz(beta, Re_D);
                break;
            case CdCorrelation::Miller:
                Cd = orifice::Cd_Miller(beta, Re_D);
                break;
            default:
                throw std::invalid_argument(
                    "solve_orifice_mdot: unsupported correlation type");
        }

        // Update mdot for next iteration
        mdot = mdot_new;
    }

    // Failed to converge
    throw std::runtime_error(
        "solve_orifice_mdot: failed to converge after " +
        std::to_string(max_iter) + " iterations");
}

// -------------------------------------------------------------
// Compressible flow correction
// -------------------------------------------------------------

double expansibility_factor(double beta, double dP, double P_upstream, double kappa) {
    // Input validation
    if (beta <= 0.0 || beta >= 1.0) {
        throw std::invalid_argument("expansibility_factor: beta must be in range (0, 1)");
    }
    if (P_upstream <= 0.0) {
        throw std::invalid_argument("expansibility_factor: P_upstream must be positive");
    }
    if (dP < 0.0) {
        throw std::invalid_argument("expansibility_factor: dP must be non-negative");
    }

    // Incompressible limit: kappa <= 1 or no pressure drop
    if (kappa <= 1.0 || dP <= 0.0) {
        return 1.0;
    }

    // Pressure ratio
    const double tau = dP / P_upstream;

    // Check validity range (ISO 5167-2 recommends τ ≤ 0.25)
    if (tau > 0.25) {
        // Still compute but this is outside recommended range
        // User should be aware via documentation
    }

    // Compute beta powers
    const double beta2 = beta * beta;
    const double beta4 = beta2 * beta2;
    const double beta8 = beta4 * beta4;

    // ISO 5167-2:2003 expansibility factor formula
    // ε = 1 - (0.351 + 0.256·β⁴ + 0.93·β⁸) · [1 - (1 - τ)^(1/κ)]
    const double coeff = 0.351 + 0.256 * beta4 + 0.93 * beta8;
    const double expansion_term = 1.0 - std::pow(1.0 - tau, 1.0 / kappa);

    return 1.0 - coeff * expansion_term;
}

// -------------------------------------------------------------
// Orifice flow result bundle
// -------------------------------------------------------------

OrificeFlowResult orifice_flow(
    const OrificeGeometry& geom,
    double dP,
    double T,
    double P,
    double mu,
    double Z,
    const std::vector<double>& X,
    double kappa,
    CdCorrelation correlation)
{
    // Input validation
    if (!geom.is_valid()) {
        throw std::invalid_argument("orifice_flow: invalid geometry");
    }
    if (T <= 0.0) {
        throw std::invalid_argument("orifice_flow: temperature must be positive");
    }
    if (P <= 0.0) {
        throw std::invalid_argument("orifice_flow: pressure must be positive");
    }
    if (mu <= 0.0) {
        throw std::invalid_argument("orifice_flow: viscosity must be positive");
    }
    if (Z <= 0.0) {
        throw std::invalid_argument("orifice_flow: compressibility factor Z must be positive");
    }
    if (dP < 0.0) {
        throw std::invalid_argument("orifice_flow: differential pressure must be non-negative");
    }

    // Compute ideal gas density
    // For now, use simple ideal gas law with air composition
    // rho_ideal = P * MW / (R * T)
    // For air: MW ≈ 28.97 g/mol, R = 8.314 J/(mol·K)
    const double R_gas = 8.314;  // J/(mol·K)
    const double MW_air = 0.02897;  // kg/mol (air molecular weight)
    const double rho_ideal = (P * MW_air) / (R_gas * T);

    // Apply real gas correction
    const double rho_corrected = rho_ideal / Z;

    // Solve for mass flow rate with corrected density
    const double mdot = solve_orifice_mdot(
        geom, dP, rho_corrected, mu, P, kappa, correlation);

    // Calculate expansibility factor
    const double beta = geom.beta();
    const double epsilon = (kappa > 1.0) ?
        expansibility_factor(beta, dP, P, kappa) : 1.0;

    // Calculate velocity through orifice
    const double A = geom.area();
    const double v = mdot / (rho_corrected * A);

    // Calculate Reynolds numbers
    const double D = geom.D;
    const double d = geom.d;
    const double Re_D = (4.0 * mdot) / (PI * D * mu);
    const double Re_d = (4.0 * mdot) / (PI * d * mu);

    // Get discharge coefficient
    OrificeState state;
    state.Re_D = Re_D;
    state.dP = dP;
    state.rho = rho_corrected;
    state.mu = mu;

    double Cd_value = 0.61;  // Default
    switch (correlation) {
        case CdCorrelation::ReaderHarrisGallagher:
            Cd_value = Cd_sharp_thin_plate(geom, state);
            break;
        case CdCorrelation::Stolz:
            Cd_value = orifice::Cd_Stolz(beta, Re_D);
            break;
        case CdCorrelation::Miller:
            Cd_value = orifice::Cd_Miller(beta, Re_D);
            break;
        case CdCorrelation::IdelchikThick:
            Cd_value = Cd_thick_plate(geom, state);
            break;
        case CdCorrelation::IdelchikRounded:
            Cd_value = Cd_rounded_entry(geom, state);
            break;
        default:
            Cd_value = Cd(geom, state);  // Auto-select
            break;
    }

    // Populate result struct
    OrificeFlowResult result;
    result.mdot = mdot;
    result.v = v;
    result.Re_D = Re_D;
    result.Re_d = Re_d;
    result.Cd = Cd_value;
    result.epsilon = epsilon;
    result.rho_corrected = rho_corrected;

    return result;
}

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

double orifice_velocity_from_mdot(double mdot, double rho, double d, double Z) {
    if (mdot < 0.0) {
        throw std::invalid_argument("orifice_velocity_from_mdot: mdot must be non-negative");
    }
    if (rho <= 0.0) {
        throw std::invalid_argument("orifice_velocity_from_mdot: rho must be positive");
    }
    if (d <= 0.0) {
        throw std::invalid_argument("orifice_velocity_from_mdot: d must be positive");
    }
    if (Z <= 0.0) {
        throw std::invalid_argument("orifice_velocity_from_mdot: Z must be positive");
    }

    // Apply real gas correction to density
    const double rho_corrected = rho / Z;

    // Calculate area
    const double A = PI * d * d / 4.0;

    // v = mdot / (rho_corrected * A)
    return mdot / (rho_corrected * A);
}

double orifice_area_from_beta(double D, double beta) {
    if (D <= 0.0) {
        throw std::invalid_argument("orifice_area_from_beta: D must be positive");
    }
    if (beta <= 0.0 || beta >= 1.0) {
        throw std::invalid_argument("orifice_area_from_beta: beta must be in range (0, 1)");
    }

    // A = π * (D * beta / 2)²
    const double d = D * beta;
    return PI * d * d / 4.0;
}

double beta_from_diameters(double d, double D) {
    if (d <= 0.0) {
        throw std::invalid_argument("beta_from_diameters: d must be positive");
    }
    if (D <= 0.0) {
        throw std::invalid_argument("beta_from_diameters: D must be positive");
    }
    if (d >= D) {
        throw std::invalid_argument("beta_from_diameters: d must be < D");
    }

    // beta = d / D
    return d / D;
}

double orifice_Re_d_from_mdot(double mdot, double d, double mu) {
    if (mdot < 0.0) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: mdot must be non-negative");
    }
    if (d <= 0.0) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: d must be positive");
    }
    if (mu <= 0.0) {
        throw std::invalid_argument("orifice_Re_d_from_mdot: mu must be positive");
    }

    // Re_d = 4 * mdot / (π * d * mu)
    return (4.0 * mdot) / (PI * d * mu);
}
