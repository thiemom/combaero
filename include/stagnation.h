#ifndef STAGNATION_H
#define STAGNATION_H

#include <cstddef>
#include <vector>

namespace combaero {

// -------------------------------------------------------------
// Isentropic stagnation / static conversion utilities
// -------------------------------------------------------------
// These functions convert between static and stagnation (total) conditions
// for a thermally perfect gas using the variable-cp thermo backend (thermo.h).
//
// Conventions used throughout CombAero:
//   - thermo.h functions h(T,X), cp(T,X), s(T,P,X) operate on STATIC conditions.
//   - Stagnation enthalpy:  h0 = h(T_static) + v^2/2  [J/kg]
//   - Stagnation temperature T0: defined by h(T0) = h0  (iterative for variable cp)
//   - Stagnation pressure P0:    isentropic relation s(T0,P0) = s(T,P)
//
// Adiabatic wall temperature (for convective heat transfer at high Mach):
//   q = h_conv * (T_aw - T_wall)          <- correct driving temperature
//   T_aw = T_static + r * v^2 / (2*cp)
//   r = Pr^(1/3)  turbulent  (default)
//   r = Pr^(1/2)  laminar
//
//   T_aw is always between T_static and T_total:
//     T_static < T_aw < T_total  (since r < 1 for air)
//   Using T_static underestimates heat load at high M.
//   Using T_total overestimates (r < 1, so T_aw < T0 always).
//
// References:
//   Anderson, J.D. (2003) Modern Compressible Flow, 3rd ed., Ch. 3
//   Incropera et al. (2007) Fundamentals of Heat and Mass Transfer, 6th ed.

// -------------------------------------------------------------
// Kinetic energy and bulk velocity
// -------------------------------------------------------------

// Specific kinetic energy: v^2/2  [J/kg]
inline double kinetic_energy(double v) { return 0.5 * v * v; }

// Bulk velocity from mass flow, density, and cross-sectional area.
// v = m_dot / (rho * area)  [m/s]
inline double bulk_velocity(double m_dot, double rho, double area) {
    return m_dot / (rho * area);
}

// -------------------------------------------------------------
// Stagnation enthalpy (exact, no iteration)
// -------------------------------------------------------------

// Stagnation enthalpy from static enthalpy and velocity.
// h0 = h_static + v^2/2  [J/kg]
// h_static must be in [J/kg] (use h_mass() from thermo.h).
double h0_from_static(double h_static_J_per_kg, double v);

// Velocity from stagnation and static enthalpy.
// v = sqrt(2*(h0 - h_static))  [m/s]
// Returns 0 if h0 <= h_static.
double v_from_h0(double h0, double h_static_J_per_kg);

// -------------------------------------------------------------
// Mach number
// -------------------------------------------------------------

// Mach number from velocity and static state.
// M = v / a(T, X)
double mach_number(double v, double T, const std::vector<double>& X);

// -------------------------------------------------------------
// Stagnation temperature (iterative — variable cp)
// -------------------------------------------------------------

// Stagnation temperature T0 from static T and Mach number.
// Solves h(T0) = h(T) + v^2/2 iteratively using Newton on h(T).
// v = M * a(T, X)  is computed internally.
//
// Parameters:
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation temperature T0 [K]
double T0_from_static(double T, double M, const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: T0 from static T and velocity v [m/s]
double T0_from_static_v(double T, double v, const std::vector<double>& X,
                        double tol = 1e-8, std::size_t max_iter = 50);

// Static temperature from stagnation T0 and Mach number.
// Solves h(T0) = h(T) + M^2*a(T)^2/2 iteratively.
// Returns: static temperature T [K]
double T_from_stagnation(double T0, double M, const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Stagnation pressure (isentropic — uses entropy conservation)
// -------------------------------------------------------------

// Stagnation pressure P0 from static (P, T) and Mach number.
// Uses isentropic relation: s(T0, P0) = s(T, P)
// Solves for P0 given T0 = T0_from_static(T, M, X).
//
// Parameters:
//   P   : static pressure [Pa]
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation pressure P0 [Pa]
double P0_from_static(double P, double T, double M,
                      const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Static pressure from stagnation (P0, T0) and Mach number.
// Uses isentropic relation: s(T, P) = s(T0, P0)
// Returns: static pressure P [Pa]
double P_from_stagnation(double P0, double T0, double M,
                         const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Adiabatic wall temperature
// -------------------------------------------------------------
// For convective heat transfer, the correct hot-side driving temperature is
// T_aw, not T_static or T_total.  Pass T_aw as T_hot to:
//   heat_transfer.h : htc_pipe(), nusselt_pipe()
//   cooling_correlations.h : cooled_wall_heat_flux()
//
// For M < 0.3: T_aw ~ T_static (< 1% error, correction negligible).
// For M = 0.3-0.8 (combustor liner, turbine cooling): correction is 2-15%.

// Adiabatic wall temperature from static temperature and velocity.
// T_aw = T_static + r * v^2 / (2 * cp(T, X))
// r = Pr^(1/3) for turbulent flow (default), Pr^(1/2) for laminar.
//
// Parameters:
//   T_static  : bulk static temperature [K]
//   v         : bulk flow velocity [m/s]
//   T         : temperature for property evaluation [K] (use T_static)
//   P         : pressure [Pa]
//   X         : mole fractions [-]
//   turbulent : true = turbulent (r = Pr^1/3), false = laminar (r = Pr^1/2)
// Returns: adiabatic wall temperature T_aw [K]
double T_adiabatic_wall(double T_static, double v,
                        double T, double P, const std::vector<double>& X,
                        bool turbulent = true);

// Combined stagnation state from static conditions and velocity.
// Returns (Tt, Pt) = (T0, P0) in a single call.
// Equivalent to calling T0_from_static_v then P0_from_static, but avoids
// recomputing Mach number internally.
std::tuple<double, double> stagnation_from_static(
    double T, double P, double v,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: T_aw from static temperature and Mach number.
// Computes v = M * a(T, X) internally.
double T_adiabatic_wall_mach(double T_static, double M,
                              double T, double P, const std::vector<double>& X,
                              bool turbulent = true);

// Recovery factor r [-]
// r = Pr^(1/3) for turbulent, Pr^(1/2) for laminar
double recovery_factor(double Pr, bool turbulent = true);

// -------------------------------------------------------------
// Static state of a duct station from stagnation T0 and mass flux (#576)
// -------------------------------------------------------------
// Static state of a chamber from its STAGNATION temperature T0, static P and
// mass flux G = |m_dot| / A (issue #576). Energy: h(T0) = h(T) + v^2 / 2 with
// v = G / rho(T, P) = G R T / P, i.e.
//   f(T) = h(T0) - h(T) - (G R T / P)^2 / 2 = 0.
// f is strictly decreasing (f' = -cp - G^2 R^2 T / P^2) and f(T0) <= 0, so
// the root is unique and lies in (0, T0]; with cp rising in T, f is concave
// and Newton from T0 approaches it monotonically from above. The root exists
// at any G -- subsonic or supersonic -- so the closure itself no longer caps
// the flow; the choking limit comes out of P0(P) at fixed G as it should.
//
// Returns T and the ln(P0 / P) of the isentropic stagnation state, with their
// partial derivatives with respect to T0, P and G (closed form, chain rule
// through the implicit f).
struct ChamberStaticState {
  double T;
  double dT_dT0, dT_dP, dT_dG;
  double lnP0P;
  double dln_dT0, dln_dP, dln_dG;
  double M;
};

ChamberStaticState chamber_static_state(double T0, double P, double G,
                                        const std::vector<double>& X);

// chamber_static_state's gradients in the MOLE fractions at fixed T0, P and
// G (#743), X unnormalised (map with mole_gradient_to_mass): the static T from
// f = h(T0) - h(T) - (G R T / P)^2 / 2 = 0 by the implicit function theorem,
// and ln(P0 / P) = sum_j X_j (s_j(T0) - s_j(T)) / R_u at that T.
struct ChamberStaticStateDX {
  std::vector<double> dT_dX, dlnP0P_dX;
};
ChamberStaticStateDX chamber_static_state_dX(double T0, double P, double G,
                                            const std::vector<double>& X);

// -------------------------------------------------------------
// The SONIC static state of a flow at stagnation temperature T0 and mass flux
// G (#637): T* from h(T0) - h(T*) = a(T*)^2 / 2 (variable cp, a^2 = gamma R T),
// and the static pressure at which G is sonic, P* = G sqrt(R T* / gamma*)
// (rho* a* = G with rho* = P* / (R T*)). A chamber of that flux cannot sit below P*: a
// duct discharging into a reservoir below it leaves choked, under-expanded,
// at P*. Derivatives in T0 and G are closed form (implicit function on T*).
// G = 0 gives P* = 0.
struct SonicStaticState {
  double T;
  double P;
  double dT_dT0;
  double dP_dT0;
  double dP_dG;
};

SonicStaticState sonic_static_state(double T0, double G, const std::vector<double>& X);

// -------------------------------------------------------------
// Static state at a duct's inlet from the STAGNATION state (T0, P0) it is
// fed from and its mass flux G = |m_dot| / A (#636): the isentropic
// acceleration from rest, i.e. chamber_static_state at the static P whose
// isentropic stagnation pressure is P0,
//   L(ln P) = ln(P0 / P)|isentropic(T0, P, G) - ln(P0 / P) = 0.
// L is monotone on the subsonic branch [P*, P0], P* the sonic pressure, so a
// bracketed Newton finds the unique subsonic root. A G above the choked
// mass flux has no subsonic state: the sonic state at P* is returned, with
// choked = true and no G dependence. A heat-transfer correlation needs this
// state: a plenum's T is the stagnation T, and reading it as the duct's
// static T put T_aw above T0 by about (1 + r) (T0 - T).
// Partials with respect to T0, P0 and G, by the implicit function theorem.
struct DuctStaticState {
  double T, P, rho, M;
  double dT_dT0, dT_dP0, dT_dG;
  double dP_dT0, dP_dP0, dP_dG;
  double drho_dT0, drho_dP0, drho_dG;
  bool choked;
  // The mass flux the state was evaluated at, and its partials: G itself
  // (dGeff_dG = 1) from duct_static_state; the bent flux from
  // duct_htc_static_state.
  double G_eff, dGeff_dT0, dGeff_dP0, dGeff_dG;
};

DuctStaticState duct_static_state(double T0, double P0, double G,
                                  const std::vector<double>& X);

// The duct state a heat-transfer correlation reads (#551): duct_static_state
// at G_eff, which is G itself up to (1 - 2w) G* and bends C1 onto (1 - w) G*
// above it (w = kDuctHtcSonicMargin, G* the sonic flux), so M stays below
// ~0.97. At G* the isentropic dM/dG is infinite and the plain state has a
// kink there; a duct choked at its INLET -- possible once cooling is in the
// momentum (Fanno-Rayleigh) -- sits exactly on it, and its wall coupling then
// could not converge. Below (1 - 2w) G* it is duct_static_state exactly.
// The partials include G_eff's dependence on (T0, P0) through G*.
constexpr double kDuctHtcSonicMargin = 1e-3;

DuctStaticState duct_htc_static_state(double T0, double P0, double G,
                                      const std::vector<double>& X);

// duct_htc_static_state's gradients in the MASS fractions at fixed T0, P0
// and G (#717): the static T, P, rho and the bent flux G_eff. The isentropic
// energy and entropy balances,
//   f = h(T0) - h(T) - (G R T / P)^2 / 2 = 0,
//   L = [sum_j X_j (s_j(T0) - s_j(T))] / R_u - ln(P0 / P) = 0,
// differentiated in X at the root (implicit function theorem; h, R and the
// species s_j are closed form). Near choke G_eff = ge(G / G*) G* also moves
// with G*(X), differenced like its T0 slope (smooth there). Zero for a
// floored T0 or P0.
struct DuctStaticStateDY {
  std::vector<double> dT_dY, dP_dY, drho_dY, dGeff_dY;
};

// `state` is duct_htc_static_state(T0, P0, G, X), already solved by the
// caller: below the bend (state.G_eff == G) nothing is solved again.
DuctStaticStateDY duct_htc_static_state_dY(double T0, double P0, double G,
                                           const std::vector<double>& X,
                                           const DuctStaticState& state);

} // namespace combaero

#endif // STAGNATION_H
