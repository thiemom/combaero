#ifndef STAGNATION_H
#define STAGNATION_H

#include <cstddef>
#include <tuple>
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
//   - Stagnation temperature Tt: defined by h(Tt) = h0  (iterative for variable cp)
//   - Stagnation pressure Pt:    isentropic relation s(Tt,Pt) = s(T,P)
//
// Adiabatic wall temperature (for convective heat transfer at high Mach):
//   q = h_conv * (T_aw - T_wall)          <- correct driving temperature
//   T_aw = T_static + r * v^2 / (2*cp)
//   r = Pr^(1/3)  turbulent  (default)
//   r = Pr^(1/2)  laminar
//
//   T_aw is always between T_static and Tt:
//     T_static < T_aw < Tt  (since r < 1 for air)
//   Using T_static underestimates heat load at high M.
//   Using Tt overestimates (r < 1, so T_aw < Tt always).
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
// Stagnation temperature (iterative -- variable cp)
// -------------------------------------------------------------

// Stagnation temperature Tt from static T and Mach number.
// Solves h(Tt) = h(T) + v^2/2 iteratively using Newton on h(T).
// v = M * a(T, X)  is computed internally.
//
// Parameters:
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation temperature Tt [K]
double Tt_from_static(double T, double M, const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Convenience: Tt from static T and velocity v [m/s]
double Tt_from_static_v(double T, double v, const std::vector<double>& X,
                        double tol = 1e-8, std::size_t max_iter = 50);

// Static temperature from stagnation Tt and Mach number.
// Solves h(Tt) = h(T) + M^2*a(T)^2/2 iteratively.
// Returns: static temperature T [K]
double T_from_stagnation(double Tt, double M, const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Stagnation pressure (isentropic -- uses entropy conservation)
// -------------------------------------------------------------

// Stagnation pressure Pt from static (P, T) and Mach number.
// Uses isentropic relation: s(Tt, Pt) = s(T, P)
// Solves for Pt given Tt = Tt_from_static(T, M, X).
//
// Parameters:
//   P   : static pressure [Pa]
//   T   : static temperature [K]
//   M   : Mach number [-]
//   X   : mole fractions [-]
// Returns: stagnation pressure Pt [Pa]
double Pt_from_static(double P, double T, double M,
                      const std::vector<double>& X,
                      double tol = 1e-8, std::size_t max_iter = 50);

// Static pressure from stagnation (Pt, Tt) and Mach number.
// Uses isentropic relation: s(T, P) = s(Tt, Pt)
// Returns: static pressure P [Pa]
double P_from_stagnation(double Pt, double Tt, double M,
                         const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);

// Deprecated names (#782): the T0/P0 spelling of the functions above. They
// forward unchanged and are removed in 0.9.0.
[[deprecated("renamed Tt_from_static (#782); removed in 0.9.0")]]
inline double T0_from_static(double T, double M, const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50) {
  return Tt_from_static(T, M, X, tol, max_iter);
}
[[deprecated("renamed Tt_from_static_v (#782); removed in 0.9.0")]]
inline double T0_from_static_v(double T, double v, const std::vector<double>& X,
                               double tol = 1e-8, std::size_t max_iter = 50) {
  return Tt_from_static_v(T, v, X, tol, max_iter);
}
[[deprecated("renamed Pt_from_static (#782); removed in 0.9.0")]]
inline double P0_from_static(double P, double T, double M, const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50) {
  return Pt_from_static(P, T, M, X, tol, max_iter);
}

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
// Returns (Tt, Pt) in a single call.
// Equivalent to calling Tt_from_static_v then Pt_from_static, but avoids
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
// Static state of a duct station from stagnation Tt and mass flux (#576)
// -------------------------------------------------------------
// Static state of a chamber from its STAGNATION temperature Tt, static P and
// mass flux G = |m_dot| / A (issue #576). Energy: h(Tt) = h(T) + v^2 / 2 with
// v = G / rho(T, P) = G R T / P, i.e.
//   f(T) = h(Tt) - h(T) - (G R T / P)^2 / 2 = 0.
// f is strictly decreasing (f' = -cp - G^2 R^2 T / P^2) and f(Tt) <= 0, so
// the root is unique and lies in (0, Tt]; with cp rising in T, f is concave
// and Newton from Tt approaches it monotonically from above. The root exists
// at any G -- subsonic or supersonic -- so the closure itself no longer caps
// the flow; the choking limit comes out of Pt(P) at fixed G as it should.
//
// Returns T and lnPtP = ln(Pt / P) of the isentropic stagnation state, with
// their partial derivatives with respect to Tt, P and G (closed form, chain
// rule through the implicit f).
struct ChamberStaticState {
  double T;
  double dT_dTt, dT_dP, dT_dG;
  double lnPtP;
  double dln_dTt, dln_dP, dln_dG;
  double M;
};

ChamberStaticState chamber_static_state(double Tt, double P, double G,
                                        const std::vector<double>& X);

// chamber_static_state's gradients in the MOLE fractions at fixed Tt, P and
// G (#743), X unnormalised (map with mole_gradient_to_mass): the static T from
// f = h(Tt) - h(T) - (G R T / P)^2 / 2 = 0 by the implicit function theorem,
// and ln(Pt / P) = sum_j X_j (s_j(Tt) - s_j(T)) / R_u at that T.
struct ChamberStaticStateDX {
  std::vector<double> dT_dX, dlnPtP_dX;
};
ChamberStaticStateDX chamber_static_state_dX(double Tt, double P, double G,
                                            const std::vector<double>& X);

// -------------------------------------------------------------
// The SONIC static state of a flow at stagnation temperature Tt and mass flux
// G (#637): T* from h(Tt) - h(T*) = a(T*)^2 / 2 (variable cp, a^2 = gamma R T),
// and the static pressure at which G is sonic, P* = G sqrt(R T* / gamma*)
// (rho* a* = G with rho* = P* / (R T*)). A chamber of that flux cannot sit below P*: a
// duct discharging into a reservoir below it leaves choked, under-expanded,
// at P*. Derivatives in Tt and G are closed form (implicit function on T*).
// G = 0 gives P* = 0.
struct SonicStaticState {
  double T;
  double P;
  double dT_dTt;
  double dP_dTt;
  double dP_dG;
};

SonicStaticState sonic_static_state(double Tt, double G, const std::vector<double>& X);

// -------------------------------------------------------------
// Static state at a duct's inlet from the STAGNATION state (Tt, Pt) it is
// fed from and its mass flux G = |m_dot| / A (#636): the isentropic
// acceleration from rest, i.e. chamber_static_state at the static P whose
// isentropic stagnation pressure is Pt,
//   L(ln P) = ln(Pt / P)|isentropic(Tt, P, G) - ln(Pt / P) = 0.
// L is monotone on the subsonic branch [P*, Pt], P* the sonic pressure, so a
// bracketed Newton finds the unique subsonic root. A G above the choked
// mass flux has no subsonic state: the sonic state at P* is returned, with
// choked = true and no G dependence. A heat-transfer correlation needs this
// state: a plenum's T is the stagnation T, and reading it as the duct's
// static T put T_aw above Tt by about (1 + r) (Tt - T).
// Partials with respect to Tt, Pt and G, by the implicit function theorem.
struct DuctStaticState {
  double T, P, rho, M;
  double dT_dTt, dT_dPt, dT_dG;
  double dP_dTt, dP_dPt, dP_dG;
  double drho_dTt, drho_dPt, drho_dG;
  bool choked;
  // The mass flux the state was evaluated at, and its partials: G itself
  // (dGeff_dG = 1) from duct_static_state; the bent flux from
  // duct_htc_static_state.
  double G_eff, dGeff_dTt, dGeff_dPt, dGeff_dG;
};

DuctStaticState duct_static_state(double Tt, double Pt, double G,
                                  const std::vector<double>& X);

// The duct state a heat-transfer correlation reads (#551): duct_static_state
// at G_eff, which is G itself up to (1 - 2w) G* and bends C1 onto (1 - w) G*
// above it (w = kDuctHtcSonicMargin, G* the sonic flux), so M stays below
// ~0.97. At G* the isentropic dM/dG is infinite and the plain state has a
// kink there; a duct choked at its INLET -- possible once cooling is in the
// momentum (Fanno-Rayleigh) -- sits exactly on it, and its wall coupling then
// could not converge. Below (1 - 2w) G* it is duct_static_state exactly.
// The partials include G_eff's dependence on (Tt, Pt) through G*.
constexpr double kDuctHtcSonicMargin = 1e-3;

DuctStaticState duct_htc_static_state(double Tt, double Pt, double G,
                                      const std::vector<double>& X);

// duct_htc_static_state's gradients in the MASS fractions at fixed Tt, Pt
// and G (#717): the static T, P, rho and the bent flux G_eff. The isentropic
// energy and entropy balances,
//   f = h(Tt) - h(T) - (G R T / P)^2 / 2 = 0,
//   L = [sum_j X_j (s_j(Tt) - s_j(T))] / R_u - ln(Pt / P) = 0,
// differentiated in X at the root (implicit function theorem; h, R and the
// species s_j are closed form). Near choke G_eff = ge(G / G*) G* also moves
// with G*(X), differenced like its Tt slope (smooth there). Zero for a
// floored Tt or Pt.
struct DuctStaticStateDY {
  std::vector<double> dT_dY, dP_dY, drho_dY, dGeff_dY;
};

// `state` is duct_htc_static_state(Tt, Pt, G, X), already solved by the
// caller: below the bend (state.G_eff == G) nothing is solved again.
DuctStaticStateDY duct_htc_static_state_dY(double Tt, double Pt, double G,
                                           const std::vector<double>& X,
                                           const DuctStaticState& state);

} // namespace combaero

#endif // STAGNATION_H
