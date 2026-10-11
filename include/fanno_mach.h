#pragma once

#include <string>
#include <vector>

namespace combaero {

// Fanno flow integrated in Mach number.
//
// Adiabatic, constant-area duct flow with wall friction. At a fixed mass flux
// G = m_dot / A and stagnation temperature Tt (h0 = h(Tt), conserved), the
// state at a station is ALGEBRAIC in its Mach number M:
//
//   energy      h(T) + M^2 a(T)^2 / 2 = h0      -> T(M)
//   kinematics  u = M a(T),  rho = G / u,  P = rho R T
//
// so the only differential relation left is momentum,
//
//   dP + G du = -(f / 2D) G u dx   ->   dx/dM = -(P' + G u') 2D / (f G u),
//
// and f depends on the station only through Re = G D / mu(T(M)). The length
// between two Mach numbers is therefore a plain integral. Unlike marching in
// x, which runs into the 1 / (1 - M^2) singularity at sonic, dx/dM is finite
// and smooth everywhere and goes to zero at M = 1: the choking length
// L* = integral up to M = 1 is exact, with no cutoff Mach, no gradient floor
// and no step control. Thermally perfect gas (cp(T)), frozen composition.

struct FannoMachState {
  double M = 0.0;    // Mach number [-]
  double T = 0.0;    // static temperature [K]
  double P = 0.0;    // static pressure [Pa]
  double rho = 0.0;  // density [kg/m^3]
  double u = 0.0;    // velocity [m/s]
  double a = 0.0;    // speed of sound [m/s]
};

struct FannoDuctResult {
  bool choked = false;   // L >= L*: the duct cannot pass G subsonically
  double G = 0.0;        // mass flux [kg/(m^2 s)]
  double L_star = 0.0;   // inlet to sonic [m]
  FannoMachState inlet;  // static state at the inlet face
  FannoMachState exit;   // static state at L (at L* when choked)
  double Pt_exit = 0.0;  // stagnation pressure at the exit station [Pa]
};

// Relative tolerance of the Mach and flux roots. The length quadrature is a
// fixed rule (see fanno_mach.cpp), accurate to ~1e-10 at constant f and to the
// ~1e-8 floor that kinks in mu(T) set once f depends on Re.
constexpr double kFannoMachRelTol = 1e-10;

// Station state at Mach M for mass flux G and stagnation temperature Tt.
// P is rho R T with rho = G / u; at M = 0 the flow is at rest (u = 0) and P is
// undefined, so 0 < M <= 1 is required.
FannoMachState fanno_state_at_mach(double G, double Tt, double M,
                                   const std::vector<double>& X);

// dx/dM at Mach M [m]. friction_model is any friction_and_jacobian tag, or
// "fixed" for a constant Darcy f = f_multiplier.
double fanno_dx_dmach(double G, double Tt, double M,
                      const std::vector<double>& X, double D,
                      double roughness, const std::string& friction_model,
                      double f_multiplier = 1.0);

// Duct length [m] over which the flow goes from M1 to M2 (0 < M1 <= M2 <= 1).
double fanno_length_between(double G, double Tt, double M1, double M2,
                            const std::vector<double>& X, double D,
                            double roughness, const std::string& friction_model,
                            double f_multiplier = 1.0);

// Subsonic Mach number at which a stagnation state (Pt, Tt) accelerated
// isentropically carries mass flux G. Returns 1 when G is at or above the
// sonic flux (the inlet itself chokes).
double fanno_inlet_mach(double Pt, double Tt, double G,
                        const std::vector<double>& X);

// Sonic (maximum) isentropic mass flux of the stagnation state [kg/(m^2 s)].
double fanno_sonic_mass_flux(double Pt, double Tt, const std::vector<double>& X);

// Duct of length L fed isentropically from (Pt, Tt) at mass flux G.
FannoDuctResult fanno_duct(double Pt, double Tt, double G,
                           const std::vector<double>& X, double L, double D,
                           double roughness, const std::string& friction_model,
                           double f_multiplier = 1.0);

// The mass flux at which a duct of length L chokes exactly at its exit
// (L* = L), fed isentropically from (Pt, Tt) [kg/(m^2 s)].
double fanno_choked_mass_flux(double Pt, double Tt,
                              const std::vector<double>& X, double L, double D,
                              double roughness,
                              const std::string& friction_model,
                              double f_multiplier = 1.0);

// Mass flux a duct passes from a stagnation state (Pt0, Tt0) against a
// downstream pressure, with its derivatives. The physical inverse of
// fanno_duct, and the form a network residual wants: it rises monotonically
// as P_target falls and SATURATES at the choked flux once the exit chokes,
// so it exists for every P_target (no infeasible region, no barrier).
//
// exit_total selects what P_target is matched against:
//   true  -- the exit STAGNATION pressure (the exit dynamic head is recovered
//            downstream: a diffuser, or a node that carries its own Pt);
//   false -- the exit STATIC pressure (the head is lost in a sudden expansion,
//            a bare exit into a plenum or the atmosphere).
//
// Below kFannoFlowDriveFloor of drive (Pt0 - P_target) the flux law, whose
// slope is infinite at zero drive, is blended C1 into an odd cubic through
// zero that matches its value and slope at the floor; at or below zero drive the flux is 0 (the caller reverses
// the duct for reversed flow).
struct FannoChannelFlow {
  double G = 0.0;              // mass flux [kg/(m^2 s)], >= 0
  bool choked = false;         // exit at sonic, flux independent of P_target
  double M_in = 0.0;           // inlet Mach
  double M_exit = 0.0;         // exit Mach (1 when choked)
  double dG_dPt0 = 0.0;        // [kg/(m^2 s Pa)]
  double dG_dTt0 = 0.0;        // [kg/(m^2 s K)]
  double dG_dP_target = 0.0;   // [kg/(m^2 s Pa)], 0 when choked
  // dG/ds [kg/(m^2 s)] along each of dY_directions (Y -> Y + s v, X
  // renormalised); empty when none are passed. Chained by a network element
  // through an orthonormal basis of its feeding node's dY/dx (see
  // fanno_mach.cpp d_ddir); a unit vector e_k gives dG/dY_k at fixed others.
  std::vector<double> dG_ddir;
};

// Drive [Pa] below which the flux is blended to a finite slope at zero.
constexpr double kFannoFlowDriveFloor = 1.0;

FannoChannelFlow fanno_channel_flow(double Pt0, double Tt0, const std::vector<double>& X,
                                    double P_target, bool exit_total, double L, double D,
                                    double roughness, const std::string& friction_model,
                                    double f_multiplier = 1.0,
                                    // Previous exit Mach, if any: a warm start
                                    // for the root (the result does not depend
                                    // on it, only the cost).
                                    double M_exit_guess = -1.0,
                                    // Mass-fraction directions for dG_ddir
                                    // (about two length integrals each).
                                    const std::vector<std::vector<double>>& dY_directions = {},
                                    // false: G, choked and the Mach numbers
                                    // only -- the derivative fields are not
                                    // computed (a residual-only evaluation).
                                    bool with_derivatives = true);

}  // namespace combaero
