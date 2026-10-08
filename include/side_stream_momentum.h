#pragma once

// Momentum of side streams joining a channel, compressible (#465, #471).
//
// A stream that enters a channel brings mass, and only the streamwise part of
// its momentum. The channel flow must make up the rest, which costs static
// pressure. Over one injection station of constant cross-section A, flow m_a
// arriving and m_b leaving, the impulse balance is exact for compressible
// flow:
//
//   P_a A + m_a u_a + J = P_b A + m_b u_b,    u = m / (rho A)
//
// with J the side streams' streamwise momentum (zero for normal injection).
// Density is each station's own (static P, T); nothing here assumes it
// constant. At low Mach it reduces to Florschuetz, Truman and Metzger's
// (1981) P + G^2/rho = const.

#include <vector>

namespace combaero {

// ---------------------------------------------------------------
// One station, half at a time: merge AND bleed (#471)
// ---------------------------------------------------------------
//
// A station is a node where a side stream joins (merge) or leaves (bleed) a
// channel of cross-section A; m_a arrives, m_b leaves along the channel.
// With the side stream's AXIAL velocity written as u_s = kappa * u_a,
// momentum over the station gives the static-pressure DROP
//
//   dP = (m_b|m_b| - m_a|m_a| + kappa m_a (m_a - m_b)) / (rho A^2)
//
// at the station node's own density.
//
//   kappa = 0     normal injection (merge): the jets bring no axial momentum.
//                 Florschuetz's P + G^2/rho = const; the #465 impingement term.
//   kappa = 0.75  bleed through the wall: Bassett, Winterbone & Pearson (2001)
//                 separating straight-run coefficient K2/K5, Eq. (15),
//                 K = q^2 - 1.5 q + 0.5 (q = m_b/m_a), theta- and
//                 psi-independent, rated excellent against their data. The
//                 station reproduces it for EVERY q exactly at kappa = 0.75:
//                 (Pt_a - Pt_b)/(rho u_a^2 / 2) = 2(1-q)(kappa - (1+q)/2).
//
// A network segment between two stations carries HALF of each, so every node
// sits at the middle of its station (the centred scheme of #465). The
// segment's own flow then enters both halves at their own densities, which
// keeps the term m^2 (1/rho_from - 1/rho_to)/(2 A^2) that a single-density
// form drops. One state per node makes the scheme consistent to
// O(M^2 dm/m), not exact.
constexpr double STATION_KAPPA_MERGE_NORMAL = 0.0;
constexpr double STATION_KAPPA_BLEED_BASSETT = 0.75;

struct StationHalfDrop {
  double dP = 0.0;      // half the station's static-pressure drop [Pa]
  double d_dm_a = 0.0;  // [Pa/(kg/s)]
  double d_dm_b = 0.0;  // [Pa/(kg/s)]
  double d_dP = 0.0;    // d/d(station static P) [-]
  double d_dT = 0.0;    // d/d(station T) [Pa/K]
};

StationHalfDrop station_half_drop(double m_a, double m_b, double P, double T,
                                  const std::vector<double>& X, double area,
                                  double kappa);

// Entry from a reservoir (Pt) into a duct, static at the entry node:
//   dP = (1 + K_in) m|m| / (2 rho A^2)
// the dynamic head plus the entry loss (K_in = 0.5 sharp, ~0 well rounded).
StationHalfDrop channel_entry_drop(double m, double P, double T,
                                   const std::vector<double>& X, double area,
                                   double K_in);

// ---------------------------------------------------------------
// Merge chamber: main inlet + side streams -> one outlet (#471)
// ---------------------------------------------------------------
//
// Constant-area control volume of area A, axis = the outlet direction. The
// chamber node carries the outlet state (P, T static, as its own stagnation
// closure takes it). The main stream arrives along the axis at the MAIN FACE,
// whose static temperature is the main stream's T_main. Impulse:
//
//   P_f + m_main^2 R T_main / (P_f A^2) = P + m_out^2 / (rho A^2) - J / A
//
// a quadratic in the face static pressure P_f for an ideal gas, solved in
// closed form on its SUBSONIC root. Pt_face = P0_from_static(P_f, T_main,
// M_f), the same entropy-based closure as the chamber itself (#357). No
// loss coefficient: the stagnation-pressure loss of mixing follows from
// momentum; transverse momentum is reacted by the walls. With no side streams
// (m_out = m_main, J = 0, T = T_main) the face IS the chamber.
//
// When the impulse cannot be carried subsonically (the main stream would
// have to choke at the face) the discriminant is held at a floor and
// `choked` is set: the result is then a bound, not a solution.
struct MergeFaceState {
  double P_face = 0.0;   // [Pa]
  double Pt_face = 0.0;  // [Pa]
  double M_face = 0.0;   // [-]
  bool choked = false;
  // d/d(m_main, m_out, J, P, T, T_main). Flows and pressures analytic; the
  // two temperatures by central difference inside C++, as the chamber's own
  // stagnation closure does (smooth, no inner discretisation).
  double dPf_dm_main = 0.0, dPf_dm_out = 0.0, dPf_dJ = 0.0;
  double dPf_dP = 0.0, dPf_dT = 0.0, dPf_dT_main = 0.0;
  double dPtf_dm_main = 0.0, dPtf_dm_out = 0.0, dPtf_dJ = 0.0;
  double dPtf_dP = 0.0, dPtf_dT = 0.0, dPtf_dT_main = 0.0;
};

MergeFaceState chamber_merge_face_state(double m_main, double T_main,
                                        const std::vector<double>& X_main,
                                        double m_out, double P, double T,
                                        const std::vector<double>& X,
                                        double side_momentum, double area);

// Streamwise impulse of a jet discharged from stagnation (Pt, Tt) into static
// P, per unit cos(theta): J = m w.
//   unchoked: w = u_is, the isentropic velocity to P (variable cp);
//   choked (P < P*): w = u* + (P* - P) / (rho* u*), the sonic jet's momentum
//                    plus its pressure thrust on the vena-contracta area.
// Cd does not enter: it sets the area, not the jet velocity. w is C1 across
// choking. P >= Pt gives w = 0 (no outflow to carry momentum).
struct JetImpulse {
  double J = 0.0;       // [N]
  double dJ_dm = 0.0;   // = w [m/s]
  double dJ_dPt = 0.0, dJ_dTt = 0.0, dJ_dP = 0.0;
  bool choked = false;
};

JetImpulse jet_impulse(double m, double Pt, double Tt, double P,
                       const std::vector<double>& X);

}  // namespace combaero
