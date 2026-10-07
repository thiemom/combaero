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
// Crossflow segment, centred (impingement crossflow, #465)
// ---------------------------------------------------------------
//
// A network segment runs between two stations. Splitting each station's
// drop half-and-half between the segments either side puts every node at the
// MIDDLE of its station, where a continuous model evaluates the pressure the
// side stream sees. With normal injection (J = 0) one segment carries
//
//   dP = (m_out|m_out| / rho_out - m_arr|m_arr| / rho_arr) / (2 A^2)
//
// m_arr arriving at its upstream station, m_out leaving its downstream
// station, each density from that station's own static (P, T, X). Putting the whole drop
// downstream instead over-fed the downstream rows: -13% in Gc/Gj at row 10
// of Florschuetz's strongest-crossflow geometry, against -2.5% centred.
struct SideStreamMomentum {
  double dP = 0.0;          // static-pressure drop along the segment [Pa]
  double d_dm_out = 0.0;    // [Pa/(kg/s)]
  double d_dm_arr = 0.0;    // [Pa/(kg/s)]
  // With respect to each station's own static state, through its density
  // (C++'s density_and_jacobians, the mixture's own equation of state).
  double d_dP_arr = 0.0, d_dT_arr = 0.0;   // [-], [Pa/K]
  double d_dP_out = 0.0, d_dT_out = 0.0;   // [-], [Pa/K]
};

SideStreamMomentum side_stream_momentum_drop(
    double m_arr, double m_out, double P_arr, double T_arr,
    const std::vector<double>& X_arr, double P_out, double T_out,
    const std::vector<double>& X_out, double area);

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
