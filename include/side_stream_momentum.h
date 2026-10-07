#pragma once

// Momentum of side streams joining a channel with no streamwise momentum.
//
// A stream that enters a channel normal to its axis (impingement jets turning
// into the crossflow, #465; any perpendicular injection) brings mass but no
// streamwise momentum, so the channel flow must accelerate it. A momentum
// balance over one injection station, flow m_a arriving and m_b leaving
// through cross-section A, gives the static-pressure drop
//
//   dP_station = (m_b |m_b| - m_a |m_a|) / (rho A^2)
//
// -- the discrete form of P + G^2/rho = const behind Florschuetz, Truman and
// Metzger's (1981) jet-array flow distribution (their Eqs. 2-6).
//
// CENTRED SEGMENT FORM. A network channel segment runs between two stations.
// Splitting each station's drop half-and-half between the segments either side
// puts every node at the MIDDLE of its station, which is where a continuous
// model evaluates the pressure the side stream sees. One segment then carries
//
//   dP = (m_out |m_out| - m_arr |m_arr|) / (2 rho A^2)
//
// with m_arr the flow arriving at its upstream station and m_out the flow
// leaving its downstream station. Putting the whole drop downstream instead
// over-fed the downstream rows: -13% in Gc/Gj at row 10 of Florschuetz's
// strongest-crossflow geometry, against -2.5% centred.

namespace combaero {

struct SideStreamMomentum {
  double dP = 0.0;         // static-pressure drop along the segment [Pa]
  double d_dm_out = 0.0;   // d(dP)/d(m_out) [Pa/(kg/s)]
  double d_dm_arr = 0.0;   // d(dP)/d(m_arr) [Pa/(kg/s)]
  double d_drho = 0.0;     // d(dP)/d(rho) [Pa/(kg/m^3)]
};

SideStreamMomentum side_stream_momentum_drop(double m_arr, double m_out,
                                             double rho, double area);

// ---------------------------------------------------------------
// Merge chamber: main inlet + side streams -> one outlet (#471)
// ---------------------------------------------------------------
//
// A constant-area control volume of area A, axis = the outlet direction.
// The main stream arrives along the axis; each side stream brings axial
// momentum J_s = m_s u_jet cos(theta) (zero at 90 deg) and discharges at the
// chamber static pressure. Momentum along the axis,
//
//   P_face A + m_main u_main + S = P A + m_out u_out,   S = sum J_s,
//
// with u = m / (rho A) at the chamber density, gives the main-inlet face
// relative to the chamber (outlet) state:
//
//   P_face  - P  = (m_out|m_out| - m_main|m_main|) / (rho A^2) - S / A
//   Pt_face - Pt = (m_out|m_out| - m_main|m_main|) / (2 rho A^2) - S / A
//
// No loss coefficient: the stagnation-pressure loss of mixing follows from
// momentum. Transverse momentum is reacted by the walls and its kinetic
// energy is dissipated. With no side streams (m_out = m_main, S = 0) both
// offsets vanish, which is the single-inlet chamber.
struct ChamberMergeOffset {
  double dP_face = 0.0;    // P_face - P [Pa]
  double dPt_face = 0.0;   // Pt_face - Pt [Pa]
  double dP_dm_main = 0.0, dP_dm_out = 0.0, dP_dS = 0.0, dP_drho = 0.0;
  double dPt_dm_main = 0.0, dPt_dm_out = 0.0, dPt_dS = 0.0, dPt_drho = 0.0;
};

ChamberMergeOffset chamber_merge_face_offset(double m_main, double m_out,
                                             double side_momentum,
                                             double rho, double area);

}  // namespace combaero
