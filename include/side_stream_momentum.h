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

}  // namespace combaero
