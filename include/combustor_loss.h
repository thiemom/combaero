#pragma once

// Combustor total-pressure loss correlations, as a FRACTION of the inlet
// total pressure (Pt_out = Pt_in (1 - xi)), with every partial derivative.
//
// The Python correlation classes (combaero.network.pressure_loss) call these,
// and PressureLossElement builds its Jacobian from the partials instead of
// differencing the correlation (#718). A user-supplied Python callable keeps
// the differenced path; that interface is the prototype for loss maps.
//
//   fraction family:  xi = (k theta + xi0) g(|mdot| / m_ref)
//   head family:      xi = (k theta + zeta0) q,  q = 0.5 rho v^2 / Pt_in
//                     v = mdot / (rho area),  rho = Pt_in MW / (R Tt_in)
//                     =>  q = mdot^2 R Tt_in / (2 area^2 Pt_in^2 MW)
//
// The head family is Lefebvre's pressure-loss factor, TOTAL-TOTAL at low
// Mach: dPt/Pt3 = PLF (R/2) (mdot sqrt(Tt3) / (A Pt3))^2, with the
// compressor-delivery (inlet) stagnation state (#782 audit). At a plenum
// inlet Tt = T and Pt = P; a coupled Fanno-Rayleigh form that reduces to it
// at low Mach is #830.
//
// theta = T_burned / T_unburned - 1. g(r) = 1 - (1 - r^2)^3 below r = 1 and 1
// above (C2 at r = 1): a fraction loss quoted at its design flow m_ref fades
// with the dynamic head and vanishes at zero flow (#550). m_ref <= 0 means no
// fade (a constant fraction at every flow).
//
// The composition partial is with respect to the inlet MASS fractions Y_k,
// each taken with the others held fixed (the "{node}.Y[k]" columns the
// network relays): MW = sum(Y) / sum(Y_k / W_k).

#include <vector>

namespace combaero {

struct CombustorLossResult {
  double xi = 0.0;         // fractional total-pressure loss [-]
  double dxi_dtheta = 0.0; // [-]
  double dxi_dmdot = 0.0;  // [1/(kg/s)]
  double dxi_dTt_in = 0.0; // [1/K]
  double dxi_dPt_in = 0.0; // [1/Pa]
  std::vector<double> dxi_dY; // per species [-]; empty for the fraction family
};

// xi = (k theta + xi0) g(|mdot| / m_ref).
CombustorLossResult combustor_fraction_loss(double theta, double mdot, double k,
                                            double xi0, double m_ref);

// xi = (k theta + zeta0) mdot^2 R Tt_in / (2 area^2 Pt_in^2 MW(Y_in)). Zero
// (with zero partials) for mdot <= 0, area <= 0 or Pt_in <= 0.
CombustorLossResult combustor_head_loss(double theta, double mdot, double Tt_in,
                                        double Pt_in,
                                        const std::vector<double> &Y_in,
                                        double k, double zeta0, double area);

} // namespace combaero
