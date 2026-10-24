#include "../include/combustor_loss.h"

#include "../include/thermo.h"
#include "../include/thermo_transport_data.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace combaero {

CombustorLossResult combustor_fraction_loss(double theta, double mdot, double k,
                                            double xi0, double m_ref) {
  CombustorLossResult r;
  const double base = k * theta + xi0;
  double g = 1.0;
  double dg_dm = 0.0;
  if (m_ref > 0.0) {
    const double rr = std::abs(mdot) / m_ref;
    if (rr < 1.0) {
      const double s = 1.0 - rr * rr;
      g = 1.0 - s * s * s;
      const double sign = mdot >= 0.0 ? 1.0 : -1.0;
      dg_dm = 6.0 * rr * s * s * sign / m_ref;
    }
  }
  r.xi = base * g;
  r.dxi_dtheta = k * g;
  r.dxi_dmdot = base * dg_dm;
  return r;
}

CombustorLossResult combustor_head_loss(double theta, double mdot, double T_in,
                                        double P_in,
                                        const std::vector<double> &Y_in,
                                        double k, double zeta0, double area) {
  if (Y_in.size() != molar_masses.size()) {
    throw std::invalid_argument("combustor_head_loss: Y_in size mismatch");
  }
  CombustorLossResult r;
  r.dxi_dY.assign(Y_in.size(), 0.0);
  if (mdot <= 0.0 || area <= 0.0 || P_in <= 0.0 || T_in <= 0.0) {
    return r;
  }
  // MW = S / U with S = sum(Y), U = sum(Y_k / W_k): mass_to_mole normalises,
  // so this is the MW the density sees [g/mol].
  double S = 0.0;
  double U = 0.0;
  for (std::size_t i = 0; i < Y_in.size(); ++i) {
    S += Y_in[i];
    U += Y_in[i] / molar_masses[i];
  }
  if (S <= 0.0 || U <= 0.0) {
    return r;
  }
  const double MW = S / U;
  // R in J/(mol K) over MW in g/mol: the specific gas constant is R*1000/MW.
  const double q = mdot * mdot * thermo::R_GAS * 1000.0 * T_in /
                   (2.0 * area * area * P_in * P_in * MW);
  const double c = k * theta + zeta0;
  r.xi = c * q;
  r.dxi_dtheta = k * q;
  r.dxi_dmdot = 2.0 * r.xi / mdot;
  r.dxi_dT_in = r.xi / T_in;
  r.dxi_dP_in = -2.0 * r.xi / P_in;
  // dMW/dY_k = 1/U - S/(U^2 W_k); dq/dMW = -q/MW.
  for (std::size_t i = 0; i < Y_in.size(); ++i) {
    const double dMW = 1.0 / U - S / (U * U * molar_masses[i]);
    r.dxi_dY[i] = -r.xi / MW * dMW;
  }
  return r;
}

} // namespace combaero
