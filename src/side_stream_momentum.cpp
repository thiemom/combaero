#include "side_stream_momentum.h"

#include <cmath>
#include <stdexcept>

namespace combaero {

SideStreamMomentum side_stream_momentum_drop(double m_arr, double m_out,
                                             double rho, double area) {
  if (!(rho > 0.0) || !(area > 0.0)) {
    throw std::invalid_argument(
        "side_stream_momentum_drop: rho and area must be positive");
  }
  const double k = 0.5 / (rho * area * area);
  SideStreamMomentum out;
  out.dP = (m_out * std::abs(m_out) - m_arr * std::abs(m_arr)) * k;
  out.d_dm_out = 2.0 * std::abs(m_out) * k;
  out.d_dm_arr = -2.0 * std::abs(m_arr) * k;
  out.d_drho = -out.dP / rho;
  return out;
}

}  // namespace combaero
