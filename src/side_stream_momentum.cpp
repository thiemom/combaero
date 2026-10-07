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

ChamberMergeOffset chamber_merge_face_offset(double m_main, double m_out,
                                             double side_momentum,
                                             double rho, double area) {
  if (!(rho > 0.0) || !(area > 0.0)) {
    throw std::invalid_argument(
        "chamber_merge_face_offset: rho and area must be positive");
  }
  const double k = 1.0 / (rho * area * area);
  const double dyn = (m_out * std::abs(m_out) - m_main * std::abs(m_main)) * k;
  ChamberMergeOffset out;
  out.dP_face = dyn - side_momentum / area;
  out.dPt_face = 0.5 * dyn - side_momentum / area;
  out.dP_dm_out = 2.0 * std::abs(m_out) * k;
  out.dP_dm_main = -2.0 * std::abs(m_main) * k;
  out.dP_dS = -1.0 / area;
  out.dP_drho = -dyn / rho;
  out.dPt_dm_out = 0.5 * out.dP_dm_out;
  out.dPt_dm_main = 0.5 * out.dP_dm_main;
  out.dPt_dS = -1.0 / area;
  out.dPt_drho = -0.5 * dyn / rho;
  return out;
}

}  // namespace combaero
