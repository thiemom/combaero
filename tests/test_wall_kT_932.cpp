// A k(T) wall (#932): each material's analytic dk/dT, on the data and on the
// #882 envelope, against central differences of material_conductivity, and
// the total derivatives of a wall's resistance at its k(T) fixed point.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "heat_transfer.h"
#include "materials.h"

using combaero::wall_resistance_fixed_point_derivatives;
using combaero::wall_temperature_profile;
namespace mat = combaero::materials;

TEST(MaterialSlope, IsTheConductivitysDerivative) {
  for (const std::string& name : mat::list_materials()) {
    const auto& m = mat::get_material(name);
    const double lo = 0.5 * mat::envelope_T_min(m);
    const double hi = mat::envelope_T_max(m) + 150.0;
    for (bool extrapolate : {false, true}) {
      for (int i = 0; i <= 200; ++i) {
        // Odd offsets keep off the tables' printed points.
        const double T = lo + (hi - lo) * i / 200.0 + 0.3137;
        double k_ref = 0.0;
        try {
          k_ref = mat::material_conductivity(name, T, extrapolate);
        } catch (const std::runtime_error&) {
          continue;  // a continuation to k <= 0
        }
        const auto r = mat::material_conductivity_and_slope(name, T, extrapolate);
        EXPECT_EQ(r.k, k_ref) << name << " " << T;
        const double h = 1e-4;
        double left = 0.0;
        double right = 0.0;
        try {
          left = (k_ref - mat::material_conductivity(name, T - h, extrapolate)) / h;
          right = (mat::material_conductivity(name, T + h, extrapolate) - k_ref) / h;
        } catch (const std::runtime_error&) {
          continue;
        }
        if (std::abs(left - right) > 1e-4 * std::abs(left) + 1e-6) {
          continue;  // a printed point or the fade's end: C0 only
        }
        EXPECT_NEAR(r.dk_dT, 0.5 * (left + right), 1e-4 * std::abs(r.dk_dT) + 1e-6)
            << name << " " << T << " " << extrapolate;
      }
    }
  }
}

namespace {

constexpr double kRf = 1e-4;
const std::vector<double> kT = {3e-4, 2e-3, 1e-3};

std::array<double, 2> kd(std::size_t i, double T) {
  switch (i) {
    case 0:
      return {1.0 + 2e-3 * T, 2e-3};
    case 1:
      return {30.0 - 0.01 * T, -0.01};
    default:
      return {0.5 + 1e-6 * T * T, 2e-6 * T};
  }
}

std::vector<double> fixed_point(const std::array<double, 4>& p) {
  std::vector<double> r(kT.size());
  for (std::size_t i = 0; i < kT.size(); ++i) {
    r[i] = kT[i] / kd(i, 0.5 * (p[0] + p[1]))[0];
  }
  for (int it = 0; it < 200; ++it) {
    const auto prof = wall_temperature_profile(p[0], p[1], p[2], p[3], r, kRf);
    std::vector<double> next(kT.size());
    double change = 0.0;
    for (std::size_t i = 0; i < kT.size(); ++i) {
      next[i] = kT[i] / kd(i, 0.5 * (prof[i] + prof[i + 1]))[0];
      change = std::max(change, std::abs(next[i] - r[i]));
    }
    r = next;
    if (change < 1e-18) {
      break;
    }
  }
  return r;
}

double sum(const std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) {
    s += x;
  }
  return s;
}

}  // namespace

TEST(WallKT, FixedPointDerivativesMatchCentralDifferences) {
  const std::array<double, 4> p = {1600.0, 500.0, 900.0, 300.0};
  const auto r = fixed_point(p);
  const auto prof = wall_temperature_profile(p[0], p[1], p[2], p[3], r, kRf);
  std::vector<double> dr(kT.size());
  for (std::size_t i = 0; i < kT.size(); ++i) {
    const auto k = kd(i, 0.5 * (prof[i] + prof[i + 1]));
    dr[i] = -kT[i] * k[1] / (k[0] * k[0]);
  }
  const auto got = wall_resistance_fixed_point_derivatives(p[0], p[1], p[2], p[3], r, kRf, dr);
  for (std::size_t j = 0; j < 4; ++j) {
    const double e = 1e-6 * p[j];
    auto hi = p;
    auto lo = p;
    hi[j] += e;
    lo[j] -= e;
    const double fd = (sum(fixed_point(hi)) - sum(fixed_point(lo))) / (2.0 * e);
    EXPECT_NEAR(got[j], fd, 1e-6 * std::abs(fd)) << j;
  }
  const auto none = wall_resistance_fixed_point_derivatives(
      p[0], p[1], p[2], p[3], r, kRf, std::vector<double>(kT.size(), 0.0));
  for (double v : none) {
    EXPECT_EQ(v, 0.0);
  }
}
