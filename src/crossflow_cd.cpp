#include "crossflow_cd.h"

#include <cmath>
#include <stdexcept>

namespace combaero {

namespace {

// d(Cd)/d(q) carried onto (G_ideal, G_c) through q = G_ideal / G_c.
CrossflowCd from_q(double Cd, double dCd_dq, double q, double G_c) {
  CrossflowCd out;
  out.Cd = Cd;
  out.dCd_dG_ideal = dCd_dq / G_c;
  out.dCd_dG_c = -dCd_dq * q / G_c;
  return out;
}

CrossflowCd capped(double Cd0) {
  CrossflowCd out;
  out.Cd = Cd0;
  out.capped = true;
  return out;
}

struct Table41Plate {
  double xn_d, yn_d, z_d;
  double cd0;       // Table 3.1, up to the break point
  double x_break;   // the printed break point (documentation; the law joins by min)
  double mult;      // multiplier on zeta
  double lin_a;     // linear piece Cd = lin_a x + lin_b (0: none)
  double lin_b;     // printed intercept (documentation; see below)
  double x_lin_end; // where the linear piece hands over to zeta
};

constexpr Table41Plate TABLE_4_1[] = {
    {5.0, 4.0, 1.0, 0.85, 0.63, 1.00, 0.0, 0.0, 0.0},
    {5.0, 4.0, 2.0, 0.85, 0.83, 1.08, 0.0, 0.0, 0.0},
    {5.0, 4.0, 3.0, 0.85, 0.90, 1.11, 0.0, 0.0, 0.0},
    {5.0, 8.0, 1.0, 0.80, 0.54, 1.00, -0.169, 0.893, 1.5},
    {10.0, 4.0, 1.0, 0.76, 0.54, 1.00, -0.128, 0.825, 1.8},
    {10.0, 8.0, 1.0, 0.76, 0.54, 1.00, -0.128, 0.825, 1.8},
};

const Table41Plate* find_plate(double xn_d, double yn_d, double z_d) {
  for (const auto& p : TABLE_4_1) {
    if (std::abs(p.xn_d - xn_d) < 1e-9 && std::abs(p.yn_d - yn_d) < 1e-9 &&
        std::abs(p.z_d - z_d) < 1e-9) {
      return &p;
    }
  }
  return nullptr;
}

}  // namespace

CrossflowCd rogers_hersh_1976_cd(double Cd0, double G_ideal, double G_c, bool suction) {
  if (!(Cd0 > 0.0)) {
    throw std::invalid_argument("rogers_hersh_1976_cd: Cd0 must be positive");
  }
  if (!(G_c > 0.0)) {
    return capped(Cd0);
  }
  if (!(G_ideal > 0.0)) {
    return CrossflowCd{};
  }
  const double a = suction ? ROGERS_HERSH_1976_SUCTION_A : ROGERS_HERSH_1976_INJECTION_A;
  const double b = suction ? ROGERS_HERSH_1976_SUCTION_B : ROGERS_HERSH_1976_INJECTION_B;
  // Cd = a (Cd q)^b  =>  Cd = a^(1/(1-b)) q^(b/(1-b)).
  const double p = 1.0 / (1.0 - b);
  const double e = b * p;
  const double q = G_ideal / G_c;
  const double cd = std::pow(a, p) * std::pow(q, e);
  // C1 join to the cap: u = cd / Cd0 is kept below u0 = 1 - w, held at 1
  // above u1 = 1 + w, and the quadratic u - (u - u0)^2 / (4 w) in between.
  const double w = ROGERS_HERSH_1976_CAP_BAND;
  const double u = cd / Cd0;
  if (u >= 1.0 + w) {
    return capped(Cd0);
  }
  const double dcd_dq = e * cd / q;
  if (u <= 1.0 - w) {
    return from_q(cd, dcd_dq, q, G_c);
  }
  const double s = u - (1.0 - w);
  const double f = u - s * s / (4.0 * w);
  const double df_du = 1.0 - s / (2.0 * w);
  return from_q(Cd0 * f, df_du * dcd_dq, q, G_c);
}

bool florschuetz_metzger_su_1982_has_plate(double xn_d, double yn_d, double z_d) {
  return find_plate(xn_d, yn_d, z_d) != nullptr;
}

CrossflowCd florschuetz_metzger_su_1982_cd(double xn_d, double yn_d, double z_d, double G_ideal,
                                           double G_c) {
  const Table41Plate* plate = find_plate(xn_d, yn_d, z_d);
  if (plate == nullptr) {
    throw std::invalid_argument(
        "florschuetz_metzger_su_1982_cd: CR-3630 Table 4.1 covers the B(5,4,1), "
        "B(5,4,2), B(5,4,3), B(5,8,1), B(10,4,1) and B(10,8,1) plates only");
  }
  if (!(G_c > 0.0)) {
    return capped(plate->cd0);
  }
  if (!(G_ideal > 0.0)) {
    return CrossflowCd{};
  }
  const double q = G_ideal / G_c;
  const double m = plate->mult * FLORSCHUETZ_METZGER_SU_1982_ZETA_C;
  const double x0 = FLORSCHUETZ_METZGER_SU_1982_ZETA_X0;
  const double n = FLORSCHUETZ_METZGER_SU_1982_ZETA_N;
  const bool has_line = plate->lin_a != 0.0;
  const double A = plate->lin_a;
  // The line through zeta at its end, with the printed slope (its intercept
  // moves 0.0015 / 0.0006 from the printed 0.893 / 0.825), so the pieces
  // meet. Within the printed digits.
  const double end = plate->x_lin_end;
  const double B = has_line ? m * std::pow(end + x0, -n) - A * end : 0.0;
  // The falling part resolved in q, unbounded by cd0: the line's root,
  // or zeta's (Newton; dF/dCd = 1 - n x/(x + x0) > 0, a unique root).
  double cd_f = 0.0;
  double dcdf_dq = 0.0;
  bool on_line = false;
  if (has_line) {
    // Cd = A / (Cd q) + B  =>  q Cd^2 - B q Cd - A = 0, the larger root.
    const double disc = B * B * q * q + 4.0 * A * q;
    if (disc >= 0.0) {
      const double cd = (B * q + std::sqrt(disc)) / (2.0 * q);
      if (cd > 0.0 && 1.0 / (cd * q) <= end) {
        const double F_cd = 1.0 + A / (cd * cd * q);
        const double F_q = A / (cd * q * q);
        cd_f = cd;
        dcdf_dq = -F_q / F_cd;
        on_line = true;
      }
    }
  }
  if (!on_line) {
    double cd = plate->cd0;
    const double cd_max = 4.0 * plate->cd0;
    for (int it = 0; it < 100; ++it) {
      const double x = 1.0 / (cd * q);
      const double z = m * std::pow(x + x0, -n);
      const double dz_dx = -n * z / (x + x0);
      const double F_cd = 1.0 + dz_dx / (cd * cd * q);  // dx/dCd = -x/Cd
      double next = cd - (cd - z) / F_cd;
      if (!(next > 0.0)) {
        next = 0.5 * cd;
      }
      if (next > cd_max) {
        next = cd_max;
      }
      const bool done = std::abs(next - cd) <= 1e-15 * plate->cd0;
      cd = next;
      if (done) {
        break;
      }
    }
    const double x = 1.0 / (cd * q);
    const double dz_dx = -n * m * std::pow(x + x0, -n) / (x + x0);
    const double F_cd = 1.0 + dz_dx / (cd * cd * q);
    const double F_q = dz_dx / (cd * q * q);  // dx/dq = -x/q
    cd_f = cd;
    dcdf_dq = -F_q / F_cd;
  }
  // The constant joins C1 as for Rogers & Hersh: the falling part below
  // (1 - w) cd0, cd0 above (1 + w) cd0, the quadratic in between.
  const double w = FLORSCHUETZ_METZGER_SU_1982_CAP_BAND;
  const double u = cd_f / plate->cd0;
  if (u >= 1.0 + w) {
    return capped(plate->cd0);
  }
  if (u <= 1.0 - w) {
    return from_q(cd_f, dcdf_dq, q, G_c);
  }
  const double sb = u - (1.0 - w);
  const double f = u - sb * sb / (4.0 * w);
  return from_q(plate->cd0 * f, (1.0 - sb / (2.0 * w)) * dcdf_dq, q, G_c);
}

}  // namespace combaero
