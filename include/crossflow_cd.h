#pragma once

// Discharge coefficient of a wall hole under a grazing crossflow (#800).
//
// A hole whose exit (or entrance) faces a fast stream passing along the wall
// discharges far less than the same hole into still air: the stream bends
// the jet over and blocks part of the hole. Both laws below are written in
//
//   q = G_ideal / G_c
//
// G_ideal = sqrt(2 rho dP) the ideal (Cd = 1) mass velocity through the hole
// from the supply to the crossflow's static pressure, G_c the crossflow's
// mass velocity. The sources state their laws in the ACTUAL jet velocity
// (Vi = Cd Vi_ideal, Gj = Cd G_ideal), which makes the hole's own flow an
// argument of its Cd: a flow-dependent Cd whose derivative at zero flow is
// infinite. Both are resolved here so that Cd is a function of q alone; the
// hole flow then goes as dP^1.07 under a strong crossflow (a linear
// resistance), smooth through dP = 0, instead of sqrt(dP).
//
// Each law reports where it holds the hole's own (crossflow-free) Cd:
// ``capped``. With no crossflow (G_c = 0) every law is capped.

namespace combaero {

// Rogers & Hersh (1976), NASA CR-2681, Fig. 11 (p. 41): single square-edged
// holes, D 12.7 mm, L/D 1, grazing V_inf 24-85 m/s, Vi/V_inf 0.004-0.9, Cd
// referenced to the grazing STATIC pressure:
//
//   injection (into the grazing stream):  Cd = 0.97 (Vi / V_inf)^0.534
//   suction (out of the grazing stream):  Cd = 0.80 (Vi / V_inf)^0.47
//
// The paper calls injection "outflow" and suction "inflow" (from the cavity's
// side). Resolved in q (equal densities, Vi / V_inf = Gj / G_c):
// Cd = a^(1/(1-b)) q^(b/(1-b)), capped at the hole's own Cd0. The cap is
// joined C1 over +/-ROGERS_HERSH_1976_CAP_BAND of Cd0 (a min() kink there
// stalled Newton on cold CR-3630 arrays): the fit itself below 0.9 Cd0,
// exactly Cd0 above 1.1 Cd0 (every Florschuetz 1981 row sits above 2), at
// most 2.5% under the fit in between.
constexpr double ROGERS_HERSH_1976_CAP_BAND = 0.1;
constexpr double ROGERS_HERSH_1976_INJECTION_A = 0.97;
constexpr double ROGERS_HERSH_1976_INJECTION_B = 0.534;
constexpr double ROGERS_HERSH_1976_SUCTION_A = 0.80;
constexpr double ROGERS_HERSH_1976_SUCTION_B = 0.47;

// Florschuetz, Metzger & Su (1982), NASA CR-3630, Table 4.1 and Eq. 4.2:
// the special Cd tests on the B-size jet plates under crossflow. In
// x = Gc/Gj, Cd is the plate's Table 3.1 value up to a break point, then
//
//   zeta(x) = 1.06 / (x + 0.806)^0.602            (Eq. 4.2)
//
//   B(5,4,1): 0.85 to 0.63, then zeta
//   B(5,4,2): 0.85 to 0.83, then 1.08 zeta
//   B(5,4,3): 0.85 to 0.90, then 1.11 zeta
//   B(5,8,1): 0.80 to 0.54, then -0.169 x + 0.893 to 1.5, then zeta
//   B(10,4,1), B(10,8,1): 0.76 to 0.54, then -0.128 x + 0.825 to 1.8, then
//             zeta (INFERRED by the authors from the (5,x) plates, not tested)
//
// The printed pieces meet only to their printed digits (B(5,8,1)'s line
// gives 0.802 at its 0.54 break against the constant 0.80, and lands 0.0015
// below zeta at 1.5). Resolved implicitly, a jump would leave a band of q
// with no Cd at all, so the law is made continuous within those digits: the
// line keeps its printed slope but passes through zeta at its end (intercept
// 0.8945 for 0.893, 0.8256 for 0.825), and the constant joins the falling
// part C1 over +/-FLORSCHUETZ_METZGER_SU_1982_CAP_BAND of Cd0 (as Rogers &
// Hersh's cap; a min() kink stalled Newton on B(5,4,1) at mc/mj 0.99). Gj = Cd G_ideal makes x depend on Cd:
// x = 1 / (Cd q), solved here (closed form on the line, Newton on zeta) with
// the implicit derivative. The Cd is the table's own, not the caller's: Cd0
// is not an input.
constexpr double FLORSCHUETZ_METZGER_SU_1982_CAP_BAND = 0.1;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_C = 1.06;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_X0 = 0.806;
constexpr double FLORSCHUETZ_METZGER_SU_1982_ZETA_N = 0.602;

struct CrossflowCd {
  double Cd = 0.0;
  double dCd_dG_ideal = 0.0;  // [1/(kg/(m^2 s))]
  double dCd_dG_c = 0.0;      // [1/(kg/(m^2 s))]
  bool capped = false;        // the crossflow-free Cd holds
};

// Rogers & Hersh (1976). ``suction``: the hole draws FROM the grazing stream
// (a reversed jet). G_ideal <= 0 gives Cd = 0 unless G_c <= 0 (capped).
CrossflowCd rogers_hersh_1976_cd(double Cd0, double G_ideal, double G_c, bool suction);

// True when Table 4.1 covers the plate (xn/d, yn/d, z/d), to 1e-9.
bool florschuetz_metzger_su_1982_has_plate(double xn_d, double yn_d, double z_d);

// Table 4.1 for a plate it covers; throws std::invalid_argument otherwise.
CrossflowCd florschuetz_metzger_su_1982_cd(double xn_d, double yn_d, double z_d, double G_ideal,
                                           double G_c);

}  // namespace combaero
