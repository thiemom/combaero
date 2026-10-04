# CombAero C++ API Reference

This document provides the technical reference for the CombAero C++ library. For Python usage, see [API_PYTHON.md](API_PYTHON.md).

## Table of Contents
- [Thermodynamics (thermo.h)](#thermodynamics-thermoh)
- [Transport Properties (transport.h)](#transport-properties-transporth)
- [Combustion (combustion.h)](#combustion-combustionh)
- [Chemical Equilibrium (equilibrium.h)](#chemical-equilibrium-equilibriumh)
- [Stagnation / Static Conversions (stagnation.h)](#stagnation--static-conversions-stagnationh)
- [Compressible Flow (compressible.h)](#compressible-flow-compressibleh)
- [Incompressible Flow (incompressible.h)](#incompressible-flow-incompressibleh)
- [Friction Factor Correlations (friction.h)](#friction-factor-correlations-frictionh)
- [Heat Transfer Correlations (heat_transfer.h)](#heat-transfer-correlations-heat_transferh)
  - [Nusselt Numbers](#nusselt-numbers)
  - [Overall Heat Transfer](#overall-heat-transfer)
  - [Channel Flow Models](#channel-flow-models)
  - [Wall Coupling](#wall-coupling)
  - [Data Structures](#data-structures)
- [Geometry Utilities (geometry.h)](#geometry-utilities-geometryh)
- [Orifice Flow (orifice.h)](#orifice-flow-orificeh)
- [Acoustics (acoustics.h)](#acoustics-acousticsh)
- [Humid Air (humidair.h)](#humid-air-humidairh)
- [Network Solver Interface (solver_interface.h)](#network-solver-interface-solver_interfaceh)
  - [Tee Junction Components](#tee-junction-components-tee_junctionh--solver_interfaceh)

---

## Thermodynamics (thermo.h)

> [!NOTE]
> **Base Units**: Temperature [K], Pressure [Pa], Mass [kg], Energy [J], Power [W].
> Compositions are `std::vector<double>` of mole fractions (sum to 1.0) unless explicitly marked as `_mass` (mass fractions).

### Composition Utilities (composition.h)
```cpp
// Mass/Mole fraction conversions
double mwmix(const std::vector<double>& X);
std::vector<double> mole_to_mass(const std::vector<double>& X);
std::vector<double> mass_to_mole(const std::vector<double>& Y);
std::vector<double> normalize_fractions(const std::vector<double>& X);
std::vector<double> convert_to_dry_fractions(const std::vector<double>& X);
```

### Species and Mixture Metadata (thermo.h)
```cpp
// Get species metadata
std::string species_name(std::size_t index);
std::size_t species_index_from_name(const std::string& name);
double species_molar_mass(std::size_t index);
double species_molar_mass_from_name(const std::string& name);
std::size_t num_species();
```

### Thermodynamic Properties (Molar Basis)
```cpp
double cp(double T, const std::vector<double>& X);
double cv(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);
double s(double T, const std::vector<double>& X, double P, double P_ref = 101325.0);
```

### Thermodynamic Properties (Mass Basis)
```cpp
double cp_mass(double T, const std::vector<double>& X);
double cv_mass(double T, const std::vector<double>& X);
double h_mass(double T, const std::vector<double>& X);
double u_mass(double T, const std::vector<double>& X);
double s_mass(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
double density(double T, double P, const std::vector<double>& X);
double speed_of_sound(double T, const std::vector<double>& X);
double isentropic_expansion_coefficient(double T, const std::vector<double>& X);
```

### Inverse Solvers
```cpp
double calc_T_from_h(double h_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_cp(double cp_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u(double u_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_h_mass(double h_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_s_mass(double s_mass_target, double P, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
double calc_T_from_u_mass(double u_mass_target, const std::vector<double>& X, double T_guess = 300.0, double tol = 1.0e-6, std::size_t max_iter = 50);
```

### Complete State Functions
```cpp
struct AirProperties { /* air properties at once */ };
struct ThermoState { /* thermodynamic properties */ };
struct CompleteState { /* thermodynamic + transport properties */ };

AirProperties air_properties(double T, double P, double humidity = 0.0);
ThermoState thermo_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
CompleteState complete_state(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
```

### Derivatives
```cpp
double dh_dT(double T, const std::vector<double>& X);
double ds_dT(double T, const std::vector<double>& X);
double dcp_dT(double T, const std::vector<double>& X);
double dg_over_RT_dT(double T, const std::vector<double>& X);
```

### State-Based Overloads
```cpp
double mwmix(const State& s);
double cp(const State& s);
double h(const State& s);
double s(const State& s);
double cv(const State& s);
double u(const State& s);
double density(const State& s);
double specific_gas_constant(const State& s);
double isentropic_expansion_coefficient(const State& s);
double speed_of_sound(const State& s);
```

---

## Transport Properties (transport.h)

> [!NOTE]
> All transport properties are evaluated at the local (T, P) state.
> **Units**: Viscosity [Pa·s], Conductivity [W/(m·K)], Diffusivity [m²/s].

```cpp
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);
double reynolds_from_state(double rho, double v, double L, double mu);
```

### Transport Model

Viscosity uses Chapman-Enskog kinetic theory with the Monchick-Mason
Omega*(2,2) collision integral (37 x 8 bilinear table in T* and delta*).
Polar species (H2O, NH3, CO) use the full 2D table; non-polar species use the
delta*=0 column. Mixture viscosity uses the Wilke mixing rule.

Thermal conductivity uses the Mason-Monchick modified Eucken formula with
Parker Z_rot temperature correction, matching the Cantera GasTransport model.

Species transport parameters are stored in `transport_props` (thermo_transport_data.h):

```cpp
struct Transport_Props {
    MolecularGeometry geometry;  // Atom, Linear, or Nonlinear
    double well_depth;           // epsilon/k_B [K]
    double diameter;             // sigma [Ang]
    double polarizability;       // alpha [Ang^3]; 0.0 for non-polar
    double dipole_moment;        // mu [Debye]; 0.0 for non-polar
    double z_rot;                // rotational relaxation collision number [-]
};
```

---

## Combustion (combustion.h)

> [!NOTE]
> Combustion functions assume complete combustion to CO₂ and H₂O unless equilibrium is invoked.
> `smooth=True` enables a sigmoid-based transition at Phi=0 for better Jacobian stability.

### Stoichiometry & Equivalence Ratio
```cpp
double oxygen_required_per_mol_fuel(std::size_t fuel_index);
double oxygen_required_per_kg_fuel(std::size_t fuel_index);

double equivalence_ratio_mole(const std::vector<double>& X_mix,
                              const std::vector<double>& X_fuel,
                              const std::vector<double>& X_ox);
double equivalence_ratio_mass(const std::vector<double>& Y_mix,
                              const std::vector<double>& Y_fuel,
                              const std::vector<double>& Y_ox);

double fuel_lhv_molar(const std::vector<double>& X_fuel, double T_ref = 298.15);
double fuel_lhv_mass(const std::vector<double>& X_fuel, double T_ref = 298.15);
```

### Complete Combustion Solvers
```cpp
// All solvers accept an optional 'smooth' parameter for Jacobian smoothness
State complete_combustion(const State& in, bool smooth = false);
State complete_combustion_isothermal(const State& in, bool smooth = false);
std::vector<double> complete_combustion_to_CO2_H2O(const std::vector<double>& X);

struct PressureLossContext {
    double m_dot;
    double P;
    double T;
    const std::vector<double>& Y;
    const std::vector<double>& Y_products;
    double theta;
};

using PressureLossCorrelation = std::function<std::tuple<double, double>(const PressureLossContext &)>;
```

---

## Chemical Equilibrium (equilibrium.h)

```cpp
// Combustion equilibrium (T, P constant)
State combustion_equilibrium(const State& in, bool smooth = false);

// Shift and reforming reactions
State wgs_equilibrium(const State& in);
State wgs_equilibrium_adiabatic(const State& in);
State smr_wgs_equilibrium(const State& in);
State reforming_equilibrium(const State& in);
```

---

## Stagnation / Static Conversions (stagnation.h)

```cpp
double T0_from_static(double T, double M, const std::vector<double>& X);
double P0_from_static(double P, double T, double M, const std::vector<double>& X);
double T_from_stagnation(double T0, double M, const std::vector<double>& X);
double P_from_stagnation(double P0, double T0, double M, const std::vector<double>& X);

double T_adiabatic_wall(double T_static, double v, double T, double P,
                        const std::vector<double>& X, bool turbulent = true);
double recovery_factor(double Pr, bool turbulent = true);

double bulk_velocity(double m_dot, double rho, double area);
double kinetic_energy(double v);

std::tuple<double, double> stagnation_from_static(
    double T, double P, double v,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);
```

---

## Compressible Flow (compressible.h)

> [!IMPORTANT]
> All compressible solvers assume **ideal-gas** behavior with variable properties (gamma and Cp are recalculated at each state).
> **Choking**: Functions return a `choked` flag. If true, the mass flow is limited by the sonic condition at the throat.

### Results Structs
```cpp
struct CompressibleFlowSolution {
    State stagnation;
    State outlet;
    double v;
    double M;
    double mdot;
    bool choked;
};

struct FannoSolution {
    State inlet, outlet;
    double mdot, h0, L, D, f_avg;
    bool choked;
    double L_choke;
    std::vector<FannoStation> profile;
};
```

### Solvers
```cpp
CompressibleFlowSolution nozzle_flow(double T0, double P0, double P_back,
                                     double A_eff, const std::vector<double>& X,
                                     double tol = 1e-8, std::size_t max_iter = 50);

FannoSolution fanno_channel(double T_in, double P_in, double u_in, double L, double D,
                          double f, const std::vector<double>& X,
                          std::size_t n_steps = 100, bool store_profile = false);

FannoSolution fanno_channel_rough(double T_in, double P_in, double u_in, double L, double D,
                                double roughness, const std::vector<double>& X,
                                const std::string& correlation = "haaland",
                                std::size_t n_steps = 100, bool store_profile = false);
```

### Quasi-1D Nozzle Flow
```cpp
// Area function type: A(x) returning area [m²] at position x [m]
using AreaFunction = std::function<double(double)>;

struct NozzleStation {
    double x, A, P, T, rho, u, M, h;
};

struct NozzleSolution {
    State inlet, outlet;
    double mdot, h0, T0, P0;
    bool choked;
    double x_throat, A_throat;
    std::vector<NozzleStation> profile;
};

NozzleSolution nozzle_quasi1d(double T0, double P0, double P_exit,
                             const AreaFunction& area_func,
                             double x_start, double x_end,
                             const std::vector<double>& X,
                             std::size_t n_stations = 100);

NozzleSolution nozzle_quasi1d(double T0, double P0, double P_exit,
                             const std::vector<std::pair<double, double>>& area_profile,
                             const std::vector<double>& X,
                             std::size_t n_stations = 100);

NozzleSolution nozzle_cd(double T0, double P0, double P_exit,
                         double A_inlet, double A_throat, double A_exit,
                         double x_throat, double x_exit,
                         const std::vector<double>& X,
                         std::size_t n_stations = 100);
```

### Inverse Solvers
```cpp
double solve_A_eff_from_mdot(double T0, double P0, double P_back, double mdot_target,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double solve_P_back_from_mdot(double T0, double P0, double A_eff, double mdot_target,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double solve_P0_from_mdot(double T0, double P_back, double A_eff, double mdot_target,
                         const std::vector<double>& X,
                         double tol = 1e-8, std::size_t max_iter = 50);
```

### Utility Functions
```cpp
double critical_pressure_ratio(double T0, double P0, const std::vector<double>& X,
                               double tol = 1e-8, std::size_t max_iter = 50);

double mach_from_pressure_ratio(double T0, double P0, double P,
                               const std::vector<double>& X,
                               double tol = 1e-8, std::size_t max_iter = 50);

double mass_flux_isentropic(double T0, double P0, double P,
                             const std::vector<double>& X,
                             double tol = 1e-8, std::size_t max_iter = 50);

double fanno_max_length(double T_in, double P_in, double u_in,
                       double D, double f, const std::vector<double>& X,
                       double tol = 1e-6, std::size_t max_iter = 100);
```

### Rocket Nozzle Thrust
```cpp
struct ThrustResult {
    double thrust;
    double specific_impulse;
    double thrust_coefficient;
    double mdot;
    double u_exit;
    double P_exit;
};

ThrustResult nozzle_thrust(const NozzleSolution& sol, double P_amb);

ThrustResult nozzle_thrust(double T0, double P0, double P_design, double P_amb,
                           double A_inlet, double A_throat, double A_exit,
                           double x_throat, double x_exit,
                           const std::vector<double>& X,
                           std::size_t n_stations = 100);
```

---

## Incompressible Flow (incompressible.h)

```cpp
double bernoulli_P2(double P1, double v1, double v2, double rho, double dz = 0.0);
double orifice_mdot(double P1, double P2, double A, double Cd, double rho);
double channel_dP(double v, double L, double D, double f, double rho);
```

---

## Friction Factor Correlations (friction.h)

```cpp
double friction_haaland(double Re, double e_D);
double friction_serghides(double Re, double e_D);
double friction_colebrook(double Re, double e_D, double tol = 1e-10, int max_iter = 20);
double friction_petukhov(double Re);
// Petukhov held below Re 3000 by a C1 soft-max over Re 2500-3500 (#446);
// used by every smooth-pipe Gnielinski / Petukhov heat-transfer path.
double friction_petukhov_clamped(double Re);
double friction_petukhov_clamped_dRe(double Re);
```

---

## Heat Transfer Correlations (heat_transfer.h)

> [!NOTE]
> **Units**: HTC [W/(m²·K)], q [W/m²].
> **Jacobians**: High-performance solver components return `dh/dmdot` and other derivatives essential for coupled fluid-thermal networks.

### Nusselt Numbers
```cpp
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);
double nusselt_gnielinski(double Re, double Pr);
double nusselt_gnielinski(double Re, double Pr, double f);
// Value and EXACT Re-derivative (no warnings); df_dRe = slope of the f passed.
NuAndDerivative nusselt_gnielinski_with_derivative(double Re, double Pr,
                                                   double f, double df_dRe = 0.0);
NuAndDerivative nusselt_gnielinski_smooth_with_derivative(double Re, double Pr);

// Channel regime transitions (#448): every switch is a C1 smoothstep, each
// regime EXACT outside its band. channel_smooth, nusselt_circular_channel and
// htc_circular_channel (gnielinski branch) use these, so neither the value nor
// the Jacobian steps at Re 2300 or 4000.
//   laminar -> turbulent  Re 2300-3000   f: 64/Re -> turbulent; Nu: 3.66/4.36 -> Gnielinski
//   smooth  -> rough      Re 3000-4000   f: clamped Petukhov -> Colebrook (e_D > 0)
FrictionAndDerivative friction_channel_and_derivative(double Re, double e_D);   // Darcy, exact df/dRe
FrictionAndDerivative friction_turbulent_and_derivative(double Re, double e_D);
NuAndDerivative nusselt_channel_gnielinski_and_derivative(
    double Re, double Pr, double f_turb, double df_turb_dRe, double Nu_laminar);
double friction_colebrook_dRe(double Re, double e_D);  // implicit, exact (friction.h)
double nusselt_sieder_tate(double Re, double Pr, double mu_ratio);
```

### Overall Heat Transfer
```cpp
double overall_htc(const std::vector<double>& h_values, const std::vector<double>& t_over_k);
double overall_htc_wall(double h_inner, double h_outer, const std::vector<double>& t_over_k_layers);
```

### Channel Flow Models
```cpp
// Smooth channel
ChannelResult channel_smooth(double T, double P, const std::vector<double>& X,
                              double velocity, double diameter, double length,
                              double T_hot = std::numeric_limits<double>::quiet_NaN(),
                              const std::string& correlation = "gnielinski",
                              bool heating = true, double Nu_multiplier = 1.0,
                              double f_multiplier = 1.0);

// Parametrised rib correlations -- rib_correlation.h
//
//   R = C_R * (e/D / nD)^a * (p/e / nP)^b * (W/H / nW)^c * (alpha/90)^d
//   G = C_G * (same geometry terms) * (e+)^n
//
// R carries no e+ term by construction, so f is independent of Reynolds
// number and df/d(mdot) is identically zero for this family.
//
// han_park_1988_angled uses two shapes the plain power law above cannot
// express: R is a quadratic in alpha (Eq. 4.17), and G's alpha/p_e
// exponents switch on channel shape, square vs rectangular (Eq. 4.18).
// Both are genuine discontinuities the source states -- up to 62.5% in R
// at alpha=90, W/H=4; up to 27% in G at W/H=1 -- and neither is smoothed,
// since no smooth transition is stated. See
// validation/cooling/extractions/han_ribbed.md, decision D5, and
// RibCorrelationSet::RAlphaShape / GShapeModel in rib_correlation.h.
//
// han_1989_narrow_channel (Eq. 4.19, W/H < 1) adds two more shapes. R has
// NO PRINTED EQUATION in this source -- only a drawn figure line -- so its
// quadratic-in-alpha coefficients (QuadraticAlphaTwoBand, one per W/H
// sub-band) are this project's OWN FIT to that line, not an extraction;
// provenance is Fitted. G (NarrowChannelAlphaSwitch) IS text-extracted: a
// hard ~20% switch in its leading constant at alpha == 90 deg, plus a
// further W/H correction to both the constant and the e+ exponent itself
// below W/H = 1/2. See han_ribbed.md item 40 and decision D6.
//
// Rib SHAPE is part of a set (#434): continuous-angled, V, crossed and
// Lambda ribs at the same alpha are different geometries, and a set refuses
// a series of another shape rather than aliasing onto it. Unspecified (the
// default for a user-built set) binds nothing.
enum class RibShape { Unspecified, Transverse, Parallel, Crossed, V, Lambda };
// set.shape; set.C_Gbar / set.Gbar_eplus_exponent carry a source's OWN
// printed four-wall G_bar (0 = not printed -> RibResult::has_G_bar false,
// callers fall back to Han's G_bar = 1.2 G).
//
// Named sets, selected explicitly -- there is no auto-switching between
// them, since they cover disjoint regimes rather than one being an update
// of the other:
RibCorrelationSet han_1988_orthogonal();       // 90 deg, Re 10e3-60e3
RibCorrelationSet rallabandi_2009_high_re();   // 45 deg sharp ribs, Re 30e3-400e3
RibCorrelationSet han_park_1988_angled();      // 30-90 deg, W/H 1-4, Re 10e3-60e3
RibCorrelationSet han_1989_narrow_channel();   // 30-90 deg, W/H 1/4-1, Re 10e3-60e3
// Han, Zhang & Lee (1991) Table 2: Transverse at 90; Parallel, Crossed, V,
// Lambda at 45 or 60. One rig (square, e/D 0.0625, P/e 10); anything else
// throws. Prints its own G_bar per configuration.
RibCorrelationSet han_zhang_lee_1991(RibShape shape, double alpha_deg);
void validate_rib_set(const RibCorrelationSet& set);   // throws on a bad set

// Accuracy carries WHERE IT CAME FROM, because only an author's own claim
// may be used as the band a model is judged against. A figure this project
// measured is that model's own error, and scoring a model inside its own
// error answers nothing (#389).
enum class AccuracyProvenance { Unstated, Stated, Measured };
struct StatedAccuracy {
  double value = NaN;                 // NaN unless provenance says otherwise
  AccuracyProvenance provenance = AccuracyProvenance::Unstated;
  bool usable_as_band() const;        // true only for Stated
  static StatedAccuracy stated(double v);
  static StatedAccuracy measured(double v);
  static StatedAccuracy unstated();
};
// set.accuracy_R / set.accuracy_G are StatedAccuracy. Of the four shipped
// sets only han_1988_orthogonal's are Stated (Han's "95% within 6%/8%");
// han_park_1988_angled's and han_1989_narrow_channel's were measured
// through evaluate_rib itself, and rallabandi_2009_high_re's accuracy_R is
// Unstated -- previously 0.0, which read as perfect agreement.
RibResult evaluate_rib(const RibCorrelationSet& set,
                       const RibGeometry& geom, double Re);
// evaluate_rib never throws: Re may be negative or zero and the guards are
// smooth through both, because the solver probes states that are not physical.

// Ratio-form rib sets -- rib_ratio_correlation.h (#444). Nu/Nu0 and f/f0 as
// multipliers on the baseline the SOURCE fitted to, beside the R/G sets, for
// sources that publish ratios (and so that their data is not forced through
// the law-of-the-wall similarity R/G assumes).
//
//   in range:    Nu = r(Re, geom) * Nu0_source        (reproduces the paper)
//   below floor: C1 blend in ln Re over [Re_floor/2, Re_floor] to
//                k * Nu0_ext, k matching the value at Re_floor;
//                Nu0_ext = smooth-pipe Gnielinski by default (laminar 3.66
//                limit at Re -> 0), the source baseline, or a user power law.
//   friction:    ratio held at its floor value below the range.
//
// Never divides by Nu0; Re enters as sqrt(Re^2 + 1) so Nu, f are even and
// their analytic derivatives odd and finite through Re = 0.
struct RatioBaseline { double coeff, re_exponent, pr_exponent; };  // c Re^m Pr^n
enum class RatioBelowFloor { Gnielinski, SourceBaseline, User };
struct RibRatioOptions { RatioBelowFloor below_floor; RatioBaseline user_Nu0; };
RibRatioResult evaluate_rib_ratio(const RibRatioSet& set, const RibGeometry& geom,
                                  double Re, double Pr,
                                  const RibRatioOptions& options = {});
// -> Nu, dNu_dRe, f, df_dRe, ratio_Nu, ratio_f, below_floor, extrapolated
void validate_rib_ratio_set(const RibRatioSet& set);
void validate_rib_ratio_options(const RibRatioOptions& options);
// set.Re_floor_f: friction's own floor (0 = Re_floor), for sources that
// measured f over a different Re range than Nu.
RibRatioSet taslim_spring_1987(double aspect_ratio_taslim, double e_D);
// Taslim & Spring (1987), two ribbed walls, transverse, p/e 10. Seven tested
// configurations: AR 0.5 (e/D 0.125, 0.250), 1.0 (0.083, 0.167), 3.5 (0.053,
// 0.107, 0.161); any other throws. Taslim's AR is height/width: W/H = 1/AR.
// Nu/Nu_DB = C (Re/1e4)^-0.2 (fig. 9), f Fanning passage-average constant
// (fig. 11). Provenance Fitted; accuracy Unstated.

// Enhanced surfaces removed in 0.7.0 for unprovenanced correlations (issue
// #339). Ribbed is re-added on han_1988_orthogonal / rallabandi_2009_high_re
// above, wired through RibbedModel (components.py) into the GUI and network
// solver. Impingement is re-added below, wired through ImpingementModel /
// SingleJetImpingementModel (components.py) into the GUI and network solver
// the same way. Dimpled and pin-fin remain removed, tracked individually
// (#335-336).
```

### Jet Impingement Correlations
```cpp
// impingement_correlation.h -- two independent regimes, matching
// han_impingement.md's own split, wired into ImpingementModel /
// SingleJetImpingementModel (components.py). Converting an array's
// total/mean mass flow into a per-row local Re_j (Florschuetz's own Eq. 7)
// is deliberately NOT implemented here -- ImpingementModel recovers each
// row's own mass flow from its own ConvectiveSurface.area instead, leaving
// row-to-row flow distribution to whoever assembles the chain of elements.

// Single jet (Goldstein, Behbahani and Heppelmann, 1986):
//   Nu_bar = Re^0.76 * (A - |L/D - 7.75|) / (B + C*(R/D)^n)
// n switches with the boundary condition the correlation was fitted under.
enum class ImpingementThermalBC { ConstantHeatFlux, ConstantWallTemperature };
SingleJetImpingementSet goldstein_1986_single_jet();
void validate_single_jet_set(const SingleJetImpingementSet& set);
double single_jet_impingement_nu(const SingleJetImpingementSet& set,
                                 ImpingementThermalBC bc, double Re,
                                 double L_D, double R_D);
// Same formula, plus an analytic d(Nu)/d(Re) -- for a network element's
// wall-coupling Jacobian, per the project's (f, J) rule.
SingleJetImpingementResult single_jet_impingement(
    const SingleJetImpingementSet& set, ImpingementThermalBC bc, double Re,
    double L_D, double R_D);

// Jet array with crossflow (Florschuetz, Truman and Metzger, 1981):
//   Nu = A * Re_j^m * {1 - B*[(z/d)(Gc/Gj)]^n} * Pr^(1/3)
// A, m, B, n are each C*(xn/d)^nx*(yn/d)^ny*(z/d)^nz (Table 4.1). Re_j and
// Gc/Gj are the ROW's own local values, not an array mean.
// JetArrayImpingementResult::dNu_dRe_j is analytic, holding Gc/Gj and
// geometry fixed -- Gc/Gj never depends on Re_j or mdot, only on geometry.
enum class JetHolePattern { Inline, Staggered };
JetArrayCorrelationSet florschuetz_1981_inline();
JetArrayCorrelationSet florschuetz_1981_staggered();
void validate_jet_array_set(const JetArrayCorrelationSet& set);
JetArrayImpingementResult jet_array_impingement_nu(
    const JetArrayCorrelationSet& set, double Re_j, double Gc_Gj, double Pr,
    double xn_d, double yn_d, double z_d);

// Crossflow-to-jet mass flux ratio, Florschuetz's own closed form (Eq. 8),
// depending on (yn/d)(z/d) only -- not xn/d, per the source. C_D is the
// jet-plate discharge coefficient (FLORSCHUETZ_1981_DEFAULT_CD = 0.79
// absent a measured value; combaero has no correlation of its own for a
// jet-plate array yet -- see #375, distinct from orifice.h's pipe-metering
// Cd family, which does not apply to this geometry).
constexpr double FLORSCHUETZ_1981_DEFAULT_CD = 0.79;
double crossflow_to_jet_ratio_at_x(double yn_d, double z_d, double C_D,
                                   double x_over_xn);
double crossflow_to_jet_ratio_at_row(double yn_d, double z_d, double C_D,
                                     int row);  // row 1 is exactly 0
```

### Wall Coupling
```cpp
WallCouplingResult wall_coupling_and_jacobian(
    double h_a, double T_aw_a, double h_b, double T_aw_b,
    double t_over_k,   // wall_thickness / wall_conductivity [m²·K/W]
    double A = 1.0     // contact area [m²]
);
```

### Data Structures
```cpp
struct ChannelResult {
    // Primary outputs
    double h;              // Heat transfer coefficient [W/(m²·K)]
    double Nu;             // Nusselt number [-]
    double Re;             // Reynolds number [-]
    double Pr;             // Prandtl number [-]
    double f;              // Friction factor [-]
    double dP;             // Pressure drop [Pa]
    double M;              // Mach number [-]
    double T_aw;           // Adiabatic wall temperature [K]
    double q;              // Heat flux [W/m²] (nan if T_hot not supplied)

    // Jacobians
    double dh_dmdot;       // dh/dmdot [W/(m²·K·kg/s)]
    double dh_dT;          // dh/dT [W/(m²·K²)]
    double ddP_dmdot;      // d(dP)/dmdot [Pa·s/kg]
    double ddP_dvelocity;  // d(dP)/d(velocity) [Pa·s/m] - chain from this, not
                           // ddP_dmdot, which uses the correlation's own flow area
    double ddP_dT;         // d(dP)/dT [Pa/K]
    double dT_aw_dmdot;    // dT_aw/dmdot [K·s/kg]
    double dT_aw_dT;       // dT_aw/dT [-] (approx 1 at low Mach)
    double dq_dmdot;       // dq/dmdot [W·s/(m²·kg)]
    double dq_dT;          // dq/dT [W/(m²·K)]
    double dq_dT_hot;     // dq/dT_hot [W/(m²·K)]
};

struct WallCouplingResult {
    double Q;              // Heat transfer rate [W]
    double dQ_dh_a;        // ∂Q/∂h_a [m²·K/W]
    double dQ_dh_b;        // ∂Q/∂h_b [m²·K/W]
    double dQ_dT_aw_a;     // ∂Q/∂T_aw_a [W/K]
    double dQ_dT_aw_b;     // ∂Q/∂T_aw_b [W/K]
};
```

---

## Geometry Utilities (geometry.h)

```cpp
double channel_area(double D);
double hydraulic_diameter(double A, double P_wetted);
double channel_roughness(const std::string& material);
double residence_time(double V, double Q);
```

---

## Orifice Flow (orifice.h)

### Geometry and State
```cpp
enum class OrificeType {
    SharpThinPlate, ThickPlate, RoundedEntry, Conical, QuarterCircle, UserDefined
};

struct OrificeGeometry {
    double d, D, t, r, bevel;

    double beta() const;          // Diameter ratio d/D [-]
    double area() const;          // Orifice area [m²]
    double t_over_d() const;      // Thickness ratio t/d [-]
    double r_over_d() const;      // Radius ratio r/d [-]
    bool is_valid() const;
};

struct OrificeState {
    double Re_D, dP, rho, mu;

    double Re_d(double beta) const;  // Orifice Reynolds number (based on d)
};
```

### Cd Correlations
```cpp
// NORMED measurement orifice: a standardised plate in a pipe, Cd referenced
// to the TAPPING differential, every member a function of beta = d/D.
enum class MeteringCdCorrelation {
    ReaderHarrisGallagher, Stolz, Miller,   // ISO 5167 family
    Constant, UserFunction
};

double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);

// Fixed Cd from an explicit value. make_correlation(MeteringCdCorrelation::Constant)
// can only return the default, so this is how a chosen constant is pinned.
std::unique_ptr<OrificeCorrelationBase> make_constant_correlation(double Cd);
```

### Discharge holes: a separate selector

A hole in a wall that dumps coolant out of its circuit is not a metering
orifice: there is no pipe to form `beta` with, and `Cd` depends on `L/d`,
`r/d` and the approach crossflow instead. It gets its own selector and its own
geometry/state pair.

```cpp
enum class DischargeCdCorrelation {
    McGreehanSchotsch1988,   // cooling hole with inlet crossflow
    Idelchik1966Sharp,       // diagram 4-17, Re 25..1e6
    Idelchik1966Thick,       // diagram 4-18a, deep hole
    Idelchik1966Beveled,     // diagram 4-18b
    Idelchik1966Rounded,     // diagram 4-18c
    Lichtarowicz1965,        // long orifice, l/d 2-10, Re 10 to 2e4
    Constant
};

struct DischargeHoleGeometry {
    double d, L, r, bevel;   // [m]; no pipe diameter, by design
    double L_over_d() const;
    double r_over_d() const;
    double bevel_over_d() const;
    double area() const;
    bool is_valid() const;
};

struct DischargeHoleState {
    double Re;            // based on the hole diameter
    double U1_over_Vi;    // INLET tangential velocity ratio; 0 for a plenum
};

class DischargeCorrelationBase {
public:
    virtual double Cd(const DischargeHoleGeometry&, const DischargeHoleState&) const = 0;
    // (Cd, dCd/dRe, dCd/d(U1_over_Vi)), analytic
    virtual std::tuple<double, double, double> Cd_and_derivatives(
        const DischargeHoleGeometry&, const DischargeHoleState&) const = 0;
    virtual std::string name() const = 0;
};

std::unique_ptr<DischargeCorrelationBase> make_discharge_correlation(DischargeCdCorrelation id);
std::unique_ptr<DischargeCorrelationBase> make_constant_discharge_correlation(double Cd);
```

For the Idelchik members `zeta` is referenced to the hole velocity and carries
the full permanent loss, so `Cd = 1/sqrt(zeta)` exactly. See
`validation/cooling/extractions/idelchik_1966_wall_orifice.md` for the tables,
the cross-source agreement with McGreehan-Schotsch, and why the in-a-pipe
diagrams are a different conversion.

### Plenum-to-plenum holes: McGreehan and Schotsch (1988)

A different configuration from the ISO 5167 family above -- a cooling transfer
hole or jet-plate hole, not a metering orifice in a pipe run. ASME J.
Turbomachinery 110(2), 213-217, Eqs. (8)-(17). Each stage is exposed
separately so a check binds to one equation rather than the whole chain.

```cpp
namespace orifice::mcgreehan_schotsch {
double reynolds_baseline(double Re);                 // Eq. (8)
double nozzle_baseline(double Re);                   // Eq. (9)
double corner_factor(double r_over_d);               // Eq. (12), f
double length_factor(double L_over_d);               // Eq. (14), g
double cd_with_corner(double Re, double r_over_d);   // Eq. (11)
double cd_with_corner_and_length(double Re, double r_over_d,
                                 double L_over_d);   // Eqs. (13), (15), (16)
double cd(double Re, double r_over_d, double L_over_d, double U1_over_Vi,
          double eps = rv_smooth_eps);               // Eq. (17)
double cd_with_crossflow(double cd_base, double U1_over_Vi,
                         double eps = rv_smooth_eps);  // Eq. (17) alone

// Solver-facing form: (Cd, dCd/dRe, dCd/d(U1/Vi)). r/d and L/d are geometry
// and carry no partials. Two numerical treatments, both stated in the header:
// the crossflow input is regularised by rv_smooth_eps (Eq. (17) is unbounded
// in slope at U1/Vi = 0, which is the default), and below re_min the value
// stays floored while the derivative is continued from the floor.
std::tuple<double, double, double> cd_and_derivatives(double Re,
                                                      double r_over_d,
                                                      double L_over_d,
                                                      double U1_over_Vi);

// Adiabatic expansion factor, for the INCOMPRESSIBLE form only.
// regime='compressible' already solves the isentropic nozzle exactly via
// nozzle_flow; Eq. (5) reproduces that to 0.008%, so using both double-counts.
double critical_pressure_ratio(double gamma);
double expansion_orifice(double S, double gamma);                 // Eq. (4)
double expansion_nozzle(double S, double gamma);                  // Eq. (5)
double expansion_blend_weight(double cd, double eps = x_smooth_eps);  // Eq. (7)
double expansion_factor(double cd, double S, double gamma,
                        double eps = x_smooth_eps);               // Eq. (6)
}

// The full chain
double orifice::Cd_McGreehanSchotsch(double Re, double r_over_d,
                                     double L_over_d, double U1_over_Vi);
```

`U1_over_Vi` is the INLET (approach, supply-side) tangential velocity ratio,
0 for a plenum-fed plate. It is not a discharge-side crossflow ratio; see the
header and `validation/cooling/extractions/orifice_discharge_coefficient.md`.
`Cd` is not monotonic in it, and is held constant below `Re = 1e4`, its stated
validity floor -- both deliberate, both from the source.

### Film cooling: Baldauf et al. (2002)

Laterally averaged adiabatic film-cooling effectiveness downstream of one row
of **cylindrical**, streamwise-inclined holes. Valid from the ejection point
rather than only far downstream, and it carries the **adjacent jet
interaction** -- the lateral hole-spacing effect driving jet lift-off -- as a
correlated parameter rather than an excluded case.

```cpp
// eta = (T_G - T_AW)/(T_G - T_C), the same convention as
// adiabatic_wall_temperature(), so the two compose directly.
double film_effectiveness_baldauf_2002(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());

// Solver-facing (f, J): (eta, d eta/dM, d eta/dP), analytic via dual numbers.
std::tuple<double, double, double>
film_effectiveness_baldauf_2002_and_derivatives(
    double x_over_D, double M, double P, double alpha_deg, double s_over_D,
    double Tu, CorrelationStatus *status = nullptr,
    double b_0_override = std::numeric_limits<double>::quiet_NaN());
```

`b_0_override` replaces Eq. (31)'s `b_0`; NaN (the default) uses Eq. (31) as
printed. It exists so the paper's Table 4 reading can be *tested* -- the
library offers no second formula because the paper contains none. **The
discrepancy is governed by `x/D`, not by `M`**: `b_0` reaches the model only
through `b_1 = b_0 / (1 + M^-3)` (Eq. 32), the descending-branch gradient, so
it is under 1% at every `M` for `x/D <= 20` and reaches -28% by `x/D = 400`
(worst over the envelope, -75%). At effusion row spacings no effusion dataset
can arbitrate it.

Stated envelope (`baldauf2002::`): `M` 0.2-2.5, `P` 1.2-1.8, `s/D` 2-5,
`alpha` 30-90 deg, `Tu` 0.0035-0.075. The paper's own RMS deviation is 5.5%.

**Reported, not enforced.** Outside the envelope both functions still answer,
finitely and with an unchanged value -- a network solve transits odd states
during Newton iteration, and refusing there would break convergence rather
than protect anyone. They now say so, following
`nusselt_dittus_boelter`'s contract: pass a `CorrelationStatus*` to take the
flag silently, or leave it null to get a warning per out-of-range input
through the global warning handler.

`alpha_deg` is the ejection angle to the **surface** and is converted
internally; the paper's trigonometry is in radians. That is not an
assumption -- seven of its Table 4 coefficients reproduce on radians and fail
on degrees by 10-30%.

**Eq. (31) is implemented as printed and contradicts the paper's own
Table 4** by 36% in `b_0`. Under 5% effect below `M ~ 0.5`, up to 50% at
`M = 2.5`. See
`validation/cooling/extractions/baldauf_2002_film_effectiveness.md`.

### Effusion plate internal heat transfer

Andrews 86-GT-225. A coolant hole cools its wall in two places and the
paper's point is that the approach dominates at the Reynolds numbers
effusion runs at.

```cpp
// Mills' entry-length factor, Eqs. (13) and (14). Two curve fits meeting
// at L/D = 2; decays to 1 as the hole lengthens.
double mills_entry_length_factor(double L_over_D);

// Sparrow's hole approach, Eq. (18): 0.881 Re^0.476 Pr^(1/3) X/(pi L)
double effusion_approach_nusselt(double Re, double Pr, double X_over_L);

// Mills' short-hole throat, Eq. (12): 0.023 Re^0.8 Pr^(1/3) R_Nu
double effusion_throat_nusselt(double Re, double Pr, double L_over_D);

// Eq. (19): the two summed.
double effusion_internal_nusselt(double Re, double Pr, double X_over_L,
                                 double L_over_D);
```

`Re` is on the **hole diameter**, and every `Nu` is referenced to the
**hole internal area** `pi D L` -- which is what makes the two summable.
For the coefficient per unit **plate** area that an effusion element needs:

```
h_plate = Nu * k / D * A_h / A,   A_h = pi D L,  A = X^2 - pi D^2 / 4
```

That factor is 3.46 for Andrews' plate C, so it is not optional. The
conversion is left to the caller because it needs only geometry.

Scored against Andrews' own Fig. 8 at **-13.5% bias, 14.7% RMSE**, labelled
**accuracy**: the correlations are the 1986 paper and the data the 1988
one, and same lab plus different study is cross-source. See
`validation/cooling/extractions/andrews_effusion_internal_h.md`.

### Multi-row film superposition

```cpp
// Sellers, Gao Eq. (1):  eta = 1 - prod_i (1 - eta_i)
double film_superposition_sellers(const std::vector<double>& eta_rows);

// Gao Eq. (7): Sellers with a per-row mainstream temperature correction.
// alpha_between_rows has n-1 entries; all ones reproduces Sellers exactly.
double film_superposition_corrected(const std::vector<double>& eta_rows,
                                    const std::vector<double>& alpha_between_rows);

// (eta, d eta / d eta_i) for each of the two above. Built from explicit
// partial products rather than dividing the total, so a fully effective row
// (eta_j = 1) gives 0, not NaN.
std::pair<double, std::vector<double>> film_superposition_sellers_and_gradient(
    const std::vector<double>& eta_rows);
std::pair<double, std::vector<double>> film_superposition_corrected_and_gradient(
    const std::vector<double>& eta_rows,
    const std::vector<double>& alpha_between_rows);

// Gao Eq. (5): the published FORM of alpha. a and b are REQUIRED -- the
// published values 12 and 0.9465 (gao_a_case1/b_case1) do not transfer:
// alpha >= b for every r, so they cannot reach the 0.69-0.85 a tighter
// plate needs.
double mainstream_temperature_correction(double mass_flow_ratio, double a, double b);

double equivalent_slot_width(double hole_area, double pitch);           // Eq. (9)
double equivalent_blowing_ratio(double M0, double A0, double Ae);        // Eq. (10)
```

Both superposition forms also have `..._and_gradient` variants returning
`(eta, d eta/d eta_i)`, built from partial products so a fully effective row
does not divide by zero.

**Sellers overestimates, and worsens with row count** -- Gao's measurement,
and the reason effusion cannot simply reuse a few-row film model. `alpha` is
the correction: an energy balance on mainstream entrained at each injection
(Gao Eqs. 3-4), bounded in [0,1], equal to 1 for uncorrected Sellers. See
`validation/cooling/extractions/film_superposition.md`.

---

## Acoustics (acoustics.h)

```cpp
std::vector<double> tube_axial_modes(const Tube& tube, double c,
                                     BoundaryCondition bc1, BoundaryCondition bc2,
                                     int n_max);

double helmholtz_frequency(double V, double A_neck, double L_neck, double c,
                           double end_correction = 0.85);

double acoustic_impedance(double rho, double c);
double sound_pressure_level(double p_rms, double p_ref = 20e-6);
```

---

## Humid Air (humidair.h)

```cpp
double humidity_ratio(double T, double P, double RH);
double dewpoint(double T, double P, double RH);
double humid_air_density(double T, double P, double RH);

class HumidAir {
public:
    void set_TP_RH(double T, double P, double RH);
    double rh() const;
    double dewpoint() const;
    State& state();
};
```

---

## Network Solver Interface (solver_interface.h)

Fast-path native interface for network solvers with combined residual and Jacobian evaluations.

### Result Types
```cpp
enum class CorrelationValidity : std::uint8_t { VALID, EXTRAPOLATED, INVALID };

template <typename T> struct CorrelationResult {
    T result;
    CorrelationValidity status;
    std::string message;
};

struct Stream {
    double m_dot;
    double T;
    double P_total;
    std::vector<double> Y;
};

struct StreamJacobian {
    double d_mdot;
    double d_T;
    double d_P_total;
    std::vector<double> d_Y;
};

struct OrificeResult {
    double m_dot_calc;
    double d_mdot_dP_total_up;
    double d_mdot_dP_static_down;
    double d_mdot_dT_up;
    std::vector<double> d_mdot_dY_up;
};

struct ChannelResult {
    double dP_calc;
    double d_dP_d_mdot;
    double d_dP_dP_static_up;
    double d_dP_dT_up;
    std::vector<double> d_dP_dY_up;
};

struct MixerResult {
    double T_mix;
    double P_total_mix;
    std::vector<double> Y_mix;
    double dT_mix_d_delta_h;
    std::vector<StreamJacobian> dT_mix_d_stream;
    std::vector<StreamJacobian> dP_total_mix_d_stream;
    std::vector<std::vector<StreamJacobian>> dY_mix_d_stream;
};

struct AdiabaticResult {
    double T_mix;
    double P_total_mix;
    std::vector<double> Y_mix;
    std::vector<StreamJacobian> dT_mix_d_stream;
    std::vector<StreamJacobian> dP_total_mix_d_stream;
    std::vector<std::vector<StreamJacobian>> dY_mix_d_stream;
};
```

### Incompressible Flow Components
```cpp
// Orifice flow with discharge coefficient
OrificeResult orifice_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down, double Cd,
    double area, double beta = 0.0);

// Channel flow with Darcy friction
ChannelResult channel_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down, double L, double D,
    double roughness, const std::string& friction_model = "haaland");

// Flow restriction (K-factor)
ChannelResult restriction_residuals_and_jacobian(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double P_static_down, double K);
```

### Compressible Flow Components
```cpp
// Compressible orifice flow using isentropic nozzle model
std::tuple<double, double, double, double> orifice_compressible_mdot_and_jacobian(
    double T0, double P0, double P_back, const std::vector<double>& X,
    double Cd, double area, double beta);

// Full compressible orifice evaluation for network solver
OrificeResult orifice_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up, const std::vector<double>& Y_up,
    double P_static_down, double Cd, double area, double beta);

// Compressible channel flow using Fanno model with variable friction
std::tuple<double, double, double, double> channel_compressible_mdot_and_jacobian(
    double T_in, double P_in, double u_in, const std::vector<double>& X,
    double L, double D, double roughness, const std::string& friction_model);

// Full compressible channel evaluation for network solver
ChannelResult channel_compressible_residuals_and_jacobian(
    double m_dot, double P_total_up, double T_up, const std::vector<double>& Y_up,
    double P_static_down, double L, double D, double roughness,
    const std::string& friction_model);
```

### Tee Junction Components (tee_junction.h + solver_interface.h)
```cpp
// Constants
constexpr double TEE_Q_LO     = 0.0;        // lower bound of validated flow ratio
constexpr double TEE_Q_HI     = 1.0;        // upper bound of validated flow ratio
constexpr double TEE_PSI_MIN  = 0.05;       // minimum area ratio psi = A_branch / A_com
constexpr double TEE_THETA_MAX = M_PI / 2;  // maximum branch angle

// Validity check (non-throwing)
struct TeeInputStatus {
    bool q_in_range;
    bool psi_valid;
    bool theta_valid;
    bool valid() const;
    CorrelationStatus status() const;
};
TeeInputStatus tee_check_inputs(double q, double psi, double theta);

// Raw K-coefficient functions (Bassett 2001, Table 2 + Eq 33/34 angle corrections).
// Pure math, no input guards. K2 and K5 are psi/theta-independent per Eq 15.
double K1(double q, double psi, double theta);
double K2(double q);                                   // = q^2 - 1.5*q + 0.5
double K3(double q, double psi, double theta);
double K4(double q, double psi, double theta);
double K5(double q);                                   // same as K2; separating type 3
double K6(double q, double psi, double theta);
double K7(double q, double psi, double theta);
double K8(double q, double psi, double theta);
double K9(double q, double psi, double theta);
double K10(double q, double psi, double theta);
double K11(double q, double psi, double theta);
double K12(double q, double psi, double theta);
// Analytic q-derivatives for each.
double dK1_dq(double q, double psi, double theta);
double dK3_dq(double q, double psi, double theta);
double dK4_dq(double q, double psi, double theta);
double dK5_dq(double q);
double dK6_dq(double q, double psi, double theta);
double dK7_dq(double q, double psi, double theta);
double dK8_dq(double q, double psi, double theta);
double dK9_dq(double q, double psi, double theta);
double dK10_dq(double q, double psi, double theta);
double dK11_dq(double q, double psi, double theta);
double dK12_dq(double q, double psi, double theta);

// Smooth blend helpers
double soft_lower(double x, double lo);              // smooth max(x, lo)
double blend_weight(double r, double k = 30.0);      // 0.5*(1+tanh(k*r))
double blend_weight_deriv(double r, double k = 30.0); // d(blend_weight)/dr

// Blended K functions (safe for any real inputs, never NaN)
double merging_tee_K_straight(double q, double psi, double theta, double blend_k = 30.0);
double merging_tee_K_branch(double q, double psi, double theta, double blend_k = 30.0);
double branching_tee_K_straight(double q, double psi, double theta, double blend_k = 30.0);
double branching_tee_K_branch(double q, double psi, double theta, double blend_k = 30.0);

// Helper utilities
double tee_flow_ratio(double m_dot_branch, double m_dot_com);
bool   tee_topology_valid(double q, double epsilon = 0.05);

// Solver result struct (declared in solver_interface.h)
struct TeeJunctionResult {
    double R_straight, R_branch;
    double dR_straight_d_mdot_com, dR_straight_d_mdot_branch;
    double dR_straight_dP_static_com, dR_straight_dT_com;
    std::vector<double> dR_straight_dY_com;
    double dR_branch_d_mdot_com, dR_branch_d_mdot_branch;
    double dR_branch_dP_static_com, dR_branch_dT_com;
    std::vector<double> dR_branch_dY_com;
    double K_straight, K_branch, q, blend_w;
    bool topology_valid;
    CorrelationValidity status;
};

// Solver interface
TeeJunctionResult merging_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com, const std::vector<double>& Y_com,
    double theta, double psi, double F_C, double blend_k = 30.0);

TeeJunctionResult branching_tee_residuals_and_jacobian(
    double m_dot_com, double m_dot_branch,
    double dP0_straight, double dP0_branch,
    double P_static_com, double T_com, const std::vector<double>& Y_com,
    double theta, double psi, double F_C, double blend_k = 30.0);
```

Port convention: MAIN_INLET=A, BRANCH=B, MAIN_OUTLET=C; `m_dot_com` = common-port
mass flow, `m_dot_branch` = branch mass flow, `F_C` = cross-sectional area [m^2] at
the common port, `psi` = A_branch / A_com, `theta` = branch angle [rad].

### Momentum-CV Junction: whole-element (f, J)

Backs `MultiPortChamberElement`. Unlike the tee functions above, this returns the whole
element's residual vector and its full Jacobian from one seeded evaluation --
forward-mode dual numbers over every unknown, so there is no separate
derivation to keep in sync (`include/mpce_junction.h`).

```cpp
struct MpceGeometry {
    std::array<double, 3> area, theta_rad, port_sign;
    double joining_etransfer_alpha, eta_scale;
};

struct MpceResidualJacobian {
    std::array<double, 4> residual;                    // 3 port rows + mass
    std::array<std::array<double, 10>, 4> jacobian;    // [row][seed]
    bool valid;                                        // false: not a junction
    int common_port;
    double k_term_sign;
    std::array<double, 3> k_per_port;
};

MpceResidualJacobian mpce_residuals_and_jacobian(
    const std::array<double, 3>& p_static, const std::array<double, 3>& p_total,
    const std::array<double, 3>& rho,      const std::array<double, 3>& drho_dp,
    const std::array<double, 3>& outer_mdot, double pt_jct,
    const MpceGeometry& geom);
```

Seed order, matching the Jacobian's columns:

| columns | unknown |
|---|---|
| 0-2 | `p_static[i]` |
| 3-5 | `p_total[i]` |
| 6-8 | `outer_mdot[i]` |
| 9 | `pt_jct` |

`outer_mdot` is the CONNECTING element's mass flow, not the junction-convention
one; `geom.port_sign` maps between them, so the Jacobian arrives already
expressed in the solver's own unknowns.

**Physics only.** The kernel does not own the degenerate-state guards or the
wrong-direction soft barrier -- those are solver policy and live in
`MultiPortChamberElement`. It reports `valid = false` for a flow pattern that is not a
junction in any regime rather than guessing.

**Thermodynamics stays at the call site**: pass `rho` and `drho/dP` per port.
A Jacobian is first order, so seeding `rho + (drho/dP)(P - P_val)` is exact.
For combaero's mixtures `drho/dP = rho/P` holds to 1e-11.

Residual sign convention: `R_straight = dP0_straight + K_straight * q_dyn`,
`R_branch = dP0_branch + K_branch * q_dyn`.

### Multi-Port Chamber (Momentum-CV Junction)

Sanctioned successor to the K-closure tee for N-port junctions
(`docs/junction/momentum cv implementation guide.pdf`). Junction = pure
conservation, loss = separate per-port `BorderCarnotLossElement`s.

> **The junction half of this was removed in 0.6.0.** `multi_port_chamber_residuals_and_jacobian` and its result structs backed `MultiPortChamberBase`'s own impulse model, which is gone; `MultiPortChamberElement` supersedes it and computes its own whole-element `(f, J)` (see the Momentum-CV Junction section above). What remains here is the Border-Carnot loss element, which is unaffected.

```cpp
// Loss-element constants (multi_port_chamber.h)
inline constexpr double HAGER_FRACTION       = 0.75;  // sharp-edge angle correction
inline constexpr double BC_LOSS_PREFACTOR    = 4.0;   // L = 4*(1 - cos(theta_eff))^2

// Loss coefficient L = 4*(1 - cos((3/4)*delta_geom))^2  [-]
double border_carnot_L(double delta_geom);
double dborder_carnot_L_ddelta(double delta_geom);

// Result structs (solver_interface.h)
// Two-port in-line loss element: Pt_in - Pt_out - L*0.5*rho*u_in^2 = 0.
// L applies the Hager (3/4) effective-angle correction, INTENDED to reproduce
// Hager xi_l and Bassett K_inc at M -> 0 on a sharp-edged lateral -- unverified,
// Tier-1 tests xfail at 11-29% deviation (issue #272). Sign-free in mdot
// (mdot^2 in the dynamic head).
BorderCarnotLossResult border_carnot_loss_residual_and_jacobian(
    double mdot,
    double Pt_in, double Pt_out,
    double P_in, double T_in,
    const std::vector<double>& Y_in,
    double area, double delta_geom);
```

DOF accounting per junction: +1 unknown (`P_jct`), +N+1 residuals = net +N
equations. Combined with port-MCN mass-row skip in the solver (the junction's
sum-mass residual replaces the N port-MCN continuity rows), the system stays
square.

### Mixing and Combustion Components
```cpp
// Stream mixing with optional heat transfer
MixerResult mixer_from_streams_and_jacobians(
    const std::vector<Stream>& streams, double Q = 0.0, double fraction = 0.0);

// Adiabatic complete combustion
AdiabaticResult adiabatic_T_complete_and_jacobian_T_from_streams(
    const std::vector<Stream>& streams, double P, double Q = 0.0, double fraction = 0.0);

// Adiabatic equilibrium combustion
AdiabaticResult adiabatic_T_equilibrium_and_jacobians_from_streams(
    const std::vector<Stream>& streams, double P, double Q = 0.0, double fraction = 0.0);

// Combined combustion with analytical pressure loss
ChamberResult combustor_residuals_and_jacobians(
    double m_dot, double P_total_up, double P_static_up, double T_up,
    const std::vector<double>& Y_up, double Q_comb,
    CombustionMethod method, bool smooth,
    const PressureLossCorrelation& pressure_loss);
```
