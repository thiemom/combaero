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

---

## Thermodynamics (thermo.h)

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

```cpp
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);
double reynolds_from_state(double rho, double v, double L, double mu);
```

---

## Combustion (combustion.h)

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
```

---

## Heat Transfer Correlations (heat_transfer.h)

### Nusselt Numbers
```cpp
double nusselt_dittus_boelter(double Re, double Pr, bool heating = true);
double nusselt_gnielinski(double Re, double Pr);
double nusselt_gnielinski(double Re, double Pr, double f);
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

// Enhanced surfaces
ChannelResult channel_ribbed(double T, double P, const std::vector<double>& X,
                             double velocity, double diameter, double length,
                             double e_D, double p_e, double w_e,
                             double T_hot = std::numeric_limits<double>::quiet_NaN(),
                             double Nu_multiplier = 1.0, double f_multiplier = 1.0);

ChannelResult channel_dimpled(double T, double P, const std::vector<double>& X,
                              double velocity, double diameter, double length,
                              double d_Dh, double h_d, double S_d,
                              double T_hot = std::numeric_limits<double>::quiet_NaN(),
                              double Nu_multiplier = 1.0, double f_multiplier = 1.0);

ChannelResult channel_pin_fin(double T, double P, const std::vector<double>& X,
                              double velocity, double diameter, double length,
                              double L_H, double S_H, double S_L, double t_D,
                              double T_hot = std::numeric_limits<double>::quiet_NaN(),
                              double Nu_multiplier = 1.0, double f_multiplier = 1.0);

ChannelResult channel_impingement(double T, double P, const std::vector<double>& X,
                                  double velocity, double diameter, double length,
                                  double H_D, double S_n, double d_j, double N_j,
                                  double T_hot = std::numeric_limits<double>::quiet_NaN(),
                                  double Nu_multiplier = 1.0, double f_multiplier = 1.0);
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
enum class CdCorrelation {
    // Sharp thin-plate
    ReaderHarrisGallagher, Stolz, Miller,
    // Thick-plate corrections
    IdelchikThick, BohlThick,
    // Rounded-entry
    IdelchikRounded, BohlRounded,
    // Special
    Constant, UserFunction
};

// Individual correlations
double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);
double Cd_thick_plate(const OrificeGeometry& geom, const OrificeState& state);
double Cd_rounded_entry(const OrificeGeometry& geom, const OrificeState& state);

// Auto-select correlation based on geometry
double Cd(const OrificeGeometry& geom, const OrificeState& state);
```

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
