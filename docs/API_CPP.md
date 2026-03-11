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
- [Geometry Utilities (geometry.h)](#geometry-utilities-geometryh)
- [Orifice Flow (orifice.h)](#orifice-flow-orificeh)
- [Acoustics (acoustics.h)](#acoustics-acousticsh)
- [Humid Air (humidair.h)](#humid-air-humidairh)

---

## Thermodynamics (thermo.h)

### Species and Mixture Metadata
```cpp
// Get species metadata
std::string species_name(std::size_t index);
std::size_t species_index_from_name(const std::string& name);
double species_molar_mass(std::size_t index);
double species_molar_mass_from_name(const std::string& name);
std::size_t num_species();

// Mass/Mole fraction conversions
double mwmix(const std::vector<double>& X);
std::vector<double> mole_to_mass(const std::vector<double>& X);
std::vector<double> mass_to_mole(const std::vector<double>& Y);
std::vector<double> normalize_fractions(const std::vector<double>& X);
```

### Thermodynamic Properties (Molar Basis)
```cpp
double cp(double T, const std::vector<double>& X);
double cv(double T, const std::vector<double>& X);
double h(double T, const std::vector<double>& X);
double u(double T, const std::vector<double>& X);
double s(double T, double P, const std::vector<double>& X, double P_ref = 101325.0);
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
double calc_T_from_h(double h_target, const std::vector<double>& X);
double calc_T_from_s(double s_target, double P, const std::vector<double>& X);
```

---

## Transport Properties (transport.h)

```cpp
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);
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
                                     double A_eff, const std::vector<double>& X);

FannoSolution fanno_pipe(double T_in, double P_in, double u_in, double L, double D,
                         double f, const std::vector<double>& X);

FannoSolution fanno_pipe_rough(double T_in, double P_in, double u_in, double L, double D,
                               double roughness, const std::vector<double>& X,
                               const std::string& correlation = "haaland");
```

---

## Incompressible Flow (incompressible.h)

```cpp
double bernoulli_P2(double P1, double v1, double v2, double rho, double dz = 0.0);
double orifice_mdot(double P1, double P2, double A, double Cd, double rho);
double pipe_dP(double v, double L, double D, double f, double rho);
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

---

## Geometry Utilities (geometry.h)

```cpp
double pipe_area(double D);
double hydraulic_diameter(double A, double P_wetted);
double pipe_roughness(const std::string& material);
double residence_time(double V, double Q);
```

---

## Orifice Flow (orifice.h)

```cpp
struct OrificeGeometry {
    double d, D, t, r;
    double beta() const;
    double area() const;
};

struct OrificeState {
    double Re_D, dP, rho, mu;
};

double Cd_sharp_thin_plate(const OrificeGeometry& geom, const OrificeState& state);
double Cd_thick_plate(const OrificeGeometry& geom, const OrificeState& state);
double Cd_rounded_entry(const OrificeGeometry& geom, const OrificeState& state);
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
