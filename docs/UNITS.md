# combaero Units Reference

This document defines the SI-based unit system used throughout the combaero library.
All functions use consistent units to avoid conversion errors.

**Auto-generated from `include/units_data.h` - do not edit manually.**

---

## Base SI Units

| Quantity           | Unit   | Symbol |
|--------------------|--------|--------|
| Temperature        | Kelvin | K      |
| Pressure           | Pascal | Pa     |
| Mass               | kilogram | kg   |
| Length             | meter  | m      |
| Time               | second | s      |
| Amount of substance| mole   | mol    |

---

## Derived Units

| Quantity             | Unit        | Symbol     |
|----------------------|-------------|------------|
| Energy               | Joule       | J = kg*m^2/s^2 |
| Power                | Watt        | W = J/s    |
| Force                | Newton      | N = kg*m/s^2 |
| Dynamic viscosity    | Pascal-second | Pa*s     |
| Kinematic viscosity  | -           | m^2/s      |
| Thermal conductivity | -           | W/(m*K)    |
| Specific heat        | -           | J/(mol*K) or J/(kg*K) |
| Enthalpy             | -           | J/mol or J/kg |
| Entropy              | -           | J/(mol*K) or J/(kg*K) |
| Density              | -           | kg/m^3     |
| Velocity             | -           | m/s        |
| Area                 | -           | m^2        |
| Volume               | -           | m^3        |
| Mass flow rate       | -           | kg/s       |
| Volumetric flow rate | -           | m^3/s      |

---

## Physical Constants

| Constant             | Value              | Unit       |
|----------------------|--------------------|------------|
| Universal gas constant R | 8.31446261815324 | J/(mol*K)  |
| Boltzmann constant   | 1.380649e-23       | J/K        |
| Avogadro's number    | 6.02214076e23      | 1/mol      |
| Standard gravity g0  | 9.80665            | m/s^2      |

---

## Function Units by Module

### thermo.h

#### Species Data

| Function             | Input Units | Output Unit |
|----------------------|-------------|-------------|
| `species_molar_mass` | -           | g/mol       |
| `mwmix`              | X: mol/mol  | g/mol       |
| `mole_to_mass`       | X: mol/mol  | kg/kg       |
| `mass_to_mole`       | Y: kg/kg    | mol/mol     |

#### Dimensionless NASA Polynomials

| Function    | Input Units | Output Unit |
|-------------|-------------|-------------|
| `cp_R`      | T: K        | - (Cp/R)    |
| `h_RT`      | T: K        | - (H/(R*T)) |
| `s_R`       | T: K        | - (S/R)     |
| `g_over_RT` | T: K        | - (G/(R*T)) |

#### Per-Species Properties

| Function     | Input Units | Output Unit |
|--------------|-------------|-------------|
| `cp_species` | T: K        | J/(mol*K)   |
| `h_species`  | T: K        | J/mol       |
| `s_species`  | T: K        | J/(mol*K)   |

#### Mixture Properties (Molar Basis)

| Function | Input Units                        | Output Unit |
|----------|------------------------------------|-------------|
| `cp`     | T: K, X: mol/mol                   | J/(mol*K)   |
| `cv`     | T: K, X: mol/mol                   | J/(mol*K)   |
| `h`      | T: K, X: mol/mol                   | J/mol       |
| `u`      | T: K, X: mol/mol                   | J/mol       |
| `s`      | T: K, X: mol/mol, P: Pa, P_ref: Pa | J/(mol*K)   |
| `dh_dT`  | T: K, X: mol/mol                   | J/(mol*K)   |
| `ds_dT`  | T: K, X: mol/mol                   | J/(mol*K^2) |
| `dcp_dT` | T: K, X: mol/mol                   | J/(mol*K^2) |

#### Mixture Properties (Mass/Other Basis)

| Function                           | Input Units                                                                                | Output Unit            |
|------------------------------------|--------------------------------------------------------------------------------------------|------------------------|
| `density`                          | T: K, P: Pa, X: mol/mol                                                                    | kg/m^3                 |
| `molar_volume`                     | T: K, P: Pa                                                                                | m^3/mol                |
| `specific_gas_constant`            | X: mol/mol                                                                                 | J/(kg*K)               |
| `isentropic_expansion_coefficient` | T: K, X: mol/mol                                                                           | - (gamma)              |
| `speed_of_sound`                   | T: K, X: mol/mol                                                                           | m/s                    |
| `cp_mass`                          | T: K, X: mol/mol                                                                           | J/(kg*K)               |
| `cv_mass`                          | T: K, X: mol/mol                                                                           | J/(kg*K)               |
| `h_mass`                           | T: K, X: mol/mol                                                                           | J/kg                   |
| `s_mass`                           | T: K, X: mol/mol, P: Pa, P_ref: Pa                                                         | J/(kg*K)               |
| `u_mass`                           | T: K, X: mol/mol                                                                           | J/kg                   |
| `air_properties`                   | T: K, P: Pa, humidity: - (0-1)                                                             | AirProperties struct   |
| `thermo_state`                     | T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)                                        | ThermoState struct     |
| `transport_state`                  | T: K, P: Pa, X: mol/mol                                                                    | TransportState struct  |
| `complete_state`                   | T: K, P: Pa, X: mol/mol, P_ref: Pa (default 101325)                                        | CompleteState struct   |
| `combustion_state`                 | X_fuel: mol/mol, X_ox: mol/mol, phi: -, T_reactants: K, P: Pa, fuel_name: str (default '') | CombustionState struct |
| `combustion_state_from_streams`    | fuel_stream: Stream, ox_stream: Stream, fuel_name: str (default '')                        | CombustionState struct |

#### Inverse Solvers

| Function         | Input Units                            | Output Unit |
|------------------|----------------------------------------|-------------|
| `calc_T_from_h`  | h_target: J/mol, X: mol/mol            | K           |
| `calc_T_from_s`  | s_target: J/(mol*K), P: Pa, X: mol/mol | K           |
| `calc_T_from_cp` | cp_target: J/(mol*K), X: mol/mol       | K           |

---

### transport.h

#### Transport Properties

| Function               | Input Units                           | Output Unit |
|------------------------|---------------------------------------|-------------|
| `viscosity`            | T: K, P: Pa, X: mol/mol               | Pa*s        |
| `thermal_conductivity` | T: K, P: Pa, X: mol/mol               | W/(m*K)     |
| `prandtl`              | T: K, P: Pa, X: mol/mol               | - (Pr)      |
| `kinematic_viscosity`  | T: K, P: Pa, X: mol/mol               | m^2/s       |
| `thermal_diffusivity`  | T: K, P: Pa, X: mol/mol               | m^2/s       |
| `reynolds`             | T: K, P: Pa, X: mol/mol, V: m/s, L: m | - (Re)      |
| `peclet`               | T: K, P: Pa, X: mol/mol, V: m/s, L: m | - (Pe)      |

---

### compressible.h

#### Nozzle Flow

| Function                   | Input Units                                         | Output Unit              |
|----------------------------|-----------------------------------------------------|--------------------------|
| `nozzle_flow`              | T0: K, P0: Pa, P_back: Pa, A_eff: m^2, X: mol/mol   | CompressibleFlowSolution |
| `nozzle_quasi1d`           | T0: K, P0: Pa, P_exit: Pa, x: m, A: m^2, X: mol/mol | NozzleSolution           |
| `nozzle_cd`                | T0: K, P0: Pa, P_exit: Pa, A: m^2, x: m, X: mol/mol | NozzleSolution           |
| `critical_pressure_ratio`  | T0: K, P0: Pa, X: mol/mol                           | - (P*/P0)                |
| `mach_from_pressure_ratio` | T0: K, P0: Pa, P: Pa, X: mol/mol                    | - (M)                    |
| `mass_flux_isentropic`     | T0: K, P0: Pa, P: Pa, X: mol/mol                    | kg/(m^2*s)               |

#### Fanno Flow

| Function           | Input Units                                                                          | Output Unit   |
|--------------------|--------------------------------------------------------------------------------------|---------------|
| `fanno_pipe`       | T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, f: -, X: mol/mol                           | FannoSolution |
| `fanno_pipe_rough` | T_in: K, P_in: Pa, u_in: m/s, L: m, D: m, roughness: m, X: mol/mol, correlation: str | FannoSolution |
| `fanno_max_length` | T_in: K, P_in: Pa, u_in: m/s, D: m, f: -, X: mol/mol                                 | m             |

#### Thrust

| Function        | Input Units               | Output Unit  |
|-----------------|---------------------------|--------------|
| `nozzle_thrust` | NozzleSolution, P_amb: Pa | ThrustResult |

---

### incompressible.h

#### Thermo-Aware High-Level Functions

| Function              | Input Units                                                                 | Output Unit                |
|-----------------------|-----------------------------------------------------------------------------|----------------------------|
| `pipe_flow`           | T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, f: -                           | IncompressibleFlowSolution |
| `pipe_flow_rough`     | T: K, P: Pa, X: mol/mol, u: m/s, L: m, D: m, roughness: m, correlation: str | IncompressibleFlowSolution |
| `orifice_flow_thermo` | T: K, P: Pa, X: mol/mol, P_back: Pa, A: m^2, Cd: -                          | IncompressibleFlowSolution |
| `pressure_drop_pipe`  | T: K, P: Pa, X: mol/mol, v: m/s, D: m, L: m, roughness: m, correlation: str | tuple(Pa, -, -)            |

#### Bernoulli & Orifice

| Function               | Input Units                                            | Output Unit |
|------------------------|--------------------------------------------------------|-------------|
| `bernoulli_P2`         | P1: Pa, v1: m/s, v2: m/s, rho: kg/m^3, dz: m, g: m/s^2 | Pa          |
| `bernoulli_v2`         | P1: Pa, P2: Pa, v1: m/s, rho: kg/m^3, dz: m, g: m/s^2  | m/s         |
| `orifice_mdot`         | P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3             | kg/s        |
| `orifice_Q`            | P1: Pa, P2: Pa, A: m^2, Cd: -, rho: kg/m^3             | m^3/s       |
| `orifice_velocity`     | P1: Pa, P2: Pa, rho: kg/m^3                            | m/s         |
| `orifice_area`         | mdot: kg/s, P1: Pa, P2: Pa, Cd: -, rho: kg/m^3         | m^2         |
| `orifice_dP`           | mdot: kg/s, A: m^2, Cd: -, rho: kg/m^3                 | Pa          |
| `expansibility_factor` | beta: -, dP: Pa, P_upstream: Pa, kappa: -              | -           |

---

### friction.h

#### Friction Factor Correlations

| Function             | Input Units   | Output Unit |
|----------------------|---------------|-------------|
| `friction_haaland`   | Re: -, e_D: - | - (f)       |
| `friction_serghides` | Re: -, e_D: - | - (f)       |
| `friction_colebrook` | Re: -, e_D: - | - (f)       |
| `friction_petukhov`  | Re: -         | - (f)       |

---

### orifice.h

#### Discharge Coefficients

| Function              | Input Units                   | Output Unit |
|-----------------------|-------------------------------|-------------|
| `Cd_sharp_thin_plate` | OrificeGeometry, OrificeState | - (Cd)      |
| `Cd_thick_plate`      | OrificeGeometry, OrificeState | - (Cd)      |
| `Cd_rounded_entry`    | OrificeGeometry, OrificeState | - (Cd)      |
| `Cd`                  | OrificeGeometry, OrificeState | - (Cd)      |

---

### humidair.h

#### Humid Air Properties

| Function                          | Input Units              | Output Unit |
|-----------------------------------|--------------------------|-------------|
| `saturation_vapor_pressure`       | T: K                     | Pa          |
| `vapor_pressure`                  | T: K, RH: - (0-1)        | Pa          |
| `humidity_ratio`                  | T: K, P: Pa, RH: - (0-1) | kg/kg       |
| `water_vapor_mole_fraction`       | T: K, P: Pa, RH: - (0-1) | mol/mol     |
| `humid_air_composition`           | T: K, P: Pa, RH: - (0-1) | mol/mol     |
| `dewpoint`                        | T: K, P: Pa, RH: - (0-1) | K           |
| `relative_humidity_from_dewpoint` | T: K, Tdp: K, P: Pa      | - (0-1)     |
| `wet_bulb_temperature`            | T: K, P: Pa, RH: - (0-1) | K           |
| `humid_air_enthalpy`              | T: K, P: Pa, RH: - (0-1) | J/kg        |
| `humid_air_density`               | T: K, P: Pa, RH: - (0-1) | kg/m^3      |

---

### combustion.h

#### Stoichiometry

| Function                          | Input Units                               | Output Unit      |
|-----------------------------------|-------------------------------------------|------------------|
| `oxygen_required_per_mol_fuel`    | -                                         | mol O2/mol fuel  |
| `oxygen_required_per_kg_fuel`     | -                                         | mol O2/kg fuel   |
| `oxygen_required_per_mol_mixture` | X: mol/mol                                | mol O2/mol mix   |
| `oxygen_required_per_kg_mixture`  | X: mol/mol                                | mol O2/kg mix    |
| `fuel_lhv_molar`                  | X_fuel: mol/mol, reference_temperature: K | J/mol fuel       |
| `fuel_lhv_mass`                   | X_fuel: mol/mol, reference_temperature: K | J/kg fuel        |
| `dryair_required_per_mol_fuel`    | -                                         | mol air/mol fuel |
| `dryair_required_per_kg_fuel`     | -                                         | mol air/kg fuel  |
| `dryair_required_per_mol_mixture` | X: mol/mol                                | mol air/mol mix  |
| `dryair_required_per_kg_mixture`  | X: mol/mol                                | mol air/kg mix   |

#### Equivalence Ratio

| Function                     | Input Units        | Output Unit |
|------------------------------|--------------------|-------------|
| `equivalence_ratio_mole`     | X: mol/mol         | - (phi)     |
| `set_equivalence_ratio_mole` | phi: -, X: mol/mol | mol/mol     |
| `equivalence_ratio_mass`     | Y: kg/kg           | - (phi)     |
| `set_equivalence_ratio_mass` | phi: -, Y: kg/kg   | kg/kg       |

#### Mixture Fraction (Bilger)

| Function                               | Input Units      | Output Unit |
|----------------------------------------|------------------|-------------|
| `bilger_beta`                          | Y: kg/kg         | -           |
| `bilger_mixture_fraction`              | Y: kg/kg         | - (Z)       |
| `bilger_stoich_mixture_fraction_mass`  | Y: kg/kg         | - (Z_st)    |
| `equivalence_ratio_from_bilger_Z_mass` | Z: -, Y: kg/kg   | - (phi)     |
| `bilger_Z_from_equivalence_ratio_mass` | phi: -, Y: kg/kg | - (Z)       |

#### Complete Combustion

| Function                         | Input Units                     | Output Unit |
|----------------------------------|---------------------------------|-------------|
| `complete_combustion_to_CO2_H2O` | X: mol/mol                      | mol/mol     |
| `complete_combustion`            | State (T: K, P: Pa, X: mol/mol) | State       |
| `complete_combustion_isothermal` | State                           | State       |

#### Stream Solvers

| Function                          | Input Units               | Output Unit |
|-----------------------------------|---------------------------|-------------|
| `set_fuel_stream_for_phi`         | phi: -                    | Stream      |
| `set_fuel_stream_for_Tad`         | T_ad_target: K            | Stream      |
| `set_fuel_stream_for_O2`          | X_O2_target: mol/mol      | Stream      |
| `set_fuel_stream_for_O2_dry`      | X_O2_dry_target: mol/mol  | Stream      |
| `set_fuel_stream_for_CO2`         | X_CO2_target: mol/mol     | Stream      |
| `set_fuel_stream_for_CO2_dry`     | X_CO2_dry_target: mol/mol | Stream      |
| `set_oxidizer_stream_for_Tad`     | T_ad_target: K            | Stream      |
| `set_oxidizer_stream_for_O2`      | X_O2_target: mol/mol      | Stream      |
| `set_oxidizer_stream_for_O2_dry`  | X_O2_dry_target: mol/mol  | Stream      |
| `set_oxidizer_stream_for_CO2`     | X_CO2_target: mol/mol     | Stream      |
| `set_oxidizer_stream_for_CO2_dry` | X_CO2_dry_target: mol/mol | Stream      |

---

### equilibrium.h

#### Chemical Equilibrium

| Function                          | Input Units                     | Output Unit |
|-----------------------------------|---------------------------------|-------------|
| `wgs_equilibrium`                 | State (T: K, P: Pa, X: mol/mol) | State       |
| `wgs_equilibrium_adiabatic`       | State                           | State       |
| `smr_wgs_equilibrium`             | State                           | State       |
| `smr_wgs_equilibrium_adiabatic`   | State                           | State       |
| `reforming_equilibrium`           | State                           | State       |
| `reforming_equilibrium_adiabatic` | State                           | State       |
| `combustion_equilibrium`          | State                           | State       |

---

### state.h

#### State Properties

| Function       | Input Units | Output Unit |
|----------------|-------------|-------------|
| `State::T`     | -           | K           |
| `State::P`     | -           | Pa          |
| `State::X`     | -           | mol/mol     |
| `State::mw`    | -           | g/mol       |
| `State::cp`    | -           | J/(mol*K)   |
| `State::cv`    | -           | J/(mol*K)   |
| `State::h`     | -           | J/mol       |
| `State::u`     | -           | J/mol       |
| `State::s`     | -           | J/(mol*K)   |
| `State::rho`   | -           | kg/m^3      |
| `State::R`     | -           | J/(kg*K)    |
| `State::gamma` | -           | - (gamma)   |
| `State::a`     | -           | m/s         |
| `State::mu`    | -           | Pa*s        |
| `State::k`     | -           | W/(m*K)     |
| `State::nu`    | -           | m^2/s       |
| `State::Pr`    | -           | - (Pr)      |
| `State::alpha` | -           | m^2/s       |

#### Stream Properties

| Function       | Input Units | Output Unit |
|----------------|-------------|-------------|
| `Stream::mdot` | -           | kg/s        |

---

### stagnation.h

#### Stagnation/Static Conversions

| Function                | Input Units                                  | Output Unit |
|-------------------------|----------------------------------------------|-------------|
| `kinetic_energy`        | v: m/s                                       | J/kg        |
| `h0_from_static`        | h_static: J/kg, v: m/s                       | J/kg        |
| `v_from_h0`             | h0: J/kg, h_static: J/kg                     | m/s         |
| `mach_number`           | v: m/s, T: K, X: mol/mol                     | - (M)       |
| `T0_from_static`        | T: K, M: -, X: mol/mol                       | K           |
| `T0_from_static_v`      | T: K, v: m/s, X: mol/mol                     | K           |
| `T_from_stagnation`     | T0: K, M: -, X: mol/mol                      | K           |
| `P0_from_static`        | P: Pa, T: K, M: -, X: mol/mol                | Pa          |
| `P_from_stagnation`     | P0: Pa, T0: K, M: -, X: mol/mol              | Pa          |
| `recovery_factor`       | Pr: -                                        | - (r)       |
| `T_adiabatic_wall`      | T_static: K, v: m/s, T: K, P: Pa, X: mol/mol | K           |
| `T_adiabatic_wall_mach` | T_static: K, M: -, T: K, P: Pa, X: mol/mol   | K           |

### materials.h - Material Thermal Conductivity

| Function                | Input Units                                                | Output Unit |
|-------------------------|------------------------------------------------------------|-------------|
| `k_inconel718`          | T: K                                                       | W/(m*K)     |
| `k_haynes230`           | T: K                                                       | W/(m*K)     |
| `k_stainless_steel_316` | T: K                                                       | W/(m*K)     |
| `k_aluminum_6061`       | T: K                                                       | W/(m*K)     |
| `k_tbc_ysz`             | T: K, hours: h (default 0), is_ebpvd: bool (default False) | W/(m*K)     |
| `list_materials`        | -                                                          | list[str]   |

### cooling_correlations.h - Advanced Cooling Correlations

| Function                          | Input Units                                                               | Output Unit |
|-----------------------------------|---------------------------------------------------------------------------|-------------|
| `rib_enhancement_factor`          | e_D: -, pitch_to_height: -, alpha_deg: deg                                | -           |
| `rib_enhancement_factor_high_re`  | e_D: -, pitch_to_height: -, alpha_deg: deg, Re: -                         | -           |
| `rib_friction_multiplier`         | e_D: -, pitch_to_height: -                                                | -           |
| `rib_friction_multiplier_high_re` | e_D: -, pitch_to_height: -                                                | -           |
| `thermal_performance_factor`      | Nu_ratio: -, f_ratio: -                                                   | -           |
| `impingement_nusselt`             | Re_jet: -, Pr: -, z_D: -, x_D: - (default 0), y_D: - (default 0)          | -           |
| `film_cooling_effectiveness`      | x_D: -, M: -, DR: -, alpha_deg: deg                                       | -           |
| `film_cooling_effectiveness_avg`  | x_D: -, M: -, DR: -, alpha_deg: deg, s_D: - (default 3.0)                 | -           |
| `film_cooling_multirow_sellers`   | row_positions_xD: list[-], eval_xD: -, M: -, DR: -, alpha_deg: deg        | -           |
| `effusion_effectiveness`          | x_D: -, M: -, DR: -, porosity: -, s_D: -, alpha_deg: deg                  | -           |
| `pin_fin_nusselt`                 | Re_d: -, Pr: -, L_D: -, S_D: -, X_D: -, is_staggered: bool (default True) | -           |
| `dimple_nusselt_enhancement`      | Re_Dh: -, d_Dh: -, h_d: -, S_d: -                                         | -           |
| `dimple_friction_multiplier`      | Re_Dh: -, d_Dh: -, h_d: -                                                 | -           |
| `effusion_discharge_coefficient`  | Re_d: -, P_ratio: -, alpha_deg: deg, L_D: - (default 4.0)                 | -           |

### acoustics.h - Acoustic Properties

| Function                    | Input Units                                                                      | Output Unit               |
|-----------------------------|----------------------------------------------------------------------------------|---------------------------|
| `acoustic_properties`       | f: Hz, rho: kg/m^3, c: m/s, p_rms: Pa (default 20e-6), p_ref: Pa (default 20e-6) | AcousticProperties struct |
| `wavelength`                | f: Hz, c: m/s                                                                    | m                         |
| `frequency_from_wavelength` | lambda: m, c: m/s                                                                | Hz                        |
| `acoustic_impedance`        | rho: kg/m^3, c: m/s                                                              | Pa*s/m                    |
| `sound_pressure_level`      | p_rms: Pa, p_ref: Pa (default 20e-6)                                             | dB                        |
| `particle_velocity`         | p: Pa, rho: kg/m^3, c: m/s                                                       | m/s                       |

### orifice.h - Orifice Flow Utilities

| Function                     | Input Units                                                              | Output Unit              |
|------------------------------|--------------------------------------------------------------------------|--------------------------|
| `orifice_flow`               | geom: OrificeGeometry, dP: Pa, T: K, P: Pa, mu: Pa*s, Z: - (default 1.0) | OrificeFlowResult struct |
| `orifice_velocity_from_mdot` | mdot: kg/s, rho: kg/m^3, d: m, Z: - (default 1.0)                        | m/s                      |
| `orifice_area_from_beta`     | D: m, beta: -                                                            | m^2                      |
| `beta_from_diameters`        | d: m, D: m                                                               | -                        |
| `orifice_Re_d_from_mdot`     | mdot: kg/s, d: m, mu: Pa*s                                               | - (Reynolds number)      |

### geometry.h - Geometric Utilities

| Function                     | Input Units            | Output Unit  |
|------------------------------|------------------------|--------------|
| `pipe_area`                  | D: m                   | m^2          |
| `annular_area`               | D_outer: m, D_inner: m | m^2          |
| `pipe_volume`                | D: m, L: m             | m^3          |
| `pipe_roughness`             | material: str          | m            |
| `standard_pipe_roughness`    | -                      | dict[str, m] |
| `hydraulic_diameter`         | A: m^2, P_wetted: m    | m            |
| `hydraulic_diameter_rect`    | a: m, b: m             | m            |
| `hydraulic_diameter_annulus` | D_outer: m, D_inner: m | m            |

### heat_transfer.h - Heat Transfer Correlations

| Function                       | Input Units                                                                                                     | Output Unit               |
|--------------------------------|-----------------------------------------------------------------------------------------------------------------|---------------------------|
| `nusselt_dittus_boelter`       | Re: -, Pr: -, heating: bool                                                                                     | - (Nu)                    |
| `nusselt_gnielinski`           | Re: -, Pr: -, f: - (optional)                                                                                   | - (Nu)                    |
| `nusselt_sieder_tate`          | Re: -, Pr: -, mu_ratio: -                                                                                       | - (Nu)                    |
| `nusselt_petukhov`             | Re: -, Pr: -, f: - (optional)                                                                                   | - (Nu)                    |
| `htc_from_nusselt`             | Nu: -, k: W/(m*K), L: m                                                                                         | W/(m^2*K)                 |
| `nusselt_pipe`                 | State, V: m/s, D: m                                                                                             | - (Nu)                    |
| `htc_pipe`                     | State, V: m/s, D: m                                                                                             | W/(m^2*K)                 |
| `htc_pipe`                     | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, correlation: str, heating: bool, mu_ratio: -, roughness: m | tuple(W/(m^2*K), -, -)    |
| `lmtd`                         | dT1: K, dT2: K                                                                                                  | K                         |
| `lmtd_counterflow`             | T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K                                                          | K                         |
| `lmtd_parallelflow`            | T_hot_in: K, T_hot_out: K, T_cold_in: K, T_cold_out: K                                                          | K                         |
| `overall_htc`                  | h_values: W/(m^2*K), t_over_k: m^2*K/W                                                                          | W/(m^2*K)                 |
| `overall_htc_wall`             | h_inner, h_outer: W/(m^2*K), t_over_k_layers: m^2*K/W                                                           | W/(m^2*K)                 |
| `overall_htc_tube`             | h_inner, h_outer: W/(m^2*K), t_wall: m, k_wall: W/(m*K)                                                         | W/(m^2*K)                 |
| `thermal_resistance`           | h: W/(m^2*K), A: m^2                                                                                            | K/W                       |
| `thermal_resistance_wall`      | t: m, k: W/(m*K), A: m^2                                                                                        | K/W                       |
| `heat_rate`                    | U: W/(m^2*K), A: m^2, dT: K                                                                                     | W                         |
| `heat_flux`                    | U: W/(m^2*K), dT: K                                                                                             | W/m^2                     |
| `heat_transfer_area`           | Q: W, U: W/(m^2*K), dT: K                                                                                       | m^2                       |
| `heat_transfer_dT`             | Q: W, U: W/(m^2*K), A: m^2                                                                                      | K                         |
| `wall_temperature_profile`     | T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                | K (vector)                |
| `ntu`                          | U: W/(m^2*K), A: m^2, C_min: W/K                                                                                | -                         |
| `capacity_ratio`               | C_min: W/K, C_max: W/K                                                                                          | -                         |
| `effectiveness_counterflow`    | NTU: -, C_r: -                                                                                                  | -                         |
| `effectiveness_parallelflow`   | NTU: -, C_r: -                                                                                                  | -                         |
| `heat_rate_from_effectiveness` | epsilon: -, C_min: W/K, T_hot_in: K, T_cold_in: K                                                               | W                         |
| `adiabatic_wall_temperature`   | T_hot: K, T_coolant: K, eta: -                                                                                  | K                         |
| `cooled_wall_heat_flux`        | T_hot: K, T_coolant: K, h_hot: W/(m^2*K), h_coolant: W/(m^2*K), eta: -, t_wall: m, k_wall: W/(m*K)              | W/m^2                     |
| `heat_flux_from_T_at_edge`     | T_measured: K, edge_idx, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                       | W/m^2                     |
| `heat_flux_from_T_at_depth`    | T_measured: K, depth: m, T_hot: K, T_cold: K, h_hot, h_cold: W/(m^2*K), thicknesses: m, conductivities: W/(m*K) | W/m^2                     |
| `bulk_T_from_edge_T_and_q`     | T_measured: K, edge_idx, q: W/m^2, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W, solve_for: str                  | K                         |
| `dT_edge_dT_hot`               | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                           | - (dT_edge/dT_hot)        |
| `dT_edge_dT_cold`              | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                           | - (dT_edge/dT_cold)       |
| `dT_edge_dT_bulk`              | edge_idx, h_hot, h_cold: W/(m^2*K), t_over_k: m^2*K/W                                                           | - (dT/dT_hot, dT/dT_cold) |
| `dT_edge_dq`                   | edge_idx, h_hot: W/(m^2*K), t_over_k: m^2*K/W                                                                   | K*m^2/W (dT_edge/dq)      |

### acoustics.h - Acoustic Mode Analysis

| Function                             | Input Units                                                                                                                                                                                          | Output Unit               |
|--------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|
| `tube_axial_modes`                   | Tube (L: m, D: m), c: m/s, upstream: BC, downstream: BC, n_max                                                                                                                                       | Hz (vector)               |
| `annulus_axial_modes`                | Annulus (L: m, D_inner: m, D_outer: m), c: m/s, upstream: BC, downstream: BC, n_max                                                                                                                  | Hz (vector)               |
| `annulus_azimuthal_modes`            | Annulus, c: m/s, m_max                                                                                                                                                                               | Hz (vector)               |
| `annulus_modes`                      | Annulus, c: m/s, upstream: BC, downstream: BC, n_max, m_max                                                                                                                                          | Hz (vector)               |
| `modes_in_range`                     | modes, f_min: Hz, f_max: Hz                                                                                                                                                                          | Hz (vector)               |
| `closest_mode`                       | modes, f_target: Hz                                                                                                                                                                                  | AcousticMode              |
| `min_mode_separation`                | modes                                                                                                                                                                                                | Hz                        |
| `axial_mode_upstream`                | f0: Hz, M: -                                                                                                                                                                                         | Hz                        |
| `axial_mode_downstream`              | f0: Hz, M: -                                                                                                                                                                                         | Hz                        |
| `axial_mode_split`                   | f0: Hz, M: -                                                                                                                                                                                         | (Hz, Hz)                  |
| `helmholtz_frequency`                | V: m^3, A_neck: m^2, L_neck: m, c: m/s, end_correction: -                                                                                                                                            | Hz                        |
| `strouhal`                           | f: Hz, L: m, u: m/s                                                                                                                                                                                  | - (St)                    |
| `frequency_from_strouhal`            | St: -, L: m, u: m/s                                                                                                                                                                                  | Hz                        |
| `quarter_wave_frequency`             | L: m, c: m/s                                                                                                                                                                                         | Hz                        |
| `orifice_impedance_with_flow`        | freq: Hz, u_bias: m/s, u_grazing: m/s, d_orifice: m, l_orifice: m, porosity: -, Cd: -, rho: kg/m^3, c: m/s                                                                                           | complex [-]               |
| `quarter_wave_resonator_tmm`         | freq: Hz, L_tube: m, A_duct: m^2, A_tube: m^2, d_orifice: m, l_orifice: m, porosity: -, Cd: -, u_bias: m/s, u_grazing: m/s, rho: kg/m^3, c: m/s                                                      | TransferMatrix            |
| `absorption_from_impedance_norm`     | z_norm: complex [-]                                                                                                                                                                                  | - (alpha)                 |
| `liner_sdof_impedance_norm`          | freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                           | complex [-]               |
| `liner_sdof_absorption`              | freq: Hz, orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                           | - (alpha)                 |
| `sweep_liner_sdof_absorption`        | freqs: Hz (vector), orifice: LinerOrificeGeometry, cavity: LinerCavity, flow: LinerFlowState, medium: AcousticMedium                                                                                 | - (alpha vector)          |
| `liner_2dof_serial_impedance_norm`   | freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium           | complex [-]               |
| `liner_2dof_serial_absorption`       | freq: Hz, face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium           | - (alpha)                 |
| `sweep_liner_2dof_serial_absorption` | freqs: Hz (vector), face_orifice: LinerOrificeGeometry, septum_orifice: LinerOrificeGeometry, depth_1: m, depth_2: m, face_flow: LinerFlowState, septum_flow: LinerFlowState, medium: AcousticMedium | - (alpha vector)          |
| `is_whistling_risk`                  | freq: Hz, u_bias: m/s, d_orifice: m                                                                                                                                                                  | bool                      |
| `annular_duct_modes_analytical`      | geom: AnnularDuctGeometry, c: m/s, f_max: Hz, bc_ends: BoundaryCondition                                                                                                                             | AnnularMode vector        |
| `to_acoustic_geometry`               | flow_geom: CanAnnularFlowGeometry, n_cans: -, L_can: m, D_can: m                                                                                                                                     | CanAnnularGeometry struct |
| `half_wave_frequency`                | L: m, c: m/s                                                                                                                                                                                         | Hz                        |
| `stokes_layer`                       | nu: m^2/s, f: Hz                                                                                                                                                                                     | m                         |
| `thermal_layer`                      | alpha: m^2/s, f: Hz                                                                                                                                                                                  | m                         |
| `effective_viscothermal_layer`       | delta_nu: m, delta_kappa: m, gamma: -                                                                                                                                                                | m                         |
| `helmholtz_Q`                        | V: m^3, A_neck: m^2, L_neck: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz                                                                                                                             | - (Q)                     |
| `tube_Q`                             | L: m, D: m, nu: m^2/s, alpha: m^2/s, gamma: -, f: Hz                                                                                                                                                 | - (Q)                     |
| `damping_ratio`                      | Q: -                                                                                                                                                                                                 | - (zeta)                  |
| `bandwidth`                          | f0: Hz, Q: -                                                                                                                                                                                         | Hz                        |

### geometry.h - Residence Time

| Function                          | Input Units                                           | Output Unit |
|-----------------------------------|-------------------------------------------------------|-------------|
| `residence_time`                  | V: m^3, Q: m^3/s                                      | s           |
| `residence_time_tube`             | tube: Tube, Q: m^3/s                                  | s           |
| `residence_time_annulus`          | annulus: Annulus, Q: m^3/s                            | s           |
| `residence_time_can_annular`      | geom: CanAnnularFlowGeometry, Q: m^3/s                | s           |
| `residence_time_mdot`             | V: m^3, mdot: kg/s, rho: kg/m^3                       | s           |
| `residence_time_mdot_can_annular` | geom: CanAnnularFlowGeometry, mdot: kg/s, rho: kg/m^3 | s           |
| `space_velocity`                  | Q: m^3/s, V: m^3                                      | 1/s         |

### heat_transfer.h - ChannelResult fields

| Function              | Input Units | Output Unit             |
|-----------------------|-------------|-------------------------|
| `ChannelResult::h`    | -           | W/(m^2*K)               |
| `ChannelResult::Nu`   | -           | - (Nu)                  |
| `ChannelResult::Re`   | -           | - (Re)                  |
| `ChannelResult::Pr`   | -           | - (Pr)                  |
| `ChannelResult::f`    | -           | - (friction/loss coeff) |
| `ChannelResult::dP`   | -           | Pa                      |
| `ChannelResult::M`    | -           | - (M)                   |
| `ChannelResult::T_aw` | -           | K                       |
| `ChannelResult::q`    | -           | W/m^2                   |

### heat_transfer.h - Channel flow functions (HTC + pressure drop)

| Function              | Input Units                                                                                                           | Output Unit   |
|-----------------------|-----------------------------------------------------------------------------------------------------------------------|---------------|
| `channel_smooth`      | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, T_wall: K                                             | ChannelResult |
| `channel_ribbed`      | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, e_D: -, pitch_to_height: -, alpha_deg: deg, T_wall: K | ChannelResult |
| `channel_dimpled`     | T: K, P: Pa, X: mol/mol, velocity: m/s, diameter: m, length: m, d_Dh: -, h_d: -, S_d: -, T_wall: K                    | ChannelResult |
| `channel_pin_fin`     | T: K, P: Pa, X: mol/mol, velocity: m/s, channel_height: m, pin_diameter: m, S_D: -, X_D: -, N_rows: -, T_wall: K      | ChannelResult |
| `channel_impingement` | T: K, P: Pa, X: mol/mol, mdot_jet: kg/s, d_jet: m, z_D: -, x_D: -, y_D: -, A_target: m^2, T_wall: K                   | ChannelResult |

### cooling_correlations.h - Pin fin friction (new scalar)

| Function           | Input Units | Output Unit |
|--------------------|-------------|-------------|
| `pin_fin_friction` | Re_d: -     | - (f_pin)   |

### correlation_status.h - Extrapolation validity utilities

| Function              | Input Units              | Output Unit     |
|-----------------------|--------------------------|-----------------|
| `is_well_behaved`     | v: any, lo: any, hi: any | bool            |
| `set_warning_handler` | handler: callable(str)   | -               |
| `get_warning_handler` | -                        | callable(str)   |
| `suppress_warnings`   | -                        | context manager |

---

## Dimensionless Quantities

| Quantity                | Symbol | Unit     | Definition                    |
|-------------------------|--------|----------|-------------------------------|
| Mach number             | M      | -        | v / a                         |
| Reynolds number         | Re     | -        | rho*V*L / mu                  |
| Prandtl number          | Pr     | -        | mu*Cp / k                     |
| Peclet number           | Pe     | -        | V*L / alpha                   |
| Isentropic exponent     | gamma  | -        | Cp / Cv                       |
| Equivalence ratio       | phi    | -        | (F/A) / (F/A)_stoich          |
| Mixture fraction        | Z      | -        | Bilger definition             |
| Discharge coefficient   | Cd     | -        | mdot_actual / mdot_ideal      |
| Friction factor (Darcy) | f      | -        | dP / (L/D * rho*v^2/2)        |
| Pressure ratio          | -      | -        | P / P0                        |
| Diameter ratio          | beta   | -        | d / D                         |

---

## Common Reference Values

| Quantity                | Value          | Unit   | Notes                    |
|-------------------------|----------------|--------|--------------------------|
| Standard pressure       | 101325         | Pa     | 1 atm                    |
| Standard temperature    | 298.15         | K      | 25 C                     |
| Standard gravity        | 9.80665        | m/s^2  | Used for Isp             |
| Sea level air density   | ~1.225         | kg/m^3 | At 15 C, 101325 Pa       |

---

## Summary: Key Unit Conventions

1. **Temperature**: Always Kelvin (K)
2. **Pressure**: Always Pascal (Pa)
3. **Molar mass**: g/mol (historical convention; convert to kg/mol for mass-basis calculations)
4. **Thermodynamic properties**: Molar basis (J/mol, J/(mol*K)) in thermo functions
5. **Compressible flow**: Mass basis (J/kg, J/(kg*K)) in flow solutions
6. **Fractions**: mole fractions X (mol/mol), mass fractions Y (kg/kg)
7. **Relative humidity**: Fraction (0-1), not percentage
8. **Angles**: Radians (rad)
