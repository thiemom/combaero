"""Internal solver tools re-exported from _core.

These are exact Jacobian helpers and internal residual calculators used by the network solver.
They are not part of the public API and may change without notice.
"""

from . import _core

# Re-export all internal solver tools from _core
P0_from_static_and_jacobian_M = _core.P0_from_static_and_jacobian_M
T0_from_static_and_jacobian_M = _core.T0_from_static_and_jacobian_M
T_adiabatic_wall_and_jacobian_v = _core.T_adiabatic_wall_and_jacobian_v
adiabatic_T_complete_and_jacobian_T = _core.adiabatic_T_complete_and_jacobian_T
adiabatic_T_complete_and_jacobian_T_from_streams = (
    _core.adiabatic_T_complete_and_jacobian_T_from_streams
)
adiabatic_T_equilibrium_and_jacobians = _core.adiabatic_T_equilibrium_and_jacobians
adiabatic_T_equilibrium_and_jacobians_from_streams = (
    _core.adiabatic_T_equilibrium_and_jacobians_from_streams
)
area_change_residuals_and_jacobian = _core.area_change_residuals_and_jacobian
channel_compressible_mdot_and_jacobian = _core.channel_compressible_mdot_and_jacobian
channel_compressible_residuals_and_jacobian = _core.channel_compressible_residuals_and_jacobian
channel_residuals_and_jacobian = _core.channel_residuals_and_jacobian
combustor_residuals_and_jacobians = _core.combustor_residuals_and_jacobians
conical_area_change_residuals_and_jacobian = _core.conical_area_change_residuals_and_jacobian
density_and_jacobians = _core.density_and_jacobians
dimple_friction_multiplier_and_jacobian = _core.dimple_friction_multiplier_and_jacobian
dimple_nusselt_enhancement_and_jacobian = _core.dimple_nusselt_enhancement_and_jacobian
effusion_effectiveness_and_jacobian = _core.effusion_effectiveness_and_jacobian
enthalpy_and_jacobian = _core.enthalpy_and_jacobian
film_cooling_effectiveness_and_jacobian = _core.film_cooling_effectiveness_and_jacobian
friction_and_jacobian = _core.friction_and_jacobian
friction_and_jacobian_colebrook = _core.friction_and_jacobian_colebrook
friction_and_jacobian_haaland = _core.friction_and_jacobian_haaland
friction_and_jacobian_petukhov = _core.friction_and_jacobian_petukhov
friction_and_jacobian_serghides = _core.friction_and_jacobian_serghides
impingement_nusselt_and_jacobian = _core.impingement_nusselt_and_jacobian
lossless_pressure_and_jacobian = _core.lossless_pressure_and_jacobian
mach_number_and_jacobian_v = _core.mach_number_and_jacobian_v
mixer_from_streams_and_jacobians = _core.mixer_from_streams_and_jacobians
momentum_chamber_residual_and_jacobian = _core.momentum_chamber_residual_and_jacobian
nusselt_and_jacobian_dittus_boelter = _core.nusselt_and_jacobian_dittus_boelter
nusselt_and_jacobian_gnielinski = _core.nusselt_and_jacobian_gnielinski
nusselt_and_jacobian_petukhov = _core.nusselt_and_jacobian_petukhov
nusselt_and_jacobian_sieder_tate = _core.nusselt_and_jacobian_sieder_tate
orifice_compressible_mdot_and_jacobian = _core.orifice_compressible_mdot_and_jacobian
orifice_compressible_residuals_and_jacobian = _core.orifice_compressible_residuals_and_jacobian
orifice_mdot_and_jacobian = _core.orifice_mdot_and_jacobian
orifice_residuals_and_jacobian = _core.orifice_residuals_and_jacobian
pin_fin_friction_and_jacobian = _core.pin_fin_friction_and_jacobian
pin_fin_nusselt_and_jacobian = _core.pin_fin_nusselt_and_jacobian
pressure_loss_and_jacobian = _core.pressure_loss_and_jacobian
rib_enhancement_factor_high_re_and_jacobian = _core.rib_enhancement_factor_high_re_and_jacobian
viscosity_and_jacobians = _core.viscosity_and_jacobians
wall_coupling_and_jacobian = _core.wall_coupling_and_jacobian
wall_coupling_and_jacobian_multilayer = _core.wall_coupling_and_jacobian_multilayer

__all__ = [
    "P0_from_static_and_jacobian_M",
    "T0_from_static_and_jacobian_M",
    "T_adiabatic_wall_and_jacobian_v",
    "adiabatic_T_complete_and_jacobian_T",
    "adiabatic_T_complete_and_jacobian_T_from_streams",
    "adiabatic_T_equilibrium_and_jacobians",
    "adiabatic_T_equilibrium_and_jacobians_from_streams",
    "area_change_residuals_and_jacobian",
    "channel_compressible_mdot_and_jacobian",
    "channel_compressible_residuals_and_jacobian",
    "channel_residuals_and_jacobian",
    "combustor_residuals_and_jacobians",
    "conical_area_change_residuals_and_jacobian",
    "density_and_jacobians",
    "dimple_friction_multiplier_and_jacobian",
    "dimple_nusselt_enhancement_and_jacobian",
    "effusion_effectiveness_and_jacobian",
    "enthalpy_and_jacobian",
    "film_cooling_effectiveness_and_jacobian",
    "friction_and_jacobian",
    "friction_and_jacobian_colebrook",
    "friction_and_jacobian_haaland",
    "friction_and_jacobian_petukhov",
    "friction_and_jacobian_serghides",
    "impingement_nusselt_and_jacobian",
    "lossless_pressure_and_jacobian",
    "mach_number_and_jacobian_v",
    "mixer_from_streams_and_jacobians",
    "momentum_chamber_residual_and_jacobian",
    "nusselt_and_jacobian_dittus_boelter",
    "nusselt_and_jacobian_gnielinski",
    "nusselt_and_jacobian_petukhov",
    "nusselt_and_jacobian_sieder_tate",
    "orifice_compressible_mdot_and_jacobian",
    "orifice_compressible_residuals_and_jacobian",
    "orifice_mdot_and_jacobian",
    "orifice_residuals_and_jacobian",
    "pin_fin_friction_and_jacobian",
    "pin_fin_nusselt_and_jacobian",
    "pressure_loss_and_jacobian",
    "rib_enhancement_factor_high_re_and_jacobian",
    "viscosity_and_jacobians",
    "wall_coupling_and_jacobian",
    "wall_coupling_and_jacobian_multilayer",
]
