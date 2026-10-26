"""Unit catalogue for result-export column headers.

The SI unit of every result key, for the CSV produced by ``/export``
(``P [Pa]``, ``T [K]`` ...) and, through ``unit_of`` and
``scripts/gui_result_keys.py``, for the GUI: the frontend's
``QUANTITY_CATALOGUE`` display units are checked against these
(``quantities.test.ts``).

Keep this module dependency-free; it is imported by ``main.py`` during the
CSV-streaming response.
"""

import re

# ---------------------------------------------------------------------------
# Scalar quantity → unit mapping
# ---------------------------------------------------------------------------

UNIT_MAP: dict[str, str] = {
    # --- Thermodynamic state ---
    "T": "K",
    "T_total": "K",
    "T_in": "K",
    "T_out": "K",
    "Tt_in": "K",
    "Tt_out": "K",
    "T_aw": "K",
    "T_hot": "K",
    "T_interface": "K",
    "P": "Pa",
    "P_total": "Pa",
    "P_in": "Pa",
    "P_out": "Pa",
    "Pt_in": "Pa",
    "Pt_out": "Pa",
    "dP": "Pa",
    "dPt": "Pa",
    "pr_total": "-",
    "pr_static": "-",
    # --- Flow ---
    "m_dot": "kg/s",
    "velocity": "m/s",
    "v": "m/s",
    "a": "m/s",
    "mach": "-",
    "mach_in": "-",
    "mach_out": "-",
    "mach_throat": "-",
    "mach_loss_ref": "-",
    "dPt_loss": "Pa",  # AreaChangeElement: the closure's stagnation loss (#731)
    "choked": "-",
    # --- Node stagnation ---
    "Tt": "K",
    "Pt": "Pa",
    # --- Transport ---
    "mu": "Pa·s",
    "k": "W/m/K",
    "nu": "m²/s",
    "Pr": "-",
    "Re": "-",
    "Dh": "m",
    # --- Thermo per-mass (node/element state) ---
    "rho": "kg/m³",
    "h": "J/kg",
    "s": "J/kg/K",
    "u": "J/kg",
    "cp": "J/kg/K",
    "cv": "J/kg/K",
    "gamma": "-",
    "mw": "g/mol",  # State::mw (include/units_data.h)
    # --- Thermo reference-state block (Phase 2, element inlet state) ---
    "rho_ref": "kg/m³",
    "h_ref": "J/kg",
    "s_ref": "J/kg/K",
    "u_ref": "J/kg",
    "cp_ref": "J/kg/K",
    "cv_ref": "J/kg/K",
    "gamma_ref": "-",
    "a_ref": "m/s",
    "mw_ref": "g/mol",
    "mu_ref": "Pa·s",
    "k_ref": "W/m/K",
    "nu_ref": "m²/s",
    "Pr_ref": "-",
    # --- Heat transfer ---
    "Nu": "-",
    "htc": "W/m²/K",
    "h_a": "W/m²/K",
    "h_b": "W/m²/K",
    "Q": "W",
    "Q_to_boundary": "W",  # wall heat that left with a stream into a boundary
    "Q_wall_out": "W",
    "Q_withheld": "W",
    "Q_fraction": "W",
    # --- Tee junction ---
    "m_dot_com": "kg/s",
    "m_dot_straight": "kg/s",
    "m_dot_branch": "kg/s",
    "mass_flow_ratio": "kg/kg",
    "K_straight": "-",
    "K_branch": "-",
    "n_reversed_ports": "-",  # ideal momentum tee: ports run against the wiring
    # --- Element diagnostics ---
    "xi": "-",
    "theta": "-",
    "zeta": "-",
    "ratio": "-",  # AreaChangeElement (renamed to area_ratio in Phase 2)
    "area_ratio": "-",  # Phase 2 name
    "is_correlation": "-",
    "Cd_in_range": "-",  # 1 inside the discharge correlation's source range
    # --- Impingement plate / crossflow (#465) ---
    "n_holes": "-",
    "Re_j": "-",
    "Gc_Gj": "-",
    "Gc_Gj_closed_form": "-",
    "surface_extrapolated": "-",
    "dP_momentum": "Pa",
    # --- Effusion plate wall (#471) ---
    "h_internal_plate_area": "W/m²/K",
    "h_gas_unblown": "W/m²/K",
    "h_gas": "W/m²/K",
    "gas_augmentation": "-",
    "T_gas": "K",
    "eta_film": "-",
    "eta_overall": "-",
    "T_adiabatic_wall": "K",
    "T_wall_hot": "K",
    "T_wall_cold": "K",
    "q_wall": "W/m²",
    "Q_wall": "W",
    "blowing_ratio": "-",
    "velocity_ratio": "-",
    "density_ratio": "-",
    "U_gas": "m/s",
    "gas_heat_flux_imposed": "W/m²",
    # --- Effusion liner (#471) ---
    "m_dot_out": "kg/s",
    "m_dot_bleed": "kg/s",
    "bleed_fraction": "-",
    "n_segments": "-",
    "U1_over_Vi": "-",
    "U1_over_Vi_max": "-",
    "bleed_nonuniformity": "-",
    "station_psi_min": "-",
    "T_wall_hot_max": "K",
    "n_rows": "-",
    "Re_j_min": "-",
    "Re_j_max": "-",
    "Gc_Gj_max": "-",
    "Nu_mean": "-",
    "htc_mean": "W/m²/K",
    "jet_flow_nonuniformity": "-",
    # --- Impingement array bypass (#467) ---
    "m_dot_bypass": "kg/s",
    "bypass_to_jet_flow": "-",
    "bypass_Tt_minus_supply": "K",
    "heat_transfer_outside_source": "-",
    "bypass_reversed": "-",
    "f": "-",
    "Cd": "-",
    "phi": "-",
    # --- Geometry ---
    "area": "m²",
    "diameter": "m",
    "length": "m",
    "roughness": "m",
}

# Units of the result keys the export does NOT label (yet). Giving one of them
# a unit in UNIT_MAP renames its CSV / DataFrame column ("loss_ratio" becomes
# "loss_ratio [-]"), which breaks a notebook that reads the column by name, so they
# are kept apart: the GUI shows these units (resultKeys.generated.json, via
# unit_of) and the export keeps its columns. Moving them into UNIT_MAP is the
# user's call. test_gui_result_keys_sync.py fails on an emitted key with no
# unit in either map. Flags (0/1) and counts are "-".
GUI_ONLY_UNITS: dict[str, str] = {
    # Solver views
    "sigma_norm": "-",
    "sigma_norm_high": "-",
    "sigma_norm_low": "-",
    "incompressible_mach": "-",
    "incompressible_mach_limit": "-",
    # Node states
    "T_static": "K",
    # Channel / Fanno-Rayleigh duct (#551, #636)
    "M_in_duct": "-",
    "M_exit_duct": "-",
    "f_surface": "-",
    "Re_surface": "-",
    "heat_area_ratio": "-",
    "heat_area_exceeds_wetted": "-",
    "rayleigh_in_momentum": "-",
    "rayleigh_omitted": "-",
    "rayleigh_omitted_frac": "-",
    # Area change (#731)
    "F0": "m\u00b2",
    "F1": "m\u00b2",
    "Re_small": "-",
    "expansion": "-",
    # Orifice approach and ISO 5167 taps (#649/#645/#651)
    "beta": "-",
    "beta_flow": "-",
    "E": "-",
    "Cd_pipe_ratio": "-",
    "P_approach": "Pa",
    "dP_tap": "Pa",
    "dP_permanent": "Pa",
    "loss_ratio": "-",
    "pressure_recovery": "-",
    "recovery_in_range": "-",
    # Discharge holes, film and slot (#386, #714)
    "G_coolant": "kg/m\u00b2/s",
    "G_jet": "kg/m\u00b2/s",
    "G_slot": "kg/m\u00b2/s",
    "dP_drive": "Pa",
    "L_over_d": "-",
    "hole_length": "m",
    "hydraulic_diameter": "m",
    "slot_height": "m",
    # A fraction on every element (#773); the slot film reports it.
    "turbulence_intensity": "-",
    "is_ingesting": "-",
    # Effusion plate (#471)
    "Nu_approach": "-",
    "Nu_internal": "-",
    "Nu_throat": "-",
    "Re_hole": "-",
    "h_hole_area": "W/m\u00b2/K",
    "h_plate_area": "W/m\u00b2/K",
    "h_gas_unblown_Re": "-",
    "h_gas_unblown_Pr": "-",
    "h_gas_unblown_Nu": "-",
    "T_gas_from_main_inlet": "-",
    "film_extrapolated": "-",
    "gas_side_chamber": "-",
    "area_approach_total": "m\u00b2",
    "area_hole_total": "m\u00b2",
    "approach_to_hole_area_ratio": "-",  # approach area / hole area (#773)
    "coolant_crossflow_ignored": "-",
    "coolant_ht_plenum_assumed": "-",
    "hole_count_exact": "-",
    "internal_Nu_multiplier": "-",
    "pitch_over_d": "-",
    "porosity": "-",
    "wall_conductivity": "W/m/K",
    # Effusion liner (#471)
    "crossflow_cd_beyond_8pct": "-",
    "crossflow_cd_degraded": "-",
    "station_psi_beyond_bassett": "-",
    # Impingement array (#465)
    "n_rows_unresolved": "-",
    # Tee and chamber junctions (#517, #536)
    "K11": "-",
    "K12": "-",
    "P_jct": "Pa",
    "Pt_jct": "Pa",
    "n_ports": "-",
    "port_area_mismatch": "-",
    # Ejector (#270, #619)
    "m_dot_primary": "kg/s",
    "m_dot_secondary": "kg/s",
    "m_dot_outlet": "kg/s",
    "entrainment_ratio": "-",  # m_dot_secondary / m_dot_primary (#773)
    "entrainment_ratio_critical": "-",
    "critical_mode": "-",
    "s_choke": "-",
    "r_gas": "J/kg/K",
    "P_c_star": "Pa",
    "P_b0": "Pa",
    # Incompressible ejector: Sanger's jet pump (#588)
    "regime_incompressible": "-",
    "mach_jet": "-",
    "head_ratio_N": "-",
    "efficiency": "-",
    # Vortex (#661)
    "shaft_speed_rpm": "rpm",
    "v_theta_in": "m/s",
    "v_theta_out": "m/s",
    "M_theta_in": "-",
    "M_theta_out": "-",
    "P_swirl_in": "Pa",
    "P_swirl_out": "Pa",
    "dP_swirl_static": "Pa",
    # Thermal wall exchange (#627)
    "T_aw_a": "K",
    "T_aw_b": "K",
    "UA": "W/K",
    "NTU": "-",
    "effectiveness": "-",
    "counterflow": "-",
}

# Per-port keys of the junctions, ``port_<i>_<quantity>`` (i = the element's
# port index) and the tee's hidden port steps, ``step_<port>_<quantity>``.
_PORT_UNITS: dict[str, str] = {
    "P": "Pa",
    "T": "K",
    "m_dot": "kg/s",
    "area": "m\u00b2",
    "outer_area": "m\u00b2",
    "connected_area_ratio": "-",  # connected element area / port area (#773)
    "sign": "-",
    "K": "-",
    "flow_fraction": "-",
}
_STEP_UNITS: dict[str, str] = {"zeta": "-", "dPt": "Pa", "connected_area_ratio": "-"}
_PORT_RE = re.compile(r"^port_\d+_(\w+)$")
_STEP_RE = re.compile(r"^step_(?:common|straight|branch)_(\w+)$")


def unit_of(key: str) -> str | None:
    """The SI unit of a result key, or None when none is known (the GUI's view;
    the export labels UNIT_MAP keys only, see GUI_ONLY_UNITS)."""
    if key in UNIT_MAP:
        return UNIT_MAP[key]
    if key in GUI_ONLY_UNITS:
        return GUI_ONLY_UNITS[key]
    m = _PORT_RE.match(key)
    if m and m.group(1) in _PORT_UNITS:
        return _PORT_UNITS[m.group(1)]
    m = _STEP_RE.match(key)
    if m and m.group(1) in _STEP_UNITS:
        return _STEP_UNITS[m.group(1)]
    return None


# Composition keys render as ``Y[<label>] [-]``. ``<label>`` may be a species
# name (e.g. ``N2``, ``H2O``, ``CH4``) or a numeric index for forward
# compatibility with older exports.
_COMPOSITION_RE = re.compile(r"^[YX]\[[A-Za-z0-9_+\-()]+\]$")


def label_with_unit(column: str) -> str:
    """Return ``'<column> [<unit>]'`` when a unit is known, else ``column``.

    Meta columns such as ``id``, ``type``, ``label``, ``kind``, ``success``
    and ``is_boundary`` are passed through unchanged.
    """
    if column in _META_COLUMNS:
        return column
    if column in UNIT_MAP:
        return f"{column} [{UNIT_MAP[column]}]"
    if _COMPOSITION_RE.match(column):
        return f"{column} [-]"
    return column


_META_COLUMNS: frozenset[str] = frozenset(
    {
        "id",
        "label",
        "type",
        "kind",
        "success",
        "is_boundary",
        "combustion_method",
        "ref_location",
    }
)
