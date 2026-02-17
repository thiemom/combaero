"""
Combaero Acoustics Examples - Industrial Damping & Dispersion Analysis

This script demonstrates acoustic analysis for gas turbine combustors,
focusing on engineering metrics like Absorption Coefficient (alpha)
and Transmission Loss (TL).
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np

import combaero as cb

warnings.filterwarnings("ignore")


def print_header(title):
    print(f"\n{'=' * 70}\n {title}\n{'=' * 70}")


# -------------------------------------------------------------------------
# DEMO 1: Aero-Engine (Pure Annular)
# -------------------------------------------------------------------------


def demo_aero_annular():
    """Compare argument-principle annular solver against analytical reference."""
    print_header("DEMO 1: Aero-Engine Annular Combustor (Dispersion)")

    geom = cb.AnnularDuctGeometry()
    geom.length = 0.35  #
    geom.radius_inner = 0.40  #
    geom.radius_outer = 0.48  #
    geom.n_azimuthal_max = 8  #

    props = cb.air_properties(T=1700.0, P=25e5)
    c, rho = props.a, props.rho  #

    modes_ap = cb.annular_duct_eigenmodes(
        geom, c, rho, f_max=4000.0, bc_ends=cb.BoundaryCondition.Closed
    )
    modes_ref = cb.annular_duct_modes_analytical(
        geom, c, f_max=4000.0, bc_ends=cb.BoundaryCondition.Closed
    )

    if not modes_ap:
        print("Warning: No modes found. Check geometry/frequency range.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    modes_azi = [m for m in modes_ap if m.n_axial == 0]
    modes_mix = [m for m in modes_ap if m.n_axial > 0]

    ax.plot(
        [m.m_azimuthal for m in modes_azi],
        [m.frequency for m in modes_azi],
        "o-",
        color="#1f77b4",
        label="Azimuthal Modes (n=0)",
        markersize=8,
    )

    if modes_mix:
        sc = ax.scatter(
            [m.m_azimuthal for m in modes_mix],
            [m.frequency for m in modes_mix],
            c=[m.n_axial for m in modes_mix],
            cmap="autumn_r",
            label="Mixed Modes (n>0)",
            zorder=3,
            edgecolors="k",
        )
        plt.colorbar(sc, ax=ax, label="Longitudinal Order (n)")

    # Overlay analytical reference for first few modes
    ref_subset = modes_ref[: min(40, len(modes_ref))]
    ax.scatter(
        [m.m_azimuthal for m in ref_subset],
        [m.frequency for m in ref_subset],
        marker="x",
        color="black",
        alpha=0.5,
        label="Analytical reference",
    )

    ax.set_xlabel("Azimuthal Mode Order (m)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("Dispersion Relation: Annular Combustor")
    ax.legend()
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------------------------
# DEMO 2: Heavy-Duty GT (Can-Annular)
# -------------------------------------------------------------------------


def demo_industrial_can_annular():
    """Restored Bloch-Floquet passband shading."""
    print_header("DEMO 2: Industrial GT Can-Annular System (Bloch Analysis)")

    geom = cb.CanAnnularGeometry()
    geom.n_cans = 16  #
    geom.length_can = 0.7  #
    geom.area_can = 0.08  #
    geom.radius_plenum = 1.8  #
    geom.area_plenum = 0.2  #

    # Solves for modes in a periodically coupled system
    modes = cb.can_annular_eigenmodes(geom, 850.0, 550.0, 2.0, 5.0, f_max=1200.0)

    fig, ax = plt.subplots(figsize=(10, 6))
    f_vals = [m.frequency for m in modes]
    m_vals = [m.m_azimuthal for m in modes]

    # KDE Shading for Passbands
    y_range = np.linspace(0, 1200, 300)
    density = np.zeros_like(y_range)
    for f in f_vals:
        density += np.exp(-0.5 * ((y_range - f) / 10.0) ** 2)

    ax.fill_betweenx(
        y_range,
        -0.5,
        geom.n_cans / 2 + 0.5,
        where=density > 0.1,
        color="blue",
        alpha=0.1,
        label="Passbands",
    )

    ax.scatter(m_vals, f_vals, c="gray", edgecolors="k", zorder=5)
    ax.set_title(f"Bloch Spectrum: {geom.n_cans}-Can Engine")
    ax.set_xlabel("Azimuthal Index (m)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_xlim(-0.5, geom.n_cans / 2 + 0.5)
    plt.tight_layout()
    plt.show()


# -------------------------------------------------------------------------
# DEMO 3: Acoustic Damping (Liner Design)
# -------------------------------------------------------------------------


def demo_damping_liner():
    """Liner performance with bias and grazing flow using streamlined API."""
    print_header("DEMO 3: Acoustic Liner Damping Design")
    c, rho = 340.0, 1.2
    bias_velocities = [0.0, 10.0, 20.0, 30.0]
    freqs = np.linspace(100.0, 2000.0, 500)

    medium = cb.AcousticMedium()
    medium.rho = rho
    medium.c = c

    orifice = cb.LinerOrificeGeometry()
    orifice.d_orifice = 0.003
    orifice.l_orifice = 0.005
    orifice.porosity = 0.08
    orifice.Cd = 0.7

    cavity = cb.LinerCavity()
    cavity.depth = 0.05

    fig, ax = plt.subplots(figsize=(10, 6))
    for u_bias in bias_velocities:
        flow = cb.LinerFlowState()
        flow.u_bias = u_bias
        flow.u_grazing = 80.0
        alphas = cb.sweep_liner_sdof_absorption(freqs.tolist(), orifice, cavity, flow, medium)
        ax.plot(freqs, alphas, label=f"Bias = {u_bias} m/s")

    ax.set_title("Liner Damping vs Bias Flow")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Absorption Coefficient (alpha)")
    ax.legend()
    plt.show()


# -------------------------------------------------------------------------
# DEMO 4: Serial Dampers (SDOF vs 2-DOF)
# -------------------------------------------------------------------------


def demo_serial_dampers():
    """Two-Degree-of-Freedom serial liner using streamlined API."""
    print_header("DEMO 4: Serial Dampers (SDOF vs 2-DOF)")
    total_depth = 0.05
    u_bias, u_grazing, c, rho = 10.0, 80.0, 340.0, 1.2
    freqs = np.linspace(100.0, 2000.0, 1000)

    medium = cb.AcousticMedium()
    medium.rho = rho
    medium.c = c

    face_orifice = cb.LinerOrificeGeometry()
    face_orifice.d_orifice = 0.003
    face_orifice.l_orifice = 0.005
    face_orifice.porosity = 0.08
    face_orifice.Cd = 0.7

    septum_orifice = cb.LinerOrificeGeometry()
    septum_orifice.d_orifice = 0.003
    septum_orifice.l_orifice = 0.005
    septum_orifice.porosity = 0.08
    septum_orifice.Cd = 0.7

    cavity = cb.LinerCavity()
    cavity.depth = total_depth

    face_flow = cb.LinerFlowState()
    face_flow.u_bias = u_bias
    face_flow.u_grazing = u_grazing

    septum_flow = cb.LinerFlowState()
    septum_flow.u_bias = u_bias
    septum_flow.u_grazing = 0.0

    alpha_sdof = cb.sweep_liner_sdof_absorption(
        freqs.tolist(), face_orifice, cavity, face_flow, medium
    )
    alpha_2dof = cb.sweep_liner_2dof_serial_absorption(
        freqs.tolist(),
        face_orifice,
        septum_orifice,
        total_depth / 2.0,
        total_depth / 2.0,
        face_flow,
        septum_flow,
        medium,
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(freqs, alpha_sdof, "k--", label="SDOF (Single Vol)", alpha=0.6)
    ax.plot(freqs, alpha_2dof, "r-", label="2-DOF (Serial)", linewidth=2)
    ax.fill_between(
        freqs,
        alpha_2dof,
        alpha_sdof,
        where=(np.array(alpha_2dof) > np.array(alpha_sdof)),
        color="red",
        alpha=0.1,
        label="Bandwidth Gain",
    )
    ax.set_title("Bandwidth Comparison: SDOF vs 2-DOF")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Absorption Coefficient (alpha)")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    demo_aero_annular()
    demo_industrial_can_annular()
    demo_damping_liner()
    demo_serial_dampers()
