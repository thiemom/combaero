import combaero as cb
from combaero.network import components


def test_universal_pressure_block():
    """Verify that _element_pressure_block emits the required keys including Tt and Pt."""
    # Create simple mock states
    st_in = components.MixtureState(
        T=300.0,
        P=100000.0,
        Pt=100500.0,
        Tt=301.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    st_out = components.MixtureState(
        T=290.0,
        P=95000.0,
        Pt=96000.0,
        Tt=295.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    pb = components._element_pressure_block(st_in, st_out, mach_in=0.1, mach_out=0.2)
    expected_keys = {
        "P_in",
        "P_out",
        "T_in",
        "T_out",
        "Pt_in",
        "Pt_out",
        "Tt_in",
        "Tt_out",
        "mach_in",
        "mach_out",
        "dPt",
        "dP",
        "pr_total",
        "pr_static",
    }
    assert set(pb.keys()) == expected_keys

    assert pb["Pt_in"] == 100500.0
    assert pb["Tt_in"] == 301.0
    assert pb["dP"] == 5000.0
    assert pb["dPt"] == 4500.0


def test_orifice_element_diagnostics():
    """Verify OrificeElement diagnostics schema including mach_in and mach_out."""
    elem = components.OrificeElement("ori", "in", "out", area=0.1)

    st_in = components.MixtureState(
        T=300.0,
        P=100000.0,
        Pt=100500.0,
        Tt=300.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    st_out = components.MixtureState(
        T=290.0,
        P=95000.0,
        Pt=96000.0,
        Tt=290.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    diag = elem.diagnostics(st_in, st_out)
    assert "ref_location" in diag
    assert diag["ref_location"] == "inlet"
    assert "mach_in" in diag
    assert "mach_out" in diag
    assert "mach_throat" in diag
    assert "Cd" in diag
    assert "p_ratio" not in diag  # Removed in Phase 2
    assert "pr_total" in diag
    assert "pr_static" in diag


def test_area_change_element_diagnostics():
    """Verify AreaChangeElement diagnostics schema."""
    elem = components.AreaChangeElement("ac", "in", "out", F0=0.1, F1=0.05)

    st_in = components.MixtureState(
        T=300.0,
        P=100000.0,
        Pt=100500.0,
        Tt=300.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    st_out = components.MixtureState(
        T=290.0,
        P=95000.0,
        Pt=96000.0,
        Tt=290.0,
        m_dot=1.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    diag = elem.diagnostics(st_in, st_out)
    assert "ref_location" in diag
    assert diag["ref_location"] == "small_section"
    assert "mach_loss_ref" in diag
    assert "area_ratio" in diag
    assert "dP_static" in diag
    assert "zeta" in diag
    assert "mach_small" not in diag  # Renamed to mach_loss_ref
    assert "ratio" not in diag  # Renamed to area_ratio


def test_plenum_diagnostics():
    """Verify PlenumNode diagnostics emit explicit Tt and Pt."""
    node = components.PlenumNode("plenum")
    st = components.MixtureState(
        T=300.0,
        P=100000.0,
        Pt=100000.0,
        Tt=300.0,
        m_dot=0.0,
        Y=list(cb.mole_to_mass(cb.species.dry_air())),
    )

    diag = node.diagnostics(st)
    assert "Tt" in diag
    assert "Pt" in diag
    assert diag["Tt"] == 300.0
    assert diag["Pt"] == 100000.0
    assert "rho" in diag
