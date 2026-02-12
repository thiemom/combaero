"""Tests for geometry utility functions.

Verifies Python bindings for hydraulic diameter and residence time calculations.
Units: Dh in [m], τ in [s], SV in [1/s], A in [m2], V in [m3].
"""

import numpy as np
import pytest

import combaero as cb


class TestSimpleGeometryHelpers:
    """Test simple geometry helper functions (Phase 1)."""

    def test_pipe_area_basic(self):
        """Test pipe area calculation."""
        D = 0.1  # Diameter [m]
        A = cb.pipe_area(D)

        # A = π * (D/2)2
        expected = np.pi * (D / 2) ** 2
        assert abs(A - expected) < 1e-15

        # Units check: should be in m2
        assert A > 0
        assert isinstance(A, float)

    def test_pipe_area_matches_manual(self):
        """Test pipe area matches manual calculation exactly."""
        D = 0.05  # 50 mm diameter
        A = cb.pipe_area(D)

        # Manual calculation
        A_manual = np.pi * 0.025**2
        assert abs(A - A_manual) < 1e-15  # Machine precision

    def test_pipe_area_various_sizes(self):
        """Test pipe area for various diameters."""
        for D in [0.01, 0.05, 0.1, 0.5, 1.0]:
            A = cb.pipe_area(D)
            expected = np.pi * (D / 2) ** 2
            assert abs(A - expected) / expected < 1e-14

    def test_pipe_area_error_on_zero(self):
        """Test pipe area raises error for zero diameter."""
        with pytest.raises(ValueError):  # C++ invalid_argument → ValueError
            cb.pipe_area(0.0)

    def test_pipe_area_error_on_negative(self):
        """Test pipe area raises error for negative diameter."""
        with pytest.raises(ValueError):
            cb.pipe_area(-0.1)

    def test_annular_area_basic(self):
        """Test annular area calculation."""
        D_outer = 0.15  # Outer diameter [m]
        D_inner = 0.10  # Inner diameter [m]
        A = cb.annular_area(D_outer, D_inner)

        # A = π * ((D_outer/2)2 - (D_inner/2)2)
        expected = np.pi * ((D_outer / 2) ** 2 - (D_inner / 2) ** 2)
        assert abs(A - expected) < 1e-15

        # Units check: should be in m2
        assert A > 0
        assert isinstance(A, float)

    def test_annular_area_matches_manual(self):
        """Test annular area matches manual calculation exactly."""
        D_outer = 0.6  # 600 mm
        D_inner = 0.4  # 400 mm
        A = cb.annular_area(D_outer, D_inner)

        # Manual calculation
        A_manual = np.pi * (0.3**2 - 0.2**2)
        assert abs(A - A_manual) < 1e-15  # Machine precision

    def test_annular_area_thin_gap(self):
        """Test annular area for thin annular gap."""
        D_outer = 0.101
        D_inner = 0.100
        A = cb.annular_area(D_outer, D_inner)

        expected = np.pi * ((D_outer / 2) ** 2 - (D_inner / 2) ** 2)
        assert abs(A - expected) / expected < 1e-14

    def test_annular_area_error_on_invalid(self):
        """Test annular area raises error for D_outer <= D_inner."""
        with pytest.raises(ValueError):
            cb.annular_area(0.1, 0.1)  # Equal

        with pytest.raises(ValueError):
            cb.annular_area(0.1, 0.15)  # Inner > outer

    def test_annular_area_error_on_negative(self):
        """Test annular area raises error for negative diameter."""
        with pytest.raises(ValueError):
            cb.annular_area(0.15, -0.1)

    def test_pipe_volume_basic(self):
        """Test pipe volume calculation."""
        D = 0.1  # Diameter [m]
        L = 2.0  # Length [m]
        V = cb.pipe_volume(D, L)

        # V = π * (D/2)2 * L
        expected = np.pi * (D / 2) ** 2 * L
        assert abs(V - expected) < 1e-15

        # Units check: should be in m3
        assert V > 0
        assert isinstance(V, float)

    def test_pipe_volume_matches_manual(self):
        """Test pipe volume matches manual calculation exactly."""
        D = 0.05  # 50 mm diameter
        L = 10.0  # 10 m length
        V = cb.pipe_volume(D, L)

        # Manual calculation
        V_manual = np.pi * 0.025**2 * 10.0
        assert abs(V - V_manual) < 1e-15  # Machine precision

    def test_pipe_volume_various_sizes(self):
        """Test pipe volume for various dimensions."""
        test_cases = [(0.01, 1.0), (0.05, 5.0), (0.1, 10.0), (0.2, 2.0)]
        for D, L in test_cases:
            V = cb.pipe_volume(D, L)
            expected = np.pi * (D / 2) ** 2 * L
            assert abs(V - expected) / expected < 1e-14

    def test_pipe_volume_consistency_with_area(self):
        """Test pipe volume is consistent with pipe area."""
        D = 0.1
        L = 5.0

        A = cb.pipe_area(D)
        V = cb.pipe_volume(D, L)

        # V should equal A * L
        assert abs(V - A * L) < 1e-15  # Machine precision

    def test_pipe_volume_error_on_zero(self):
        """Test pipe volume raises error for zero dimensions."""
        with pytest.raises(ValueError):
            cb.pipe_volume(0.0, 1.0)

        with pytest.raises(ValueError):
            cb.pipe_volume(0.1, 0.0)

    def test_pipe_volume_error_on_negative(self):
        """Test pipe volume raises error for negative dimensions."""
        with pytest.raises(ValueError):
            cb.pipe_volume(-0.1, 1.0)

        with pytest.raises(ValueError):
            cb.pipe_volume(0.1, -1.0)


class TestHydraulicDiameter:
    """Test hydraulic diameter calculations."""

    def test_hydraulic_diameter_circular_pipe(self):
        """Test hydraulic diameter for circular pipe."""
        # For circular pipe: Dh = D
        D = 0.1  # Diameter [m]
        A = 3.14159 * (D / 2) ** 2  # Area [m2]
        P = 3.14159 * D  # Perimeter [m]

        Dh = cb.hydraulic_diameter(A, P)

        # Should equal diameter for circular pipe
        assert abs(Dh - D) < 0.001

        # Units check: should be in meters
        assert Dh > 0
        assert isinstance(Dh, float)

    def test_hydraulic_diameter_square_duct(self):
        """Test hydraulic diameter for square duct."""
        # For square duct: Dh = side length
        a = 0.1  # Side length [m]
        A = a * a  # Area [m2]
        P = 4 * a  # Perimeter [m]

        Dh = cb.hydraulic_diameter(A, P)

        # Dh = 4*A/P = 4*a2/(4*a) = a
        assert abs(Dh - a) < 1e-10
        assert Dh > 0

    def test_hydraulic_diameter_rect_square(self):
        """Test rectangular duct formula for square."""
        a = 0.1  # Side length [m]
        b = 0.1  # Same side length [m]

        Dh = cb.hydraulic_diameter_rect(a, b)

        # For square: Dh = 2*a*a/(a+a) = a
        assert abs(Dh - a) < 1e-10
        assert Dh > 0

    def test_hydraulic_diameter_rect_general(self):
        """Test rectangular duct with different sides."""
        a = 0.1  # Width [m]
        b = 0.2  # Height [m]

        Dh = cb.hydraulic_diameter_rect(a, b)

        # Dh = 2*a*b/(a+b) = 2*0.1*0.2/(0.1+0.2) = 0.04/0.3 = 0.1333...
        expected = 2 * a * b / (a + b)
        assert abs(Dh - expected) < 1e-10

        # Should be between min and max side
        assert min(a, b) < Dh < max(a, b)

    def test_hydraulic_diameter_annulus(self):
        """Test annular duct hydraulic diameter."""
        D_outer = 0.1  # Outer diameter [m]
        D_inner = 0.08  # Inner diameter [m]

        Dh = cb.hydraulic_diameter_annulus(D_outer, D_inner)

        # For annulus: Dh = D_outer - D_inner
        assert abs(Dh - (D_outer - D_inner)) < 1e-10
        assert abs(Dh - 0.02) < 1e-10

    def test_hydraulic_diameter_annulus_thin(self):
        """Test thin annular gap."""
        D_outer = 0.101  # Outer diameter [m]
        D_inner = 0.100  # Inner diameter [m]

        Dh = cb.hydraulic_diameter_annulus(D_outer, D_inner)

        # Very thin gap
        assert abs(Dh - 0.001) < 1e-10
        assert Dh > 0

    def test_hydraulic_diameter_units_meters(self):
        """Verify hydraulic diameter returns meters."""
        A = 0.01  # Area [m2]
        P = 0.4  # Perimeter [m]

        Dh = cb.hydraulic_diameter(A, P)

        # Dh = 4*0.01/0.4 = 0.1 m
        assert abs(Dh - 0.1) < 1e-10

        # Typical range for engineering applications
        assert 0.001 < Dh < 10.0


class TestResidenceTime:
    """Test residence time calculations."""

    def test_residence_time_basic(self):
        """Test basic residence time calculation."""
        V = 1.0  # Volume [m3]
        Q = 0.1  # Volumetric flow rate [m3/s]

        tau = cb.residence_time(V, Q)

        # τ = V/Q = 1.0/0.1 = 10 s
        assert abs(tau - 10.0) < 1e-10

        # Units check: should be in seconds
        assert tau > 0
        assert isinstance(tau, float)

    def test_residence_time_fast_flow(self):
        """Test residence time with fast flow."""
        V = 0.5  # Volume [m3]
        Q = 10.0  # High flow rate [m3/s]

        tau = cb.residence_time(V, Q)

        # τ = 0.5/10 = 0.05 s
        assert abs(tau - 0.05) < 1e-10
        assert tau > 0

    def test_residence_time_mdot_basic(self):
        """Test residence time from mass flow rate."""
        V = 1.0  # Volume [m3]
        mdot = 1.2  # Mass flow rate [kg/s]
        rho = 1.2  # Density (air) [kg/m3]

        tau = cb.residence_time_mdot(V, mdot, rho)

        # τ = V*ρ/ṁ = 1.0*1.2/1.2 = 1.0 s
        assert abs(tau - 1.0) < 1e-10
        assert tau > 0

    def test_residence_time_mdot_vs_volumetric(self):
        """Test consistency between mass and volumetric methods."""
        V = 2.0  # Volume [m3]
        rho = 1.2  # Density [kg/m3]
        Q = 0.5  # Volumetric flow rate [m3/s]
        mdot = rho * Q  # Mass flow rate [kg/s]

        tau_vol = cb.residence_time(V, Q)
        tau_mass = cb.residence_time_mdot(V, mdot, rho)

        # Should give same result
        assert abs(tau_vol - tau_mass) < 1e-10

    def test_residence_time_combustor(self):
        """Test realistic combustor residence time."""
        V = 0.1  # Combustor volume [m3]
        mdot = 0.5  # Air mass flow [kg/s]
        rho = 1.2  # Air density [kg/m3]

        tau = cb.residence_time_mdot(V, mdot, rho)

        # τ = 0.1*1.2/0.5 = 0.24 s = 240 ms
        assert abs(tau - 0.24) < 1e-10

        # Typical combustor residence time: 10-500 ms
        assert 0.01 < tau < 1.0

    def test_residence_time_units_seconds(self):
        """Verify residence time returns seconds."""
        V = 5.0  # Volume [m3]
        Q = 2.0  # Flow rate [m3/s]

        tau = cb.residence_time(V, Q)

        # τ = 5/2 = 2.5 s
        assert abs(tau - 2.5) < 1e-10
        assert tau > 0


class TestSpaceVelocity:
    """Test space velocity calculations."""

    def test_space_velocity_basic(self):
        """Test basic space velocity calculation."""
        Q = 0.1  # Volumetric flow rate [m3/s]
        V = 1.0  # Volume [m3]

        SV = cb.space_velocity(Q, V)

        # SV = Q/V = 0.1/1.0 = 0.1 1/s
        assert abs(SV - 0.1) < 1e-10

        # Units check: should be in 1/s
        assert SV > 0
        assert isinstance(SV, float)

    def test_space_velocity_inverse_residence_time(self):
        """Test space velocity is inverse of residence time."""
        Q = 0.5  # Flow rate [m3/s]
        V = 2.0  # Volume [m3]

        SV = cb.space_velocity(Q, V)
        tau = cb.residence_time(V, Q)

        # SV = 1/τ
        assert abs(SV * tau - 1.0) < 1e-10
        assert abs(SV - 1.0 / tau) < 1e-10

    def test_space_velocity_high_flow(self):
        """Test space velocity with high flow rate."""
        Q = 10.0  # High flow rate [m3/s]
        V = 0.5  # Small volume [m3]

        SV = cb.space_velocity(Q, V)

        # SV = 10/0.5 = 20 1/s
        assert abs(SV - 20.0) < 1e-10
        assert SV > 0

    def test_space_velocity_units_inverse_seconds(self):
        """Verify space velocity returns 1/s."""
        Q = 2.0  # Flow rate [m3/s]
        V = 5.0  # Volume [m3]

        SV = cb.space_velocity(Q, V)

        # SV = 2/5 = 0.4 1/s
        assert abs(SV - 0.4) < 1e-10

        # Typical range for reactors: 0.1-100 1/s
        assert 0.01 < SV < 1000.0


class TestIntegration:
    """Integration tests combining multiple geometry functions."""

    def test_pipe_flow_residence_time(self):
        """Test complete pipe flow residence time calculation."""
        # Circular pipe
        D = 0.05  # Diameter [m]
        L = 2.0  # Length [m]
        v = 5.0  # Velocity [m/s]

        # Calculate volume
        A = 3.14159 * (D / 2) ** 2  # Area [m2]
        V = A * L  # Volume [m3]

        # Calculate flow rate
        Q = A * v  # Volumetric flow rate [m3/s]

        # Calculate residence time
        tau = cb.residence_time(V, Q)

        # τ = L/v (for constant area pipe)
        expected_tau = L / v
        assert abs(tau - expected_tau) < 0.001

    def test_reactor_design_parameters(self):
        """Test reactor design with geometry utilities."""
        # Rectangular reactor
        a = 0.2  # Width [m]
        b = 0.3  # Height [m]
        L = 1.0  # Length [m]
        Q = 0.01  # Flow rate [m3/s]

        # Calculate hydraulic diameter
        Dh = cb.hydraulic_diameter_rect(a, b)

        # Calculate volume and residence time
        V = a * b * L
        tau = cb.residence_time(V, Q)
        SV = cb.space_velocity(Q, V)

        # Verify relationships
        assert abs(SV * tau - 1.0) < 1e-10
        assert Dh > 0
        assert tau > 0
        assert SV > 0

    def test_annular_combustor(self):
        """Test annular combustor geometry."""
        # Annular combustor
        D_outer = 0.5  # Outer diameter [m]
        D_inner = 0.4  # Inner diameter [m]
        L = 0.3  # Length [m]
        mdot = 2.0  # Mass flow [kg/s]
        rho = 1.2  # Density [kg/m3]

        # Calculate hydraulic diameter
        Dh = cb.hydraulic_diameter_annulus(D_outer, D_inner)

        # Calculate volume (annular area * length)
        A_annular = 3.14159 * ((D_outer / 2) ** 2 - (D_inner / 2) ** 2)
        V = A_annular * L

        # Calculate residence time
        tau = cb.residence_time_mdot(V, mdot, rho)

        # All should be positive and reasonable
        assert Dh > 0
        assert V > 0
        assert tau > 0
        assert 0.01 < tau < 1.0  # Typical combustor range


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
