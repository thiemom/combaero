"""Tests for geometry utility functions.

Verifies Python bindings for hydraulic diameter and residence time calculations.
Units: Dh in [m], τ in [s], SV in [1/s].
"""

import pytest

import combaero as cb


class TestHydraulicDiameter:
    """Test hydraulic diameter calculations."""

    def test_hydraulic_diameter_circular_pipe(self):
        """Test hydraulic diameter for circular pipe."""
        # For circular pipe: Dh = D
        D = 0.1  # Diameter [m]
        A = 3.14159 * (D / 2) ** 2  # Area [m²]
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
        A = a * a  # Area [m²]
        P = 4 * a  # Perimeter [m]

        Dh = cb.hydraulic_diameter(A, P)

        # Dh = 4*A/P = 4*a²/(4*a) = a
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
        A = 0.01  # Area [m²]
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
        V = 1.0  # Volume [m³]
        Q = 0.1  # Volumetric flow rate [m³/s]

        tau = cb.residence_time(V, Q)

        # τ = V/Q = 1.0/0.1 = 10 s
        assert abs(tau - 10.0) < 1e-10

        # Units check: should be in seconds
        assert tau > 0
        assert isinstance(tau, float)

    def test_residence_time_fast_flow(self):
        """Test residence time with fast flow."""
        V = 0.5  # Volume [m³]
        Q = 10.0  # High flow rate [m³/s]

        tau = cb.residence_time(V, Q)

        # τ = 0.5/10 = 0.05 s
        assert abs(tau - 0.05) < 1e-10
        assert tau > 0

    def test_residence_time_mdot_basic(self):
        """Test residence time from mass flow rate."""
        V = 1.0  # Volume [m³]
        mdot = 1.2  # Mass flow rate [kg/s]
        rho = 1.2  # Density (air) [kg/m³]

        tau = cb.residence_time_mdot(V, mdot, rho)

        # τ = V*ρ/ṁ = 1.0*1.2/1.2 = 1.0 s
        assert abs(tau - 1.0) < 1e-10
        assert tau > 0

    def test_residence_time_mdot_vs_volumetric(self):
        """Test consistency between mass and volumetric methods."""
        V = 2.0  # Volume [m³]
        rho = 1.2  # Density [kg/m³]
        Q = 0.5  # Volumetric flow rate [m³/s]
        mdot = rho * Q  # Mass flow rate [kg/s]

        tau_vol = cb.residence_time(V, Q)
        tau_mass = cb.residence_time_mdot(V, mdot, rho)

        # Should give same result
        assert abs(tau_vol - tau_mass) < 1e-10

    def test_residence_time_combustor(self):
        """Test realistic combustor residence time."""
        V = 0.1  # Combustor volume [m³]
        mdot = 0.5  # Air mass flow [kg/s]
        rho = 1.2  # Air density [kg/m³]

        tau = cb.residence_time_mdot(V, mdot, rho)

        # τ = 0.1*1.2/0.5 = 0.24 s = 240 ms
        assert abs(tau - 0.24) < 1e-10

        # Typical combustor residence time: 10-500 ms
        assert 0.01 < tau < 1.0

    def test_residence_time_units_seconds(self):
        """Verify residence time returns seconds."""
        V = 5.0  # Volume [m³]
        Q = 2.0  # Flow rate [m³/s]

        tau = cb.residence_time(V, Q)

        # τ = 5/2 = 2.5 s
        assert abs(tau - 2.5) < 1e-10
        assert tau > 0


class TestSpaceVelocity:
    """Test space velocity calculations."""

    def test_space_velocity_basic(self):
        """Test basic space velocity calculation."""
        Q = 0.1  # Volumetric flow rate [m³/s]
        V = 1.0  # Volume [m³]

        SV = cb.space_velocity(Q, V)

        # SV = Q/V = 0.1/1.0 = 0.1 1/s
        assert abs(SV - 0.1) < 1e-10

        # Units check: should be in 1/s
        assert SV > 0
        assert isinstance(SV, float)

    def test_space_velocity_inverse_residence_time(self):
        """Test space velocity is inverse of residence time."""
        Q = 0.5  # Flow rate [m³/s]
        V = 2.0  # Volume [m³]

        SV = cb.space_velocity(Q, V)
        tau = cb.residence_time(V, Q)

        # SV = 1/τ
        assert abs(SV * tau - 1.0) < 1e-10
        assert abs(SV - 1.0 / tau) < 1e-10

    def test_space_velocity_high_flow(self):
        """Test space velocity with high flow rate."""
        Q = 10.0  # High flow rate [m³/s]
        V = 0.5  # Small volume [m³]

        SV = cb.space_velocity(Q, V)

        # SV = 10/0.5 = 20 1/s
        assert abs(SV - 20.0) < 1e-10
        assert SV > 0

    def test_space_velocity_units_inverse_seconds(self):
        """Verify space velocity returns 1/s."""
        Q = 2.0  # Flow rate [m³/s]
        V = 5.0  # Volume [m³]

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
        A = 3.14159 * (D / 2) ** 2  # Area [m²]
        V = A * L  # Volume [m³]

        # Calculate flow rate
        Q = A * v  # Volumetric flow rate [m³/s]

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
        Q = 0.01  # Flow rate [m³/s]

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
        rho = 1.2  # Density [kg/m³]

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
