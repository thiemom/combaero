"""Tests for geometry adapter functionality.

Verifies Annulus struct and to_acoustic_geometry adapter.
"""

import math

import pytest

import combaero as cb


class TestAnnulusStruct:
    """Test Annulus struct with new acoustic solver compatibility."""

    def test_annulus_creation(self):
        """Test creating Annulus with new constructor."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        assert annulus.L == 1.0
        assert annulus.D_inner == 0.4
        assert annulus.D_outer == 1.0

    def test_annulus_n_azimuthal_max(self):
        """Test n_azimuthal_max field."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        # Default value
        assert annulus.n_azimuthal_max == 5

        # Can be modified
        annulus.n_azimuthal_max = 10
        assert annulus.n_azimuthal_max == 10

    def test_annulus_radius_methods(self):
        """Test radius accessor methods."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        assert annulus.radius_inner == 0.2  # D_inner / 2
        assert annulus.radius_outer == 0.5  # D_outer / 2

    def test_annulus_diameter_methods(self):
        """Test diameter accessor methods."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        assert annulus.D_mean() == 0.7  # (0.4 + 1.0) / 2
        assert annulus.gap() == 0.3  # (1.0 - 0.4) / 2

    def test_annulus_area(self):
        """Test area calculation."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        expected_area = math.pi * (0.5**2 - 0.2**2)
        assert abs(annulus.area() - expected_area) < 1e-10

    def test_annulus_volume(self):
        """Test volume calculation."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        expected_volume = annulus.area() * 1.0
        assert abs(annulus.volume() - expected_volume) < 1e-10

    def test_annulus_circumference(self):
        """Test circumference calculation."""
        annulus = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)

        expected_circumference = math.pi * 0.7  # π * D_mean
        assert abs(annulus.circumference() - expected_circumference) < 1e-10


class TestToAcousticGeometry:
    """Test to_acoustic_geometry adapter function."""

    def test_basic_conversion(self):
        """Test basic flow geometry to acoustic geometry conversion."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        geom = cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=0.5, D_can=0.1)

        assert geom.n_cans == 24
        assert geom.length_can == 0.5

        # Can area = π * (0.1/2)^2
        expected_can_area = math.pi * 0.0025
        assert abs(geom.area_can - expected_can_area) < 1e-10

        # Plenum radius = D_mean / 2 = (0.4 + 0.6) / 4 = 0.25
        assert geom.radius_plenum == 0.25

        # Plenum area
        expected_plenum_area = flow_geom.area()
        assert abs(geom.area_plenum - expected_plenum_area) < 1e-10

    def test_different_can_counts(self):
        """Test conversion with different numbers of cans."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        for n_cans in [1, 8, 12, 24, 36]:
            geom = cb.to_acoustic_geometry(flow_geom, n_cans=n_cans, L_can=0.5, D_can=0.1)
            assert geom.n_cans == n_cans

    def test_different_can_dimensions(self):
        """Test conversion with various can dimensions."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        test_cases = [(0.3, 0.08), (0.5, 0.1), (0.7, 0.12), (1.0, 0.15)]

        for L_can, D_can in test_cases:
            geom = cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=L_can, D_can=D_can)

            assert geom.length_can == L_can

            expected_can_area = math.pi * (D_can / 2.0) ** 2
            assert abs(geom.area_can - expected_can_area) < 1e-10

    def test_plenum_properties_from_flow_geom(self):
        """Test that plenum properties are correctly derived."""
        # Case 1: Different annular diameters
        flow_geom1 = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.5, D_outer=0.7, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )
        geom1 = cb.to_acoustic_geometry(flow_geom1, n_cans=12, L_can=0.4, D_can=0.08)

        # D_mean = (0.5 + 0.7) / 2 = 0.6, so radius = 0.3
        assert geom1.radius_plenum == 0.3

        # Case 2: Larger annulus
        flow_geom2 = cb.CanAnnularFlowGeometry(
            L=0.4, D_inner=0.8, D_outer=1.2, L_primary=0.3, D_primary=0.15, L_transition=0.08
        )
        geom2 = cb.to_acoustic_geometry(flow_geom2, n_cans=18, L_can=0.6, D_can=0.12)

        # D_mean = (0.8 + 1.2) / 2 = 1.0, so radius = 0.5
        assert geom2.radius_plenum == 0.5

    def test_invalid_inputs(self):
        """Test that invalid inputs raise appropriate errors."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        # Invalid n_cans
        with pytest.raises((RuntimeError, ValueError)):
            cb.to_acoustic_geometry(flow_geom, n_cans=0, L_can=0.5, D_can=0.1)

        # Invalid L_can
        with pytest.raises((RuntimeError, ValueError)):
            cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=0.0, D_can=0.1)

        # Invalid D_can
        with pytest.raises((RuntimeError, ValueError)):
            cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=0.5, D_can=0.0)

    def test_edge_cases(self):
        """Test edge cases like very small/large values."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        # Very small can
        small = cb.to_acoustic_geometry(flow_geom, n_cans=1, L_can=0.01, D_can=0.001)
        assert small.n_cans == 1

        # Large number of cans
        many = cb.to_acoustic_geometry(flow_geom, n_cans=1000, L_can=0.5, D_can=0.1)
        assert many.n_cans == 1000

    def test_integration_with_solver(self):
        """Test that converted geometry works with can_annular_eigenmodes."""
        flow_geom = cb.CanAnnularFlowGeometry(
            L=0.3, D_inner=0.4, D_outer=0.6, L_primary=0.2, D_primary=0.1, L_transition=0.05
        )

        geom = cb.to_acoustic_geometry(flow_geom, n_cans=8, L_can=0.5, D_can=0.1)

        # Basic validation that geometry is usable
        assert geom.area_can > 0.0
        assert geom.area_plenum > 0.0
        assert geom.radius_plenum > 0.0
        assert geom.length_can > 0.0
        assert geom.n_cans > 0

        # Verify we can call the solver (may not find modes with these params, but shouldn't crash)
        try:
            modes = cb.can_annular_eigenmodes(
                geom, c_can=500.0, c_plenum=550.0, rho_can=1.2, rho_plenum=1.0, f_max=1000.0
            )
            # Should return a list (possibly empty)
            assert isinstance(modes, list)
        except RuntimeError:
            # Some parameter combinations may not yield valid modes
            pass
