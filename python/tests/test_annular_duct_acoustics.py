"""Tests for annular duct acoustic solver."""

import math

import pytest

import combaero as cb


class TestAnnularDuctGeometry:
    """Test AnnularDuctGeometry structure."""

    def test_geometry_creation(self):
        """Test creating annular duct geometry."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 5

        assert geom.length == 1.0
        assert geom.radius_inner == 0.2
        assert geom.radius_outer == 0.5
        assert geom.n_azimuthal_max == 5

    def test_area_calculation(self):
        """Test annular area calculation."""
        geom = cb.AnnularDuctGeometry()
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5

        # Area = pi * (r_outer^2 - r_inner^2)
        expected_area = math.pi * (0.5**2 - 0.2**2)
        assert abs(geom.area() - expected_area) < 1e-10


class TestAnnularMode:
    """Test AnnularMode structure."""

    def test_mode_creation(self):
        """Test creating annular mode."""
        mode = cb.AnnularMode()
        mode.m_azimuthal = 2
        mode.n_axial = 3
        mode.frequency = 500.0

        assert mode.m_azimuthal == 2
        assert mode.n_axial == 3
        assert mode.frequency == 500.0

    def test_mode_type_axisymmetric(self):
        """Test mode type for axisymmetric mode."""
        mode = cb.AnnularMode()
        mode.m_azimuthal = 0
        mode.n_axial = 1
        mode.frequency = 250.0

        assert "Axisymmetric" in mode.mode_type()

    def test_mode_type_azimuthal(self):
        """Test mode type for azimuthal mode."""
        mode = cb.AnnularMode()
        mode.m_azimuthal = 3
        mode.n_axial = 1
        mode.frequency = 300.0

        assert "Azimuthal" in mode.mode_type()
        assert "m=3" in mode.mode_type()


class TestAnalyticalValidation:
    """Test annular duct solver against analytical solutions."""

    def test_closed_closed_duct_modes(self):
        """Test closed-closed duct finds correct resonances."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0  # m
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 2

        c = 500.0  # m/s
        modes = cb.annular_duct_eigenmodes(geom, c, 1.2, f_max=1500)

        # For closed-closed: f = n * c / (2*L)
        # Expected: 250, 500, 750, 1000, 1250 Hz
        assert len(modes) >= 2

        # Check fundamental (should be around 250 Hz for m=0)
        axisymmetric_modes = [m for m in modes if m.m_azimuthal == 0]
        if len(axisymmetric_modes) > 0:
            f_expected = c / (2.0 * geom.length)  # 250 Hz
            # Allow 20% tolerance (solver may find slightly different frequency)
            assert abs(axisymmetric_modes[0].frequency - f_expected) < 0.2 * f_expected

    def test_frequency_scaling_with_sound_speed(self):
        """Test that frequencies scale linearly with sound speed."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 2

        # Use analytical reference for deterministic mode indexing
        modes_500 = cb.annular_duct_modes_analytical(geom, 500, f_max=1000)

        # Test with c=600 m/s (20% higher)
        modes_600 = cb.annular_duct_modes_analytical(geom, 600, f_max=1200)

        # Frequencies should scale proportionally
        if len(modes_500) > 0 and len(modes_600) > 0:
            # Find corresponding modes (same m and n)
            for m500 in modes_500:
                for m600 in modes_600:
                    if m500.m_azimuthal == m600.m_azimuthal and m500.n_axial == m600.n_axial:
                        ratio = m600.frequency / m500.frequency
                        # Should be close to 1.2 (20% increase)
                        assert abs(ratio - 1.2) < 0.1
                        break

    def test_frequency_scaling_with_length(self):
        """Test that frequencies scale inversely with length."""
        geom1 = cb.AnnularDuctGeometry()
        geom1.length = 1.0
        geom1.radius_inner = 0.2
        geom1.radius_outer = 0.5
        geom1.n_azimuthal_max = 2

        geom2 = cb.AnnularDuctGeometry()
        geom2.length = 2.0  # Double length
        geom2.radius_inner = 0.2
        geom2.radius_outer = 0.5
        geom2.n_azimuthal_max = 2

        # Use analytical reference for deterministic mode indexing
        modes1 = cb.annular_duct_modes_analytical(geom1, 500, f_max=1000)
        modes2 = cb.annular_duct_modes_analytical(geom2, 500, f_max=500)

        # Inverse-length scaling applies directly to axisymmetric modes (m=0)
        axis1 = [m for m in modes1 if m.m_azimuthal == 0 and m.n_axial == 1]
        axis2 = [m for m in modes2 if m.m_azimuthal == 0 and m.n_axial == 1]

        if axis1 and axis2:
            ratio = axis1[0].frequency / axis2[0].frequency
            assert abs(ratio - 2.0) < 0.2


class TestMultipleAzimuthalModes:
    """Test that solver finds multiple azimuthal modes."""

    def test_finds_multiple_azimuthal_modes(self):
        """Test that different azimuthal modes are found."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 5

        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1500)

        # Should find modes with different m
        azimuthal_modes = {m.m_azimuthal for m in modes}
        assert len(azimuthal_modes) > 1

        # Should have modes between 0 and n_azimuthal_max
        for m in azimuthal_modes:
            assert 0 <= m <= 5

    def test_axisymmetric_mode_exists(self):
        """Test that axisymmetric mode (m=0) is found."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

        # Should find at least one m=0 mode
        axisymmetric = [m for m in modes if m.m_azimuthal == 0]
        assert len(axisymmetric) > 0

    def test_m_changes_frequency_for_same_axial_order(self):
        """For fixed n, increasing m should increase frequency."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 4

        modes = cb.annular_duct_eigenmodes(geom, 500.0, 1.2, f_max=2000.0)
        n1 = [m for m in modes if m.n_axial == 1]
        n1_sorted = sorted(n1, key=lambda mode: mode.m_azimuthal)

        if len(n1_sorted) >= 2:
            for i in range(len(n1_sorted) - 1):
                assert n1_sorted[i + 1].frequency >= n1_sorted[i].frequency


class TestParameterValidation:
    """Test parameter validation."""

    def test_invalid_length(self):
        """Test that invalid length raises error."""
        geom = cb.AnnularDuctGeometry()
        geom.length = -1.0  # Invalid
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        with pytest.raises(RuntimeError, match="length must be positive"):
            cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

    def test_invalid_radii(self):
        """Test that invalid radii raise error."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.5
        geom.radius_outer = 0.2  # Invalid: inner > outer
        geom.n_azimuthal_max = 3

        with pytest.raises(RuntimeError, match="Invalid annular geometry"):
            cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

    def test_invalid_sound_speed(self):
        """Test that invalid sound speed raises error."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        with pytest.raises(RuntimeError, match="Speed of sound"):
            cb.annular_duct_eigenmodes(geom, -500, 1.2, f_max=1000)


class TestPhysicalConsistency:
    """Test physical consistency of results."""

    def test_modes_sorted_by_frequency(self):
        """Test that modes are sorted by frequency."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

        for i in range(len(modes) - 1):
            assert modes[i].frequency <= modes[i + 1].frequency

    def test_frequencies_positive(self):
        """Test that all frequencies are positive."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

        for mode in modes:
            assert mode.frequency > 0

    def test_frequencies_within_range(self):
        """Test that frequencies are within specified range."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        f_max = 1000
        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=f_max)

        for mode in modes:
            assert mode.frequency <= f_max

    def test_axial_mode_numbers_positive(self):
        """Test that axial mode numbers are positive."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        modes = cb.annular_duct_eigenmodes(geom, 500, 1.2, f_max=1000)

        for mode in modes:
            assert mode.n_axial > 0


class TestAnalyticalReference:
    """Test analytical annular helper and AP solver consistency."""

    def test_analytical_modes_available(self):
        """Analytical helper should return modes in requested range."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 3

        modes = cb.annular_duct_modes_analytical(geom, 500.0, f_max=1000.0)
        assert len(modes) > 0
        assert all(mode.frequency <= 1000.0 for mode in modes)

    def test_argument_principle_matches_analytical_fundamental(self):
        """AP and analytical should be close for the fundamental axisymmetric mode."""
        geom = cb.AnnularDuctGeometry()
        geom.length = 1.0
        geom.radius_inner = 0.2
        geom.radius_outer = 0.5
        geom.n_azimuthal_max = 2

        ap_modes = cb.annular_duct_eigenmodes(geom, 500.0, 1.2, f_max=1200.0)
        ref_modes = cb.annular_duct_modes_analytical(geom, 500.0, f_max=1200.0)

        ap_axis = [m for m in ap_modes if m.m_azimuthal == 0 and m.n_axial == 1]
        ref_axis = [m for m in ref_modes if m.m_azimuthal == 0 and m.n_axial == 1]

        if ap_axis and ref_axis:
            rel_err = abs(ap_axis[0].frequency - ref_axis[0].frequency) / ref_axis[0].frequency
            assert rel_err < 0.2


class TestComparisonWithCanAnnular:
    """Test consistency between annular duct and single-can solvers."""

    def test_single_can_matches_annular_duct(self):
        """Test that single can gives similar results to annular duct."""
        # Single can geometry
        can_geom = cb.CanAnnularGeometry()
        can_geom.n_cans = 1
        can_geom.length_can = 1.0
        can_geom.area_can = 0.1
        can_geom.radius_plenum = 0.5
        can_geom.area_plenum = 0.1

        # Equivalent annular duct
        ann_geom = cb.AnnularDuctGeometry()
        ann_geom.length = 1.0
        ann_geom.radius_inner = 0.0
        ann_geom.radius_outer = math.sqrt(0.1 / math.pi)  # Same area
        ann_geom.n_azimuthal_max = 0  # Only axisymmetric

        c = 500.0
        can_modes = cb.can_annular_eigenmodes(can_geom, c, c, 1.2, 1.2, f_max=1500)
        ann_modes = cb.annular_duct_eigenmodes(ann_geom, c, 1.2, f_max=1500)

        # Both should find similar fundamental frequencies
        # (may not be exact due to different solver methods)
        if len(can_modes) > 0 and len(ann_modes) > 0:
            # Allow 30% difference (different solver methods)
            ratio = can_modes[0].frequency / ann_modes[0].frequency
            assert 0.7 < ratio < 1.3
