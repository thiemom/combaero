"""Tests for can-annular and annular duct acoustic solvers."""

import pytest
import combaero as cb


class TestCanAnnularGeometry:
    """Test CanAnnularGeometry structure."""

    def test_geometry_creation(self):
        """Test creating can-annular geometry."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        assert geom.n_cans == 24
        assert geom.length_can == 0.5
        assert geom.area_can == 0.01
        assert geom.radius_plenum == 0.5
        assert geom.area_plenum == 0.02


class TestBlochMode:
    """Test BlochMode structure."""

    def test_bloch_mode_creation(self):
        """Test creating Bloch mode."""
        mode = cb.BlochMode()
        mode.m_azimuthal = 12
        mode.frequency = 250.0
        mode.n_cans = 24

        assert mode.m_azimuthal == 12
        assert mode.frequency == 250.0
        assert mode.n_cans == 24

    def test_symmetry_type(self):
        """Test symmetry type classification."""
        # Push-push (m=0)
        mode = cb.BlochMode()
        mode.m_azimuthal = 0
        mode.frequency = 250.0
        mode.n_cans = 24
        assert "Push-Push" in mode.symmetry_type()

        # Push-pull (m=N/2)
        mode.m_azimuthal = 12
        mode.n_cans = 24
        assert "Push-Pull" in mode.symmetry_type()

        # Spinning (other m)
        mode.m_azimuthal = 6
        mode.n_cans = 24
        assert "Spinning" in mode.symmetry_type()


class TestSingleCanValidation:
    """Test single can eigenmode solver against analytical solutions."""

    def test_single_can_half_wave_modes(self):
        """Test single can finds half-wave resonances (closed-closed)."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 1
        geom.length_can = 0.5  # m
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        c = 500.0  # m/s
        modes = cb.can_annular_eigenmodes(geom, c, c, 1.0, 1.0, f_max=1500)

        # Should find modes at n*c/(2L) = n*500/(2*0.5) = n*500 Hz
        # Expected: 500, 1000, 1500 Hz (within search range)
        assert len(modes) >= 2

        # Check first mode (fundamental)
        f_expected_1 = c / (2.0 * geom.length_can)  # 500 Hz
        assert abs(modes[0].frequency - f_expected_1) < 10  # Within 10 Hz

        # Check second mode
        if len(modes) >= 2:
            f_expected_2 = 2 * c / (2.0 * geom.length_can)  # 1000 Hz
            assert abs(modes[1].frequency - f_expected_2) < 10

    def test_single_can_all_modes_push_push(self):
        """Test that single can modes are all push-push (m=0)."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 1
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 500, 1.0, 1.0, f_max=1000)

        # All modes should be m=0 for single can
        for mode in modes:
            assert mode.m_azimuthal == 0
            assert "Push-Push" in mode.symmetry_type()

    def test_single_can_frequency_scaling(self):
        """Test that mode frequencies scale with speed of sound."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 1
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        # Test with c=500 m/s
        modes_500 = cb.can_annular_eigenmodes(geom, 500, 500, 1.0, 1.0, f_max=1000)

        # Test with c=600 m/s (20% higher)
        modes_600 = cb.can_annular_eigenmodes(geom, 600, 600, 1.0, 1.0, f_max=1200)

        # Frequencies should scale proportionally
        if len(modes_500) > 0 and len(modes_600) > 0:
            ratio = modes_600[0].frequency / modes_500[0].frequency
            assert abs(ratio - 1.2) < 0.05  # Within 5% of expected 20% increase


class TestMultiCanModeSplitting:
    """Test multi-can system shows mode splitting."""

    def test_24_can_finds_multiple_azimuthal_modes(self):
        """Test that 24-can system finds modes with different m."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        # Should find modes with different azimuthal numbers
        azimuthal_modes = set(mode.m_azimuthal for mode in modes)
        assert len(azimuthal_modes) > 1  # Multiple azimuthal modes

        # Should have modes between 0 and N/2 = 12
        for m in azimuthal_modes:
            assert 0 <= m <= 12

    def test_mode_splitting_exists(self):
        """Test that coupling causes mode splitting."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        # Should find multiple modes in similar frequency range (mode splitting)
        if len(modes) >= 2:
            # Check that modes are clustered (not all at exact same frequency)
            freq_range = max(m.frequency for m in modes) - min(m.frequency for m in modes)
            assert freq_range > 10  # At least 10 Hz spread due to splitting

    def test_push_pull_mode_exists(self):
        """Test that push-pull mode (m=N/2) is found."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        # Should find push-pull mode (m=12 for N=24)
        push_pull_modes = [m for m in modes if m.m_azimuthal == 12]
        assert len(push_pull_modes) > 0

        # Verify symmetry type
        assert "Push-Pull" in push_pull_modes[0].symmetry_type()


class TestParameterValidation:
    """Test parameter validation for can-annular solver."""

    def test_invalid_n_cans(self):
        """Test that invalid number of cans raises error."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 0  # Invalid
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        with pytest.raises(RuntimeError, match="cans must be positive"):
            cb.can_annular_eigenmodes(geom, 500, 500, 1.0, 1.0, f_max=1000)

    def test_invalid_geometry(self):
        """Test that invalid geometry raises error."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = -0.5  # Invalid
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        with pytest.raises(RuntimeError, match="dimensions must be positive"):
            cb.can_annular_eigenmodes(geom, 500, 500, 1.0, 1.0, f_max=1000)

    def test_invalid_sound_speed(self):
        """Test that invalid sound speed raises error."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        with pytest.raises(RuntimeError, match="Speed of sound"):
            cb.can_annular_eigenmodes(geom, -500, 500, 1.0, 1.0, f_max=1000)

    def test_invalid_density(self):
        """Test that invalid density raises error."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        with pytest.raises(RuntimeError, match="Density"):
            cb.can_annular_eigenmodes(geom, 500, 500, -1.0, 1.0, f_max=1000)


class TestPhysicalConsistency:
    """Test physical consistency of results."""

    def test_modes_sorted_by_frequency(self):
        """Test that modes are returned sorted by frequency."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        # Check sorting
        for i in range(len(modes) - 1):
            assert modes[i].frequency <= modes[i + 1].frequency

    def test_frequencies_positive(self):
        """Test that all frequencies are positive."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        for mode in modes:
            assert mode.frequency > 0

    def test_frequencies_within_range(self):
        """Test that all frequencies are within specified range."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        f_max = 600
        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=f_max)

        for mode in modes:
            assert mode.frequency <= f_max

    def test_no_duplicate_modes(self):
        """Test that no duplicate modes are returned."""
        geom = cb.CanAnnularGeometry()
        geom.n_cans = 24
        geom.length_can = 0.5
        geom.area_can = 0.01
        geom.radius_plenum = 0.5
        geom.area_plenum = 0.02

        modes = cb.can_annular_eigenmodes(geom, 500, 600, 1.0, 1.2, f_max=600)

        # Check for duplicates (same m and frequency within 1 Hz)
        for i in range(len(modes)):
            for j in range(i + 1, len(modes)):
                if modes[i].m_azimuthal == modes[j].m_azimuthal:
                    assert abs(modes[i].frequency - modes[j].frequency) >= 1.0
