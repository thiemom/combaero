"""Tests for acoustic properties utilities.

Verifies the AcousticProperties dataclass and acoustic utility functions.
"""

import numpy as np
import pytest

import combaero as cb


class TestAcousticProperties:
    """Test AcousticProperties struct/class behavior."""

    def test_acoustic_properties_returns_struct(self):
        """Test that acoustic_properties returns AcousticProperties object."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)
        assert isinstance(props, cb.AcousticProperties)

    def test_acoustic_properties_has_all_attributes(self):
        """Test that AcousticProperties has all expected attributes."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)

        assert hasattr(props, "wavelength")
        assert hasattr(props, "frequency")
        assert hasattr(props, "impedance")
        assert hasattr(props, "particle_velocity")
        assert hasattr(props, "spl")

    def test_acoustic_properties_readonly(self):
        """Test that AcousticProperties attributes are read-only."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)

        with pytest.raises(AttributeError):
            props.wavelength = 999.0

    def test_acoustic_properties_repr(self):
        """Test that AcousticProperties has nice repr."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)
        repr_str = repr(props)

        assert "AcousticProperties" in repr_str
        assert "f=" in repr_str
        assert "Hz" in repr_str


class TestAcousticPropertiesCalculations:
    """Test acoustic_properties calculations."""

    def test_wavelength_calculation(self):
        """Test that wavelength = c / f."""
        f = 1000  # Hz
        c = 340  # m/s
        props = cb.acoustic_properties(f=f, rho=1.2, c=c)

        wavelength_expected = c / f
        assert abs(props.wavelength - wavelength_expected) < 1e-15

    def test_impedance_calculation(self):
        """Test that impedance = rho * c."""
        rho = 1.2  # kg/m^3
        c = 340  # m/s
        props = cb.acoustic_properties(f=1000, rho=rho, c=c)

        impedance_expected = rho * c
        assert abs(props.impedance - impedance_expected) < 1e-15

    def test_particle_velocity_calculation(self):
        """Test that particle_velocity = p_rms / (rho * c)."""
        rho = 1.2  # kg/m^3
        c = 340  # m/s
        p_rms = 1.0  # Pa
        props = cb.acoustic_properties(f=1000, rho=rho, c=c, p_rms=p_rms)

        particle_velocity_expected = p_rms / (rho * c)
        assert abs(props.particle_velocity - particle_velocity_expected) < 1e-15

    def test_spl_calculation(self):
        """Test that SPL = 20 * log10(p_rms / p_ref)."""
        p_rms = 1.0  # Pa
        p_ref = 20e-6  # Pa (standard reference)
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=p_rms, p_ref=p_ref)

        spl_expected = 20.0 * np.log10(p_rms / p_ref)
        assert abs(props.spl - spl_expected) < 1e-12

    def test_frequency_preserved(self):
        """Test that frequency is preserved in output."""
        f = 1234.5  # Hz
        props = cb.acoustic_properties(f=f, rho=1.2, c=340)

        assert props.frequency == f


class TestAcousticPropertiesDefaultValues:
    """Test default parameter values."""

    def test_default_p_rms_gives_0dB(self):
        """Test that default p_rms = 20e-6 Pa gives 0 dB."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)
        assert abs(props.spl - 0.0) < 1e-12

    def test_custom_p_ref(self):
        """Test custom reference pressure."""
        p_rms = 1.0  # Pa
        p_ref = 1.0  # Pa (custom reference)
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=p_rms, p_ref=p_ref)

        # SPL should be 0 dB when p_rms = p_ref
        assert abs(props.spl - 0.0) < 1e-12


class TestUtilityFunctions:
    """Test utility functions."""

    def test_wavelength(self):
        """Test wavelength calculation."""
        f = 1000  # Hz
        c = 340  # m/s
        wavelength = cb.wavelength(f, c)

        wavelength_expected = c / f
        assert abs(wavelength - wavelength_expected) < 1e-15

    def test_frequency_from_wavelength(self):
        """Test frequency from wavelength."""
        wavelength = 0.34  # m
        c = 340  # m/s
        f = cb.frequency_from_wavelength(wavelength, c)

        f_expected = c / wavelength
        assert abs(f - f_expected) < 1e-15

    def test_acoustic_impedance(self):
        """Test acoustic impedance calculation."""
        rho = 1.2  # kg/m^3
        c = 340  # m/s
        Z = cb.acoustic_impedance(rho, c)

        Z_expected = rho * c
        assert abs(Z - Z_expected) < 1e-15

    def test_sound_pressure_level(self):
        """Test SPL calculation."""
        p_rms = 1.0  # Pa
        p_ref = 20e-6  # Pa
        spl = cb.sound_pressure_level(p_rms, p_ref)

        spl_expected = 20.0 * np.log10(p_rms / p_ref)
        assert abs(spl - spl_expected) < 1e-12

    def test_sound_pressure_level_default_ref(self):
        """Test SPL with default reference pressure."""
        p_rms = 1.0  # Pa
        spl = cb.sound_pressure_level(p_rms)

        spl_expected = 20.0 * np.log10(p_rms / 20e-6)
        assert abs(spl - spl_expected) < 1e-12

    def test_particle_velocity(self):
        """Test particle velocity calculation."""
        p = 1.0  # Pa
        rho = 1.2  # kg/m^3
        c = 340  # m/s
        u = cb.particle_velocity(p, rho, c)

        u_expected = p / (rho * c)
        assert abs(u - u_expected) < 1e-15

    def test_round_trip_wavelength_frequency(self):
        """Test round-trip: f → λ → f."""
        f_original = 1000  # Hz
        c = 340  # m/s

        wavelength = cb.wavelength(f_original, c)
        f_recovered = cb.frequency_from_wavelength(wavelength, c)

        assert abs(f_recovered - f_original) < 1e-12


class TestTypicalConditions:
    """Test with typical acoustic conditions."""

    def test_air_at_standard_conditions(self):
        """Test air at 20degC, 1 atm."""
        f = 1000  # Hz
        rho = 1.2  # kg/m^3
        c = 343  # m/s at 20degC
        p_rms = 0.02  # Pa (quiet room)

        props = cb.acoustic_properties(f=f, rho=rho, c=c, p_rms=p_rms)

        # Verify all properties match individual calls
        assert props.wavelength == pytest.approx(cb.wavelength(f, c), rel=1e-15)
        assert props.impedance == pytest.approx(cb.acoustic_impedance(rho, c), rel=1e-15)
        assert props.spl == pytest.approx(cb.sound_pressure_level(p_rms), rel=1e-15)

    def test_combustor_conditions(self):
        """Test hot gas in combustor."""
        f = 500  # Hz (typical combustor instability)
        T = 1500  # K
        P = 15e5  # Pa (15 bar)

        # Approximate properties for hot air
        rho = P / (287 * T)  # ~3.5 kg/m^3
        c = np.sqrt(1.3 * 287 * T)  # ~775 m/s

        props = cb.acoustic_properties(f=f, rho=rho, c=c, p_rms=1000)

        # Verify all properties match individual calls
        assert props.wavelength == pytest.approx(cb.wavelength(f, c), rel=1e-15)
        assert props.impedance == pytest.approx(cb.acoustic_impedance(rho, c), rel=1e-15)

    def test_water_acoustics(self):
        """Test underwater acoustics."""
        f = 1000  # Hz
        rho = 1000  # kg/m^3
        c = 1500  # m/s (speed of sound in water)

        props = cb.acoustic_properties(f=f, rho=rho, c=c, p_rms=1.0)

        # Verify all properties match individual calls
        assert props.wavelength == pytest.approx(cb.wavelength(f, c), rel=1e-15)
        assert props.impedance == pytest.approx(cb.acoustic_impedance(rho, c), rel=1e-15)
        assert abs(props.impedance - 1.5e6) < 1000


class TestFrequencyRange:
    """Test across frequency range."""

    def test_low_frequency(self):
        """Test low frequency (20 Hz)."""
        props = cb.acoustic_properties(f=20, rho=1.2, c=340)
        assert props.wavelength == pytest.approx(17.0, rel=1e-10)

    def test_mid_frequency(self):
        """Test mid frequency (1 kHz)."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340)
        assert props.wavelength == pytest.approx(0.34, rel=1e-10)

    def test_high_frequency(self):
        """Test high frequency (20 kHz)."""
        props = cb.acoustic_properties(f=20000, rho=1.2, c=340)
        assert props.wavelength == pytest.approx(0.017, rel=1e-10)


class TestSPLLevels:
    """Test SPL at various levels."""

    def test_threshold_of_hearing(self):
        """Test threshold of hearing (0 dB)."""
        p_rms = 20e-6  # Pa
        spl = cb.sound_pressure_level(p_rms)
        assert abs(spl - 0.0) < 1e-12

    def test_quiet_room(self):
        """Test quiet room (~40 dB)."""
        p_rms = 20e-6 * 10 ** (40 / 20)  # Pa
        spl = cb.sound_pressure_level(p_rms)
        assert abs(spl - 40.0) < 1e-10

    def test_normal_conversation(self):
        """Test normal conversation (~60 dB)."""
        p_rms = 20e-6 * 10 ** (60 / 20)  # Pa
        spl = cb.sound_pressure_level(p_rms)
        assert abs(spl - 60.0) < 1e-10

    def test_loud_music(self):
        """Test loud music (~100 dB)."""
        p_rms = 20e-6 * 10 ** (100 / 20)  # Pa
        spl = cb.sound_pressure_level(p_rms)
        assert abs(spl - 100.0) < 1e-10

    def test_threshold_of_pain(self):
        """Test threshold of pain (~130 dB)."""
        p_rms = 20e-6 * 10 ** (130 / 20)  # Pa
        spl = cb.sound_pressure_level(p_rms)
        assert abs(spl - 130.0) < 1e-10


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_negative_frequency_raises_error(self):
        """Test that negative frequency raises error."""
        with pytest.raises(ValueError, match="frequency must be positive"):
            cb.acoustic_properties(f=-100, rho=1.2, c=340)

    def test_zero_frequency_raises_error(self):
        """Test that zero frequency raises error."""
        with pytest.raises(ValueError, match="frequency must be positive"):
            cb.acoustic_properties(f=0, rho=1.2, c=340)

    def test_negative_density_raises_error(self):
        """Test that negative density raises error."""
        with pytest.raises(ValueError, match="density must be positive"):
            cb.acoustic_properties(f=1000, rho=-1.2, c=340)

    def test_negative_speed_of_sound_raises_error(self):
        """Test that negative speed of sound raises error."""
        with pytest.raises(ValueError, match="speed of sound must be positive"):
            cb.acoustic_properties(f=1000, rho=1.2, c=-340)

    def test_negative_p_rms_raises_error(self):
        """Test that negative p_rms raises error."""
        with pytest.raises(ValueError, match="p_rms must be non-negative"):
            cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=-1.0)

    def test_zero_p_rms_gives_minus_infinity_spl(self):
        """Test that zero p_rms gives -∞ dB."""
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=0.0)
        assert props.spl == -np.inf

    def test_negative_p_ref_raises_error(self):
        """Test that negative p_ref raises error."""
        with pytest.raises(ValueError, match="p_ref must be positive"):
            cb.acoustic_properties(f=1000, rho=1.2, c=340, p_ref=-20e-6)


class TestConsistency:
    """Test consistency between composite and utility functions."""

    def test_wavelength_matches_utility(self):
        """Test that props.wavelength matches wavelength()."""
        f = 1000
        c = 340
        props = cb.acoustic_properties(f=f, rho=1.2, c=c)
        wavelength_util = cb.wavelength(f, c)

        assert props.wavelength == wavelength_util

    def test_impedance_matches_utility(self):
        """Test that props.impedance matches acoustic_impedance()."""
        rho = 1.2
        c = 340
        props = cb.acoustic_properties(f=1000, rho=rho, c=c)
        impedance_util = cb.acoustic_impedance(rho, c)

        assert props.impedance == impedance_util

    def test_particle_velocity_matches_utility(self):
        """Test that props.particle_velocity matches particle_velocity()."""
        rho = 1.2
        c = 340
        p_rms = 1.0
        props = cb.acoustic_properties(f=1000, rho=rho, c=c, p_rms=p_rms)
        u_util = cb.particle_velocity(p_rms, rho, c)

        assert props.particle_velocity == u_util

    def test_spl_matches_utility(self):
        """Test that props.spl matches sound_pressure_level()."""
        p_rms = 1.0
        p_ref = 20e-6
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=p_rms, p_ref=p_ref)
        spl_util = cb.sound_pressure_level(p_rms, p_ref)

        assert props.spl == spl_util


class TestPhysicalRelationships:
    """Test physical relationships between properties."""

    def test_wavelength_frequency_product(self):
        """Test that wavelength * frequency = c."""
        f = 1000
        c = 340
        props = cb.acoustic_properties(f=f, rho=1.2, c=c)

        assert abs(props.wavelength * props.frequency - c) < 1e-12

    def test_impedance_particle_velocity_product(self):
        """Test that impedance * particle_velocity = p_rms."""
        p_rms = 1.0
        props = cb.acoustic_properties(f=1000, rho=1.2, c=340, p_rms=p_rms)

        assert abs(props.impedance * props.particle_velocity - p_rms) < 1e-12


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
