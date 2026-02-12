"""Tests for pipe roughness database functions (Phase 6)."""

import pytest
import combaero as cb


class TestPipeRoughness:
    """Test pipe roughness database functions."""

    def test_smooth_pipe(self):
        """Test smooth pipe (theoretical)."""
        eps = cb.pipe_roughness('smooth')
        assert eps == 0.0

    def test_commercial_steel(self):
        """Test commercial steel roughness."""
        eps = cb.pipe_roughness('commercial_steel')
        # White (2011): 0.045 mm = 4.5e-5 m
        assert eps == 4.5e-5
        assert abs(eps - 0.000045) < 1e-10

    def test_galvanized_iron(self):
        """Test galvanized iron roughness."""
        eps = cb.pipe_roughness('galvanized_iron')
        # Moody (1944), Crane (2009): 0.15 mm = 1.5e-4 m
        assert eps == 1.5e-4

    def test_cast_iron(self):
        """Test cast iron roughness."""
        eps = cb.pipe_roughness('cast_iron')
        # Moody (1944), White (2011): 0.26 mm = 2.6e-4 m
        assert eps == 2.6e-4

    def test_concrete(self):
        """Test concrete roughness."""
        eps = cb.pipe_roughness('concrete')
        # Moody (1944): 0.3 mm = 3.0e-4 m
        assert eps == 3.0e-4

    def test_drawn_tubing(self):
        """Test drawn tubing (very smooth)."""
        eps = cb.pipe_roughness('drawn_tubing')
        # White (2011): 0.0015 mm = 1.5e-6 m
        assert eps == 1.5e-6

    def test_pvc_plastic(self):
        """Test PVC/plastic pipes."""
        eps_pvc = cb.pipe_roughness('pvc')
        eps_plastic = cb.pipe_roughness('plastic')
        # Both should be very smooth
        assert eps_pvc == 1.5e-6
        assert eps_plastic == 1.5e-6
        assert eps_pvc == eps_plastic

    def test_riveted_steel(self):
        """Test riveted steel (very rough)."""
        eps = cb.pipe_roughness('riveted_steel')
        # Moody (1944), Crane (2009): 0.9 mm = 9.0e-4 m
        assert eps == 9.0e-4

    def test_case_insensitive(self):
        """Test case-insensitive lookup."""
        eps_lower = cb.pipe_roughness('commercial_steel')
        eps_upper = cb.pipe_roughness('COMMERCIAL_STEEL')
        eps_mixed = cb.pipe_roughness('Commercial_Steel')
        
        assert eps_lower == eps_upper
        assert eps_lower == eps_mixed
        assert eps_upper == eps_mixed

    def test_aliases(self):
        """Test material aliases."""
        # new_steel is alias for commercial_steel
        eps_commercial = cb.pipe_roughness('commercial_steel')
        eps_new = cb.pipe_roughness('new_steel')
        assert eps_commercial == eps_new
        
        # galvanized_steel is alias for galvanized_iron
        eps_iron = cb.pipe_roughness('galvanized_iron')
        eps_steel = cb.pipe_roughness('galvanized_steel')
        assert eps_iron == eps_steel

    def test_unknown_material_error(self):
        """Test error handling for unknown material."""
        with pytest.raises(ValueError, match="unknown material"):
            cb.pipe_roughness('nonexistent_material')
        
        with pytest.raises(ValueError, match="unknown material"):
            cb.pipe_roughness('aluminum')  # Not in database

    def test_standard_pipe_roughness_dict(self):
        """Test standard_pipe_roughness returns complete dictionary."""
        roughness_db = cb.standard_pipe_roughness()
        
        # Should be a dictionary
        assert isinstance(roughness_db, dict)
        
        # Should contain all expected materials
        expected_materials = [
            'smooth', 'drawn_tubing', 'pvc', 'plastic',
            'commercial_steel', 'new_steel', 'wrought_iron',
            'galvanized_iron', 'galvanized_steel', 'rusted_steel',
            'cast_iron', 'asphalted_cast_iron',
            'concrete', 'rough_concrete',
            'riveted_steel', 'wood_stave', 'corrugated_metal'
        ]
        
        for material in expected_materials:
            assert material in roughness_db
            assert roughness_db[material] >= 0.0

    def test_roughness_ordering(self):
        """Test that roughness values are in expected order."""
        # Smooth < drawn tubing < commercial steel < galvanized < cast iron < concrete < riveted
        eps_smooth = cb.pipe_roughness('smooth')
        eps_drawn = cb.pipe_roughness('drawn_tubing')
        eps_steel = cb.pipe_roughness('commercial_steel')
        eps_galv = cb.pipe_roughness('galvanized_iron')
        eps_cast = cb.pipe_roughness('cast_iron')
        eps_concrete = cb.pipe_roughness('concrete')
        eps_riveted = cb.pipe_roughness('riveted_steel')
        
        assert eps_smooth < eps_drawn
        assert eps_drawn < eps_steel
        assert eps_steel < eps_galv
        assert eps_galv < eps_cast
        assert eps_cast < eps_concrete
        assert eps_concrete < eps_riveted

    def test_realistic_values(self):
        """Test that all roughness values are realistic."""
        roughness_db = cb.standard_pipe_roughness()
        
        for material, eps in roughness_db.items():
            # All values should be non-negative
            assert eps >= 0.0
            
            # All values should be less than 10 cm (very rough)
            assert eps < 0.1
            
            # Most values should be in typical range (0 to 10 mm)
            if material != 'corrugated_metal':
                assert eps < 0.01

    def test_usage_with_friction_factor(self):
        """Test using pipe_roughness with friction factor calculation."""
        # Get roughness for commercial steel
        eps = cb.pipe_roughness('commercial_steel')
        D = 0.1  # 100 mm pipe
        e_D = eps / D  # Relative roughness
        
        # Calculate friction factor for turbulent flow
        Re = 1e5
        f = cb.friction_haaland(Re, e_D)
        
        # Should be reasonable value for turbulent flow
        assert 0.01 < f < 0.05
        
        # Smooth pipe should have lower friction
        eps_smooth = cb.pipe_roughness('smooth')
        e_D_smooth = eps_smooth / D
        f_smooth = cb.friction_haaland(Re, e_D_smooth)
        
        assert f_smooth < f

    def test_all_materials_accessible(self):
        """Test that all materials in database are accessible."""
        roughness_db = cb.standard_pipe_roughness()
        
        # Every material in the database should be retrievable
        for material in roughness_db.keys():
            eps_from_dict = roughness_db[material]
            eps_from_func = cb.pipe_roughness(material)
            assert eps_from_dict == eps_from_func

    def test_roughness_units(self):
        """Test that roughness values are in meters."""
        # Commercial steel: 0.045 mm = 4.5e-5 m
        eps = cb.pipe_roughness('commercial_steel')
        assert 4.0e-5 < eps < 5.0e-5  # Should be around 45 microns
        
        # Galvanized iron: 0.15 mm = 1.5e-4 m
        eps = cb.pipe_roughness('galvanized_iron')
        assert 1.4e-4 < eps < 1.6e-4  # Should be around 150 microns
