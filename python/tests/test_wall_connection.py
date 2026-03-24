"""Test WallConnection dataclass."""

import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel, WallConnection


def test_wall_connection_basic():
    """Test basic WallConnection instantiation and compute_coupling."""
    wall = WallConnection(
        id="test_wall",
        element_a="hot_pipe",
        element_b="cold_pipe",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )

    assert wall.id == "test_wall"
    assert wall.element_a == "hot_pipe"
    assert wall.element_b == "cold_pipe"
    assert wall.wall_thickness == 0.002
    assert wall.wall_conductivity == 25.0
    assert wall.contact_area is None


def test_wall_connection_compute_coupling():
    """Test WallConnection.compute_coupling method."""
    wall = WallConnection(
        id="test_wall",
        element_a="hot_pipe",
        element_b="cold_pipe",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )

    # Hot side: h=500 W/(m^2*K), T_aw=1000 K, A=0.1 m^2
    # Cold side: h=2000 W/(m^2*K), T_aw=400 K, A=0.15 m^2
    Q, T_wall = wall.compute_coupling(
        h_a=500.0,
        T_aw_a=1000.0,
        A_conv_a=0.1,
        h_b=2000.0,
        T_aw_b=400.0,
        A_conv_b=0.15,
    )

    # Should use min(A_conv_a, A_conv_b) = 0.1 m^2
    assert Q > 0.0  # Heat flows from hot to cold
    assert T_wall > 400.0 and T_wall < 1000.0

    # Verify against direct C++ call
    t_over_k = 0.002 / 25.0
    result = cb.wall_coupling_and_jacobian(500.0, 1000.0, 2000.0, 400.0, t_over_k, 0.1)
    assert abs(Q - result.Q) < 1e-10
    assert abs(T_wall - result.T_wall) < 1e-10


def test_wall_connection_with_contact_area():
    """Test WallConnection with explicit contact_area override."""
    wall = WallConnection(
        id="test_wall",
        element_a="hot_pipe",
        element_b="cold_pipe",
        wall_thickness=0.002,
        wall_conductivity=25.0,
        contact_area=0.05,  # Override area
    )

    Q, T_wall = wall.compute_coupling(
        h_a=500.0,
        T_aw_a=1000.0,
        A_conv_a=0.1,
        h_b=2000.0,
        T_aw_b=400.0,
        A_conv_b=0.15,
    )

    # Should use contact_area = 0.05 m^2 instead of min(0.1, 0.15)
    t_over_k = 0.002 / 25.0
    result = cb.wall_coupling_and_jacobian(500.0, 1000.0, 2000.0, 400.0, t_over_k, 0.05)
    assert abs(Q - result.Q) < 1e-10


def test_wall_connection_with_convective_surfaces():
    """Test WallConnection with ConvectiveSurface objects."""
    # Create two convective surfaces
    hot_surface = ConvectiveSurface(
        area=0.1, model=SmoothModel(correlation="gnielinski"), heating=False
    )

    cold_surface = ConvectiveSurface(
        area=0.12, model=SmoothModel(correlation="gnielinski"), heating=True
    )

    # Compute HTCs from flow conditions
    X_air = cb.species.dry_air()

    # Hot side
    ch_hot = hot_surface.htc_and_T(
        T=1000.0,
        P=2e5,
        X=X_air,
        velocity=50.0,
        diameter=0.02,
        length=0.5,
        T_wall=800.0,
    )

    # Cold side
    ch_cold = cold_surface.htc_and_T(
        T=400.0,
        P=2.5e5,
        X=X_air,
        velocity=60.0,
        diameter=0.015,
        length=0.5,
        T_wall=500.0,
    )

    # Create wall connection
    wall = WallConnection(
        id="test_wall",
        element_a="hot_pipe",
        element_b="cold_pipe",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )

    # Compute coupling
    Q, T_wall = wall.compute_coupling(
        ch_hot.h,
        ch_hot.T_aw,
        hot_surface.area,
        ch_cold.h,
        ch_cold.T_aw,
        cold_surface.area,
    )

    assert Q > 0.0  # Heat flows from hot to cold
    assert T_wall > ch_cold.T_aw and T_wall < ch_hot.T_aw
    assert ch_hot.h > 0.0 and ch_cold.h > 0.0


def test_wall_connection_zero_temperature_difference():
    """Test WallConnection with equal temperatures on both sides."""
    wall = WallConnection(
        id="test_wall",
        element_a="pipe_a",
        element_b="pipe_b",
        wall_thickness=0.002,
        wall_conductivity=25.0,
    )

    Q, T_wall = wall.compute_coupling(
        h_a=500.0,
        T_aw_a=700.0,
        A_conv_a=0.1,
        h_b=2000.0,
        T_aw_b=700.0,  # Same temperature
        A_conv_b=0.1,
    )

    assert abs(Q) < 1e-10  # No heat transfer
    assert abs(T_wall - 700.0) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
