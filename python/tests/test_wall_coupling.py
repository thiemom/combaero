"""Test wall_coupling_and_jacobian function and WallCouplingResult."""

import pytest

import combaero as cb


def test_hot_coupling_basic():
    """Test basic wall coupling calculation."""
    # Hot side: h=500 W/(m^2*K), T_aw=1000 K
    # Cold side: h=2000 W/(m^2*K), T_aw=400 K
    # Wall: t/k = 0.002/25 = 8e-5 m^2*K/W
    # Area: 0.1 m^2

    h_a = 500.0
    T_aw_a = 1000.0
    h_b = 2000.0
    T_aw_b = 400.0
    t_over_k = 0.002 / 25.0  # 8e-5
    A = 0.1

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    # Overall HTC: U = 1 / (1/500 + 8e-5 + 1/2000) = 1 / 0.00258 = 387.6 W/(m^2*K)
    U_expected = 1.0 / (1.0 / h_a + t_over_k + 1.0 / h_b)
    Q_expected = U_expected * A * (T_aw_a - T_aw_b)

    assert abs(result.Q - Q_expected) < 1e-6
    assert result.Q > 0.0  # Heat flows from hot to cold
    assert result.T_hot > T_aw_b and result.T_hot < T_aw_a


def test_hot_coupling_jacobians():
    """Test wall coupling Jacobians with finite differences."""
    h_a = 500.0
    T_aw_a = 1000.0
    h_b = 2000.0
    T_aw_b = 400.0
    t_over_k = 8e-5
    A = 0.1

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    # Finite difference validation
    eps = 1e-6

    # dQ/dh_a
    result_plus = cb.wall_coupling_and_jacobian(h_a + eps, T_aw_a, h_b, T_aw_b, t_over_k, A)
    result_minus = cb.wall_coupling_and_jacobian(h_a - eps, T_aw_a, h_b, T_aw_b, t_over_k, A)
    dQ_dh_a_fd = (result_plus.Q - result_minus.Q) / (2 * eps)
    assert abs(result.dQ_dh_a - dQ_dh_a_fd) / abs(dQ_dh_a_fd) < 1e-5

    # dQ/dh_b
    result_plus = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b + eps, T_aw_b, t_over_k, A)
    result_minus = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b - eps, T_aw_b, t_over_k, A)
    dQ_dh_b_fd = (result_plus.Q - result_minus.Q) / (2 * eps)
    assert abs(result.dQ_dh_b - dQ_dh_b_fd) / abs(dQ_dh_b_fd) < 1e-5

    # dQ/dT_aw_a
    result_plus = cb.wall_coupling_and_jacobian(h_a, T_aw_a + eps, h_b, T_aw_b, t_over_k, A)
    result_minus = cb.wall_coupling_and_jacobian(h_a, T_aw_a - eps, h_b, T_aw_b, t_over_k, A)
    dQ_dT_aw_a_fd = (result_plus.Q - result_minus.Q) / (2 * eps)
    assert abs(result.dQ_dT_aw_a - dQ_dT_aw_a_fd) / abs(dQ_dT_aw_a_fd) < 1e-5

    # dQ/dT_aw_b
    result_plus = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b + eps, t_over_k, A)
    result_minus = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b - eps, t_over_k, A)
    dQ_dT_aw_b_fd = (result_plus.Q - result_minus.Q) / (2 * eps)
    assert abs(result.dQ_dT_aw_b - dQ_dT_aw_b_fd) / abs(dQ_dT_aw_b_fd) < 1e-5


def test_hot_coupling_zero_temperature_difference():
    """Test wall coupling with zero temperature difference."""
    h_a = 500.0
    T_aw_a = 700.0
    h_b = 2000.0
    T_aw_b = 700.0  # Same temperature
    t_over_k = 8e-5
    A = 0.1

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    assert abs(result.Q) < 1e-10  # No heat transfer
    assert abs(result.T_hot - T_aw_a) < 1e-10  # Wall at same temperature


def test_hot_coupling_high_resistance():
    """Test wall coupling with high wall resistance."""
    h_a = 500.0
    T_aw_a = 1000.0
    h_b = 2000.0
    T_aw_b = 400.0
    t_over_k = 0.01  # High resistance wall
    A = 0.1

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    # With high wall resistance, U should be dominated by wall resistance
    U_expected = 1.0 / (1.0 / h_a + t_over_k + 1.0 / h_b)
    assert U_expected < 100.0  # Low overall HTC due to wall resistance
    assert result.Q > 0.0


def test_hot_coupling_asymmetric_htcs():
    """Test wall coupling with very different HTCs on each side."""
    h_a = 50.0  # Low HTC (gas side)
    T_aw_a = 1200.0
    h_b = 5000.0  # High HTC (liquid side)
    T_aw_b = 300.0
    t_over_k = 1e-4
    A = 0.05

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    # Overall HTC should be dominated by the lower HTC (gas side)
    U_expected = 1.0 / (1.0 / h_a + t_over_k + 1.0 / h_b)
    assert U_expected < h_a * 1.1  # U should be close to h_a

    # T_hot is the hot-side surface temperature: T_aw_a - Q / (h_a * A)
    T_hot_expected = T_aw_a - result.Q / (h_a * A)
    assert abs(result.T_hot - T_hot_expected) < 1e-6
    assert result.T_hot < T_aw_a  # Below the hot-side adiabatic wall temperature


def test_hot_coupling_repr():
    """Test WallCouplingResult string representation."""
    h_a = 500.0
    T_aw_a = 1000.0
    h_b = 2000.0
    T_aw_b = 400.0
    t_over_k = 8e-5
    A = 0.1

    result = cb.wall_coupling_and_jacobian(h_a, T_aw_a, h_b, T_aw_b, t_over_k, A)

    repr_str = repr(result)
    assert "WallCouplingResult" in repr_str
    assert "Q=" in repr_str
    assert "T_hot=" in repr_str


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
