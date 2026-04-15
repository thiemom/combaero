"""Tests for combustor pressure-loss correlation classes."""

from __future__ import annotations

import pytest

from combaero.network import (
    ConstantFractionLoss,
    ConstantLoss,
    LinearThetaFractionLoss,
    LinearThetaLoss,
)

# ---------------------------------------------------------------------------
# Minimal duck-typed context for pure-Python tests (no C++ required)
# ---------------------------------------------------------------------------


class _Ctx:
    def __init__(self, theta: float = 0.0) -> None:
        self.theta = theta


# ---------------------------------------------------------------------------
# ConstantFractionLoss (and alias ConstantLoss)
# ---------------------------------------------------------------------------


class TestConstantFractionLoss:
    def test_returns_xi_and_zero_derivative(self) -> None:
        loss = ConstantFractionLoss(xi=0.03)
        xi, dxi_dtheta = loss(_Ctx(theta=0.0))
        assert xi == pytest.approx(0.03)
        assert dxi_dtheta == pytest.approx(0.0)

    def test_independent_of_theta(self) -> None:
        loss = ConstantFractionLoss(xi=0.05)
        for theta in [-1.0, 0.0, 0.5, 2.0, 10.0]:
            xi, dxi = loss(_Ctx(theta=theta))
            assert xi == pytest.approx(0.05)
            assert dxi == pytest.approx(0.0)

    def test_default_xi(self) -> None:
        loss = ConstantFractionLoss()
        xi, _ = loss(_Ctx())
        assert xi == pytest.approx(0.03)

    def test_zero_loss(self) -> None:
        loss = ConstantFractionLoss(xi=0.0)
        xi, dxi = loss(_Ctx(theta=1.5))
        assert xi == pytest.approx(0.0)
        assert dxi == pytest.approx(0.0)

    def test_alias_constant_loss(self) -> None:
        assert ConstantLoss is ConstantFractionLoss


# ---------------------------------------------------------------------------
# LinearThetaFractionLoss (and alias LinearThetaLoss)
# ---------------------------------------------------------------------------


class TestLinearThetaFractionLoss:
    def test_at_zero_theta(self) -> None:
        loss = LinearThetaFractionLoss(k=0.5, xi0=0.02)
        xi, dxi = loss(_Ctx(theta=0.0))
        assert xi == pytest.approx(0.02)
        assert dxi == pytest.approx(0.5)

    def test_linear_slope(self) -> None:
        loss = LinearThetaFractionLoss(k=0.4, xi0=0.01)
        theta = 2.0
        xi, dxi = loss(_Ctx(theta=theta))
        assert xi == pytest.approx(0.4 * theta + 0.01)
        assert dxi == pytest.approx(0.4)

    def test_derivative_constant_wrt_theta(self) -> None:
        loss = LinearThetaFractionLoss(k=0.6, xi0=0.03)
        for theta in [0.0, 0.5, 1.0, 3.0]:
            _, dxi = loss(_Ctx(theta=theta))
            assert dxi == pytest.approx(0.6)

    def test_default_values(self) -> None:
        loss = LinearThetaFractionLoss()
        xi, dxi = loss(_Ctx(theta=1.0))
        assert xi == pytest.approx(0.5 * 1.0 + 0.02)
        assert dxi == pytest.approx(0.5)

    def test_fd_consistency(self) -> None:
        """Analytical derivative matches finite difference."""
        loss = LinearThetaFractionLoss(k=0.35, xi0=0.015)
        theta = 1.5
        eps = 1e-6
        xi0_val, dxi = loss(_Ctx(theta=theta))
        xi_p, _ = loss(_Ctx(theta=theta + eps))
        fd = (xi_p - xi0_val) / eps
        assert dxi == pytest.approx(fd, rel=1e-5)

    def test_alias_linear_theta_loss(self) -> None:
        assert LinearThetaLoss is LinearThetaFractionLoss


# ---------------------------------------------------------------------------
# Integration test: CombustorNode + LinearThetaFractionLoss reduces P_total
# ---------------------------------------------------------------------------


class TestCombustorNodeWithPressureLoss:
    def _make_streams(self):
        import combaero as cb
        from combaero.network import MixtureState

        ns = cb.species.num_species
        fuel_Y = [0.0] * ns
        air_Y = [0.0] * ns
        for k in range(ns):
            name = cb.species_name(k)
            if name == "CH4":
                fuel_Y[k] = 1.0
            elif name == "N2":
                air_Y[k] = 0.767
            elif name == "O2":
                air_Y[k] = 0.233

        fuel = MixtureState(
            P=150000.0, P_total=150000.0, T=300.0, T_total=300.0, m_dot=0.05, Y=fuel_Y
        )
        air = MixtureState(P=150000.0, P_total=150000.0, T=600.0, T_total=600.0, m_dot=1.0, Y=air_Y)
        return fuel, air

    def test_no_loss_baseline(self) -> None:
        from combaero.network import CombustorNode

        fuel, air = self._make_streams()
        node = CombustorNode("comb", method="complete")
        _, _, mix_res = node.compute_derived_state([fuel, air])
        assert mix_res is not None
        assert mix_res.P_total_mix == pytest.approx(150000.0, rel=1e-6)

    def test_constant_loss_reduces_P_total(self) -> None:
        from combaero.network import CombustorNode

        fuel, air = self._make_streams()
        node = CombustorNode("comb", method="complete", pressure_loss=ConstantFractionLoss(xi=0.05))
        _, _, mix_res = node.compute_derived_state([fuel, air])
        assert mix_res.P_total_mix == pytest.approx(150000.0 * (1.0 - 0.05), rel=1e-4)

    def test_linear_theta_loss_reduces_P_total(self) -> None:
        from combaero.network import CombustorNode

        fuel, air = self._make_streams()
        node_base = CombustorNode("base", method="complete")
        _, _, mix_base = node_base.compute_derived_state([fuel, air])

        node_loss = CombustorNode(
            "loss", method="complete", pressure_loss=LinearThetaFractionLoss(k=0.5, xi0=0.02)
        )
        _, _, mix_loss = node_loss.compute_derived_state([fuel, air])

        assert mix_loss.P_total_mix < 150000.0
        assert mix_loss.P_total_mix < mix_base.P_total_mix

    def test_outlet_T_unchanged_by_pressure_loss(self) -> None:
        """Pressure loss does not affect outlet temperature (ideal gas)."""
        from combaero.network import CombustorNode

        fuel, air = self._make_streams()
        node_base = CombustorNode("base", method="complete")
        T_base, _, _ = node_base.compute_derived_state([fuel, air])

        node_loss = CombustorNode(
            "loss", method="complete", pressure_loss=LinearThetaFractionLoss(k=0.5, xi0=0.02)
        )
        T_loss, _, _ = node_loss.compute_derived_state([fuel, air])

        assert T_loss == pytest.approx(T_base, rel=1e-6)

    def test_constant_head_loss_reduces_P_total(self) -> None:
        from combaero.network import CombustorNode, ConstantHeadLoss

        fuel, air = self._make_streams()
        node = CombustorNode(
            "head", method="complete", pressure_loss=ConstantHeadLoss(zeta=5.0, area=0.1)
        )
        _, _, mix_res = node.compute_derived_state([fuel, air])
        assert mix_res.P_total_mix < 150000.0

    def test_head_loss_smaller_than_full_p_in(self) -> None:
        """Head loss with realistic zeta should not exceed P_in."""
        from combaero.network import CombustorNode, ConstantHeadLoss

        fuel, air = self._make_streams()
        node = CombustorNode(
            "head", method="complete", pressure_loss=ConstantHeadLoss(zeta=10.0, area=0.05)
        )
        _, _, mix_res = node.compute_derived_state([fuel, air])
        assert mix_res.P_total_mix > 0.0
