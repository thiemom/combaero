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
# Integration test: PressureLossElement residual enforces correct ratio
# ---------------------------------------------------------------------------


class TestPressureLossElementResidual:
    """Unit-level checks on ``PressureLossElement.residuals`` for each correlation.

    Verifies the element residual ``P_total_out - P_total_in * (1 - xi) = 0``
    vanishes at the physically correct P ratio, for every supported correlation.
    """

    def _make_air_state(self, P_total: float = 150000.0, T: float = 700.0, m_dot: float = 1.0):
        """Helper: build an air-composition MixtureState."""
        import combaero as cb
        from combaero.network import MixtureState

        ns = cb.species.num_species
        Y = [0.0] * ns
        for k in range(ns):
            name = cb.species_name(k)
            if name == "N2":
                Y[k] = 0.767
            elif name == "O2":
                Y[k] = 0.233
        return MixtureState(P=P_total, P_total=P_total, T=T, T_total=T, m_dot=m_dot, Y=Y)

    def test_constant_fraction_residual_zero_at_exact_ratio(self) -> None:
        """res == 0 iff P_total_out / P_total_in == (1 - xi)."""
        from combaero.network import PressureLossElement

        xi = 0.05
        P_in = 150000.0
        state_in = self._make_air_state(P_total=P_in)
        state_out = self._make_air_state(P_total=P_in * (1.0 - xi))

        elem = PressureLossElement("loss", "a", "b", correlation=ConstantFractionLoss(xi=xi))
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)

    def test_constant_fraction_residual_nonzero_off_ratio(self) -> None:
        from combaero.network import PressureLossElement

        xi = 0.05
        P_in = 150000.0
        state_in = self._make_air_state(P_total=P_in)
        state_out = self._make_air_state(P_total=P_in)  # no drop

        elem = PressureLossElement("loss", "a", "b", correlation=ConstantFractionLoss(xi=xi))
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(P_in * xi, rel=1e-6)

    def test_linear_theta_residual_uses_cold_flow_without_source(self) -> None:
        """Without a theta source, linear-theta correlation degrades to xi0."""
        import warnings as _w

        from combaero.network import PressureLossElement

        k, xi0 = 0.5, 0.02
        P_in = 150000.0
        state_in = self._make_air_state(P_total=P_in)
        state_out = self._make_air_state(P_total=P_in * (1.0 - xi0))

        elem = PressureLossElement(
            "loss", "a", "b", correlation=LinearThetaFractionLoss(k=k, xi0=xi0)
        )
        # residuals() itself does not warn (warning fires in resolve_topology);
        # here we bypass resolve_topology and hit the no-graph path.
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)

    def test_constant_head_residual_xi_from_dynamic_head(self) -> None:
        """ConstantHeadLoss: xi = zeta * 0.5 * rho * v^2 / P_in."""
        from combaero.network import ConstantHeadLoss, PressureLossElement

        zeta = 5.0
        area = 0.1
        P_in = 150000.0
        mdot = 1.0
        state_in = self._make_air_state(P_total=P_in, m_dot=mdot)

        corr = ConstantHeadLoss(zeta=zeta, area=area)
        elem = PressureLossElement("loss", "a", "b", correlation=corr)
        # Evaluate xi manually via the same ctx builder used by the element
        ctx = elem._build_ctx(state_in, state_in, None)
        xi, _ = corr(ctx)
        assert 0.0 < xi < 1.0

        state_out = self._make_air_state(P_total=P_in * (1.0 - xi), m_dot=mdot)
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)
