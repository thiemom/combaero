"""Tests for acoustic mode analysis."""

import pytest

import combaero as cb


class TestTransferMatrixMethod:
    """Test advanced thermoacoustic transfer matrix tools."""

    def test_orifice_impedance_baseline(self):
        """Test orifice impedance at typical conditions."""
        Z = cb.orifice_impedance_with_flow(
            freq=100,
            u_bias=10,
            u_grazing=5,
            d_orifice=0.01,
            l_orifice=0.005,
            porosity=0.1,
            Cd=0.61,
            rho=1.2,
            c=343,
        )
        # Should have positive real (resistance) and imaginary (reactance) parts
        assert Z.real > 0
        assert Z.imag > 0

    def test_orifice_impedance_no_flow(self):
        """Test that zero flow gives only reactance."""
        Z = cb.orifice_impedance_with_flow(
            freq=100,
            u_bias=0,
            u_grazing=0,
            d_orifice=0.01,
            l_orifice=0.005,
            porosity=0.1,
            Cd=0.61,
            rho=1.2,
            c=343,
        )
        # No flow -> no resistance
        assert Z.real == pytest.approx(0, abs=1e-10)
        assert Z.imag > 0  # Still has inertance

    def test_quarter_wave_tmm_structure(self):
        """Test quarter-wave resonator returns valid transfer matrix."""
        tm = cb.quarter_wave_resonator_tmm(
            freq=100,
            L_tube=0.85,
            A_duct=0.01,
            A_tube=0.005,
            d_orifice=0.01,
            l_orifice=0.005,
            porosity=0.1,
            Cd=0.61,
            u_bias=10,
            u_grazing=5,
            rho=1.2,
            c=343,
        )
        # Should be a shunt element: T11=1, T12=0, T22=1
        assert pytest.approx(1.0 + 0j, abs=1e-10) == tm.T11
        assert pytest.approx(0.0 + 0j, abs=1e-10) == tm.T12
        assert pytest.approx(1.0 + 0j, abs=1e-10) == tm.T22
        # T21 should be non-zero (branch admittance)
        assert abs(tm.T21) > 0

    def test_transfer_matrix_multiplication(self):
        """Test transfer matrix multiplication."""
        tm1 = cb.TransferMatrix()
        tm1.T11 = 1.0 + 0j
        tm1.T12 = 0.5 + 0j
        tm1.T21 = 0.0 + 0j
        tm1.T22 = 1.0 + 0j

        tm2 = cb.TransferMatrix()
        tm2.T11 = 1.0 + 0j
        tm2.T12 = 0.0 + 0j
        tm2.T21 = 2.0 + 0j
        tm2.T22 = 1.0 + 0j

        tm_result = tm1 * tm2
        # Verify matrix multiplication: [1 0.5] * [1 0] = [1+0.5*2  0+0.5*1] = [2  0.5]
        #                                [0  1 ]   [2 1]   [0+1*2    0+1*1  ]   [2  1  ]
        assert pytest.approx(2.0 + 0j) == tm_result.T11
        assert pytest.approx(0.5 + 0j) == tm_result.T12
        assert pytest.approx(2.0 + 0j) == tm_result.T21
        assert pytest.approx(1.0 + 0j) == tm_result.T22

    def test_whistling_risk_detection(self):
        """Test Strouhal-based whistling risk."""
        # St = f*d/U = 500*0.01/50 = 0.1 (safe)
        assert cb.is_whistling_risk(500, 50, 0.01) is False

        # St = 1000*0.01/50 = 0.2 (critical)
        assert cb.is_whistling_risk(1000, 50, 0.01) is True

        # St = 2000*0.01/50 = 0.4 (critical)
        assert cb.is_whistling_risk(2000, 50, 0.01) is True

        # St = 3000*0.01/50 = 0.6 (safe)
        assert cb.is_whistling_risk(3000, 50, 0.01) is False

    def test_whistling_no_flow(self):
        """Test that zero flow means no whistling risk."""
        assert cb.is_whistling_risk(1000, 0, 0.01) is False

    def test_orifice_impedance_parameter_validation(self):
        """Test parameter validation for orifice impedance."""
        # Mach number too high
        with pytest.raises(RuntimeError, match="Mach"):
            cb.orifice_impedance_with_flow(100, 150, 5, 0.01, 0.005, 0.1, 0.61, 1.2, 343)

        # Invalid porosity
        with pytest.raises(RuntimeError, match="Porosity"):
            cb.orifice_impedance_with_flow(100, 10, 5, 0.01, 0.005, 0.6, 0.61, 1.2, 343)

    def test_qw_tmm_cutoff_validation(self):
        """Test that quarter-wave TMM rejects frequencies above cutoff."""
        # kL > pi should fail
        with pytest.raises(RuntimeError, match="cutoff"):
            cb.quarter_wave_resonator_tmm(
                freq=1000,
                L_tube=0.85,
                A_duct=0.01,
                A_tube=0.005,
                d_orifice=0.01,
                l_orifice=0.005,
                porosity=0.1,
                Cd=0.61,
                u_bias=10,
                u_grazing=5,
                rho=1.2,
                c=343,
            )


class TestLinerHelperApi:
    """Test streamlined liner helper APIs for demos and user code."""

    @staticmethod
    def _build_common_objects():
        medium = cb.AcousticMedium()
        medium.rho = 1.2
        medium.c = 340.0

        orifice = cb.LinerOrificeGeometry()
        orifice.d_orifice = 0.003
        orifice.l_orifice = 0.005
        orifice.porosity = 0.08
        orifice.Cd = 0.7

        flow = cb.LinerFlowState()
        flow.u_bias = 10.0
        flow.u_grazing = 80.0

        cavity = cb.LinerCavity()
        cavity.depth = 0.05

        return medium, orifice, flow, cavity

    def test_sdof_absorption_range(self):
        """SDOF helper should return bounded physical absorption."""
        medium, orifice, flow, cavity = self._build_common_objects()
        alpha = cb.liner_sdof_absorption(400.0, orifice, cavity, flow, medium)
        assert 0.0 <= alpha <= 1.0

    def test_sdof_sweep_length(self):
        """Sweep helper returns one value per frequency."""
        medium, orifice, flow, cavity = self._build_common_objects()
        freqs = [200.0, 400.0, 600.0]
        alpha = cb.sweep_liner_sdof_absorption(freqs, orifice, cavity, flow, medium)
        assert len(alpha) == len(freqs)

    def test_serial_2dof_differs_from_sdof(self):
        """2-DOF serial response should differ from baseline SDOF response."""
        medium, orifice, flow, cavity = self._build_common_objects()

        alpha_sdof = cb.sweep_liner_sdof_absorption(
            [300.0, 600.0, 900.0], orifice, cavity, flow, medium
        )
        alpha_2dof = cb.sweep_liner_2dof_serial_absorption(
            [300.0, 600.0, 900.0],
            orifice,
            orifice,
            0.025,
            0.025,
            flow,
            flow,
            medium,
        )

        assert any(abs(a - b) > 1e-6 for a, b in zip(alpha_sdof, alpha_2dof, strict=True))
