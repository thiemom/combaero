"""Central-difference accuracy checks for the pybind-exposed ejector C++
(f, J) functions (combaero._solver_tools.ejector_*_and_jacobian, i.e.
include/ejector.h / src/ejector.cpp).

Mirrors test_jacobian_migration.py: call each binding once for the analytic
result, then perturb each of the four thermodynamic inputs (p_g, t_g, p_e,
t_e -- plus p0/t0 for the choked-mass-flow function) and confirm the
analytic derivative fields match central finite differences. The C++ ctests
(tests/test_ejector_jacobians.cpp) already validate the same Jacobians
against the baked reference data; this is the Python-binding-surface guard,
proving the pybind wiring and struct field names round-trip correctly.
"""

import pytest

import combaero as cb

# A representative critical-mode operating point (Huang 1999 case T3-EH),
# air-like gamma/r_gas stand-ins are unnecessary -- these are the reference
# R141b numbers, and the physics is input-parameterized on gamma/r_gas.
P_G = 604000.0
T_G = 368.15
P_E = 40000.0
T_E = 281.15
ARN = 3.271
ARM = 10.87
GAMMA = 1.1121
R_GAS = 71.094
RECOVERY = 1.0
THROAT_AREA = 3.14e-5
ETA = 0.95

REL = 1e-4


def _cd(f, x0, i, eps_rel=1e-6):
    """Central difference of scalar f(vec) w.r.t. component i of x0."""
    eps = eps_rel * abs(x0[i]) if x0[i] != 0.0 else eps_rel
    xp = list(x0)
    xm = list(x0)
    xp[i] += eps
    xm[i] -= eps
    return (f(xp) - f(xm)) / (2.0 * eps)


def test_choked_mass_flow_jacobian_accuracy():
    st = cb._solver_tools
    x0 = [P_G, T_G]  # (p0, t0)

    def mdot(x):
        return st.ejector_choked_mass_flow_and_jacobian(x[0], x[1], THROAT_AREA, GAMMA, R_GAS, ETA)[
            0
        ]

    val, d_dp0, d_dt0 = st.ejector_choked_mass_flow_and_jacobian(
        P_G, T_G, THROAT_AREA, GAMMA, R_GAS, ETA
    )
    assert val > 0.0
    assert d_dp0 == pytest.approx(_cd(mdot, x0, 0), rel=REL)
    assert d_dt0 == pytest.approx(_cd(mdot, x0, 1), rel=REL)


def test_entrainment_ratio_jacobian_accuracy():
    st = cb._solver_tools
    x0 = [P_G, T_G, P_E, T_E]

    def omega(x):
        return st.ejector_entrainment_ratio_and_jacobian(
            x[0], x[1], x[2], x[3], ARN, ARM, GAMMA
        ).value.omega

    res = st.ejector_entrainment_ratio_and_jacobian(P_G, T_G, P_E, T_E, ARN, ARM, GAMMA)
    assert res.domega_dp_g == pytest.approx(_cd(omega, x0, 0), rel=REL)
    assert res.domega_dt_g == pytest.approx(_cd(omega, x0, 1), rel=REL)
    assert res.domega_dp_e == pytest.approx(_cd(omega, x0, 2), rel=REL)
    assert res.domega_dt_e == pytest.approx(_cd(omega, x0, 3), rel=REL)


def test_critical_back_pressure_jacobian_accuracy():
    st = cb._solver_tools
    x0 = [P_G, T_G, P_E, T_E]

    def pc(x):
        return st.ejector_critical_back_pressure_and_jacobian(
            x[0], x[1], x[2], x[3], ARN, ARM, GAMMA, R_GAS, RECOVERY
        ).value.p_c_pa

    res = st.ejector_critical_back_pressure_and_jacobian(
        P_G, T_G, P_E, T_E, ARN, ARM, GAMMA, R_GAS, RECOVERY
    )
    assert res.dpc_dp_g == pytest.approx(_cd(pc, x0, 0), rel=REL)
    assert res.dpc_dt_g == pytest.approx(_cd(pc, x0, 1), rel=REL)
    assert res.dpc_dp_e == pytest.approx(_cd(pc, x0, 2), rel=REL)
    assert res.dpc_dt_e == pytest.approx(_cd(pc, x0, 3), rel=REL)


def test_recovery_efficiency_scales_p_c_linearly():
    st = cb._solver_tools
    base = st.ejector_critical_back_pressure_and_jacobian(
        P_G, T_G, P_E, T_E, ARN, ARM, GAMMA, R_GAS, 1.0
    ).value
    half = st.ejector_critical_back_pressure_and_jacobian(
        P_G, T_G, P_E, T_E, ARN, ARM, GAMMA, R_GAS, 0.5
    ).value
    # P_c* = recovery_efficiency * p03; p03 itself is recovery-independent.
    assert half.p_c_pa == pytest.approx(0.5 * base.p_c_pa, rel=1e-12)
    assert half.p_mixed_stagnation_pa == pytest.approx(base.p_mixed_stagnation_pa, rel=1e-12)
