"""Acoustics checked against closed forms derived independently of the code.

Each reference here is either a textbook closed form or a convention-free
reduction (a symmetry argument, a 1-D finite-element model of the full
N-can graph), so a sign or factor error in the implementation cannot be
reproduced by the reference.
"""

import math

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.optimize import brentq

import combaero as cb

# ---------------------------------------------------------------------------
# Can-annular Bloch modes
# ---------------------------------------------------------------------------

L_CAN = 0.5
A_CAN = 0.01
R_PLENUM = 0.5
A_PLENUM = 0.02
C_CAN, C_PLENUM = 500.0, 600.0
RHO_CAN, RHO_PLENUM = 1.0, 1.2


def _geom(n_cans: int) -> cb.CanAnnularGeometry:
    g = cb.CanAnnularGeometry()
    g.n_cans = n_cans
    g.length_can = L_CAN
    g.area_can = A_CAN
    g.radius_plenum = R_PLENUM
    g.area_plenum = A_PLENUM
    return g


def _symmetric_residual(f: float, n_cans: int, push_pull: bool) -> float:
    # Volume-flow balance at one junction, normalised by i.
    # m = 0: by symmetry u' = 0 midway between cans -> each annulus arm of
    #        length s/2 is rigid-terminated, the same kind of branch as the
    #        closed-top can, so both admittances carry the same sign.
    # m = N/2: p' = 0 midway between cans -> each arm is pressure-released.
    w = 2.0 * math.pi * f
    s = 2.0 * math.pi * R_PLENUM / n_cans
    y_can = A_CAN / (RHO_CAN * C_CAN) * math.tan(w / C_CAN * L_CAN)
    y_arm = A_PLENUM / (RHO_PLENUM * C_PLENUM)
    half = w / C_PLENUM * s / 2.0
    if push_pull:
        return y_can - 2.0 * y_arm / math.tan(half)
    return y_can + 2.0 * y_arm * math.tan(half)


def _symmetric_roots(n_cans: int, push_pull: bool, f_max: float) -> list[float]:
    fs = np.linspace(5.0, f_max, 200001)
    v = np.array([_symmetric_residual(f, n_cans, push_pull) for f in fs])
    roots = []
    for i in range(len(fs) - 1):
        # sign change through zero, not through a tan() pole
        if np.sign(v[i]) != np.sign(v[i + 1]) and max(abs(v[i]), abs(v[i + 1])) < 1.0:
            roots.append(brentq(_symmetric_residual, fs[i], fs[i + 1], args=(n_cans, push_pull)))
    return roots


def _fem_frequencies(n_cans: int, n_sec: int = 24, n_can: int = 60) -> np.ndarray:
    # 1-D linear elements on the whole ring + N closed-top cans, pressure form
    # K p = w^2 M p. Rigid ends are natural BCs; no admittance sign enters.
    s = 2.0 * math.pi * R_PLENUM / n_cans
    n_ring = n_cans * n_sec
    n_tot = n_ring + n_cans * n_can
    K = np.zeros((n_tot, n_tot))
    M = np.zeros((n_tot, n_tot))

    def add(i: int, j: int, area: float, rho: float, c: float, h: float) -> None:
        k = area / rho / h
        m = area * h / (rho * c * c) / 6.0
        K[i, i] += k
        K[j, j] += k
        K[i, j] -= k
        K[j, i] -= k
        M[i, i] += 2.0 * m
        M[j, j] += 2.0 * m
        M[i, j] += m
        M[j, i] += m

    for i in range(n_ring):
        add(i, (i + 1) % n_ring, A_PLENUM, RHO_PLENUM, C_PLENUM, s / n_sec)
    nxt = n_ring
    for j in range(n_cans):
        prev = j * n_sec
        for _ in range(n_can):
            add(prev, nxt, A_CAN, RHO_CAN, C_CAN, L_CAN / n_can)
            prev = nxt
            nxt += 1
    w2 = eigh(K, M, eigvals_only=True)
    return np.sqrt(np.clip(w2, 0.0, None)) / (2.0 * math.pi)


class TestCanAnnularAgainstIndependentReferences:
    F_MAX = 1400.0

    def test_push_push_matches_symmetry_reduction(self) -> None:
        n = 24
        ref = _symmetric_roots(n, push_pull=False, f_max=self.F_MAX)
        modes = cb.can_annular_eigenmodes(
            _geom(n), C_CAN, C_PLENUM, RHO_CAN, RHO_PLENUM, f_max=self.F_MAX
        )
        got = sorted(m.frequency for m in modes if m.m_azimuthal == 0)
        assert len(ref) >= 2
        assert len(got) == len(ref), f"m=0 code {got} vs exact {ref}"
        for g, r in zip(got, ref, strict=True):
            assert g == pytest.approx(r, rel=2e-3)

    def test_push_pull_matches_symmetry_reduction(self) -> None:
        n = 24
        ref = _symmetric_roots(n, push_pull=True, f_max=self.F_MAX)
        modes = cb.can_annular_eigenmodes(
            _geom(n), C_CAN, C_PLENUM, RHO_CAN, RHO_PLENUM, f_max=self.F_MAX
        )
        got = sorted(m.frequency for m in modes if m.m_azimuthal == n // 2)
        assert len(ref) >= 2
        assert len(got) == len(ref), f"m=N/2 code {got} vs exact {ref}"
        for g, r in zip(got, ref, strict=True):
            assert g == pytest.approx(r, rel=2e-3)

    def test_every_mode_matches_full_network_fem(self) -> None:
        # All Bloch orders, not just the two symmetric ones.
        n = 8
        f_max = 900.0
        # FEM: linear elements, (k h)^2 / 24 < 4e-4 here.
        fem_all = _fem_frequencies(n)
        modes = cb.can_annular_eigenmodes(
            _geom(n), C_CAN, C_PLENUM, RHO_CAN, RHO_PLENUM, f_max=f_max
        )
        got = np.array([m.frequency for m in modes])
        assert len(got) > 0
        for f in got:
            assert np.min(np.abs(fem_all - f)) / f < 3e-3, f"spurious mode at {f:.2f} Hz"
        # Every FEM frequency inside the scanned band (away from its edges) is found.
        # This includes the mode with a pressure node at every junction,
        # sin(k_p s) = 0 at 764 Hz, which the bare admittance sum cannot see.
        fem = fem_all[(fem_all > 20.0) & (fem_all < f_max - 20.0)]
        for f in fem:
            assert np.min(np.abs(got - f)) / f < 3e-3, f"missed FEM mode at {f:.2f} Hz"


def test_push_pull_label_only_for_even_can_count() -> None:
    # For odd N there is no m = N/2 Bloch order; m = (N-1)/2 is a spinning order.
    mode = cb.BlochMode()
    mode.n_cans = 5
    mode.m_azimuthal = 2
    mode.frequency = 100.0
    assert "Push-Pull" not in mode.symmetry_type()
    mode.n_cans = 6
    mode.m_azimuthal = 3
    assert "Push-Pull" in mode.symmetry_type()


# ---------------------------------------------------------------------------
# Orifice / side-branch resonator
# ---------------------------------------------------------------------------


def test_orifice_inertia_uses_both_face_end_corrections() -> None:
    # Thin perforate, no flow: X = k (t + 2 * (8/(3 pi)) r) / sigma = k (t + 0.85 d)/sigma
    # (Rayleigh flanged correction on BOTH faces), the same 0.85 d that
    # helmholtz_frequency applies to a neck.
    f, d, t, sigma, c = 500.0, 1e-3, 1e-3, 0.05, 343.0
    z = cb.orifice_impedance_with_flow(f, 0.0, 0.0, d, t, sigma, 0.7, 1.2, c)
    k = 2.0 * math.pi * f / c
    assert z.imag == pytest.approx(k * (t + 0.85 * d) / sigma, rel=1e-12)


def test_orifice_impedance_refuses_nonpositive_sound_speed() -> None:
    for c in (0.0, -343.0, float("nan")):
        with pytest.raises(RuntimeError):
            cb.orifice_impedance_with_flow(500.0, 0.0, 0.0, 1e-3, 1e-3, 0.05, 0.7, 1.2, c)


@pytest.mark.parametrize("a_tube", [1e-4, 1e-3, 1e-2])
def test_quarter_wave_branch_admittance_matches_lumped_series(a_tube: float) -> None:
    # Branch = neck (referred to the panel area A_duct, as orifice_impedance_with_flow
    # normalises it) in series with a rigid-terminated tube of area A_tube:
    #   Z_branch = rho c [ z_neck / A_duct - i cot(k L) / A_tube ],  T21 = 1 / Z_branch
    freq, L, a_duct, sigma, rho, c = 300.0, 0.2, 0.01, 0.1, 1.2, 343.0
    tm = cb.quarter_wave_resonator_tmm(
        freq=freq,
        L_tube=L,
        A_duct=a_duct,
        A_tube=a_tube,
        d_orifice=0.02,
        l_orifice=0.01,
        porosity=sigma,
        Cd=0.7,
        u_bias=2.0,
        u_grazing=5.0,
        rho=rho,
        c=c,
    )
    k = 2.0 * math.pi * freq / c
    z_neck = cb.orifice_impedance_with_flow(freq, 2.0, 5.0, 0.02, 0.01, sigma, 0.7, rho, c)
    z_branch = rho * c * (z_neck / a_duct + (-1j / math.tan(k * L)) / a_tube)
    assert pytest.approx(1.0 / z_branch, rel=1e-12) == tm.T21


# ---------------------------------------------------------------------------
# Viscothermal Q
# ---------------------------------------------------------------------------

NU, ALPHA, GAMMA = 1.5e-5, 2.1e-5, 1.4


def _delta_eff(f: float) -> float:
    w = 2.0 * math.pi * f
    return math.sqrt(2.0 * NU / w) + (GAMMA - 1.0) * math.sqrt(2.0 * ALPHA / w)


def test_tube_q_matches_kirchhoff_wide_tube() -> None:
    # Kirchhoff wide-tube attenuation alpha = k delta_eff / D; a line resonator
    # has Q = k / (2 alpha) = D / (2 delta_eff), as acoustics.h documents.
    for D, f in [(0.02, 343.0), (0.05, 120.0), (0.01, 900.0)]:
        assert cb.tube_Q(0.25, D, NU, ALPHA, GAMMA, f) == pytest.approx(
            D / (2.0 * _delta_eff(f)), rel=1e-12
        )


def test_helmholtz_q_is_lumped_omega_m_over_r() -> None:
    # Lumped HR: inertance M = rho L_eff / A, neck wall resistance (Kirchhoff,
    # per unit length rho w delta / (a A)) over the physical neck length:
    #   Q = w M / R = a L_eff / (delta_eff L_neck),  a = d/2, L_eff = L + 0.85 d.
    # At a fixed frequency it cannot depend on the cavity volume.
    d, L = 0.02, 0.03
    A = math.pi * d * d / 4.0
    f = 150.0
    L_eff = L + 0.85 * d
    q_ref = (d / 2.0) * L_eff / (_delta_eff(f) * L)
    assert cb.helmholtz_Q(5e-4, A, L, NU, ALPHA, GAMMA, f) == pytest.approx(q_ref, rel=1e-12)
    assert cb.helmholtz_Q(2e-3, A, L, NU, ALPHA, GAMMA, f) == pytest.approx(q_ref, rel=1e-12)
