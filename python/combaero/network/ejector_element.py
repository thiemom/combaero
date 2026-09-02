"""EjectorElement: 3-port supersonic ejector on the momentum-CV junction
topology (``MultiPortChamberElement``).

Ports: ``primary_node`` (motive, high-pressure inlet), ``secondary_node``
(suction/entrained inlet), ``outlet_node``. Reuses the validated closed-form
physics in ``combaero.network._ejector_huang1999`` (Huang 1999's entrainment
ratio, Eqs. 1-8; Kracik & Dvorak 2016's mixing closure for the critical back
pressure, their Eqs. 7-13) -- see that module's docstring and
``validation/ejector/README.md`` for the full derivation, paper references,
and the accuracy comparison against three alternative closures that were
also evaluated.

Operating regimes. This element resolves all three physical regimes in one
Newton-solvable residual system (design + provenance in
``validation/ejector/OPERATING_REGIMES_DESIGN.md``):

  * CRITICAL (double-choked): primary choked, entrained flow critical, so
    ``omega`` is back-pressure-independent. The outlet floats freely below
    the ejector's critical back pressure ``P_c*``.
  * SUBCRITICAL droop: primary still choked but the outlet is pushed above
    ``P_c*``; the outlet pressure sets the operating point.
  * UNCHOKED jet pump: the forced primary mass flow is below the choke
    threshold ``mdot_choked(P_back)``, so the primary Pt floats to just above
    the back pressure (forward flow, ``dp >= 0``) and the device acts as a
    subsonic constant-area jet pump (Keenan-Neumann-Lustwerk mixing via
    ``ejector_jetpump_discharge``).

The regime is selected smoothly (no branching discontinuity) by two
smootherstep weights: ``s_choke = smootherstep(mp / mdot_choked(Pt_p))``
(0 = unchoked jet pump, 1 = choked) and ``s_sub`` (0 = outlet below ``P_c*``,
critical; 1 = outlet above it, subcritical). See the residual design below.

Working fluid: combaero's real-fluid EOS covers combustion/air-relevant
species only (N2, O2, Ar, CO2, H2O, light hydrocarbons, CO, H2, NH3) -- no
halocarbon refrigerants. gamma and R are evaluated LIVE from that EOS at
each residual call, at the entrained-flow choking plane (same convention
validated in the reference model, there baked from CoolProp for R141b; here
computed on the fly via ``cb.cp_mass``/``cb.cv_mass``), using the secondary
stream's own composition -- Huang's model assumes a single working fluid, so
this is the representative gamma/R for the whole closed form.

Residual design. Ports are ``[primary, secondary, outlet]``, so
``_port_signs = [-1, -1, +1]`` and the base class's ``unknowns()``/
``n_equations()`` (one owned unknown ``{id}.P_jct``, N+1 = 4 rows) are
inherited unchanged -- but the owned unknown is REPURPOSED as the physical
mixing-plane static pressure ``P_py`` (the pressure common to the primary jet
and the entrained stream at the mixing-chamber entrance), and the row content
is ejector physics instead of impulse-CV physics. With ``s = s_choke``:

    R0 = mp - cd_nozzle(Pt_p, P_py, A_t, A_e)      (primary C-D nozzle)
    R1 = ms - [s*omega_crit*mp + (1-s)*ms_sec]     (blended entrainment)
    R2 = mdot_out - mp - ms                        (mass conservation)
    R3 = [w_pin*Pt_outlet + (1-w_pin)*P_py] - recovery*discharge

``P_py`` is a genuine Newton unknown (not derived): it is the quantity the two
nozzle relations R0/R1 share, so making it free keeps BOTH live and lets the
downstream back pressure select the operating regime. ``R0``'s C-D nozzle
smooth-mins the subsonic-flux and sonic-cap branches, so it is choked or
unchoked automatically as ``P_py`` moves. ``discharge = recovery * (s*P_c* +
(1-s)*P03_jetpump)`` blends the critical and jet-pump mixing closures; each is
evaluated ONLY where its weight is nonzero (the critical closures return NaN
for an unchoked primary, and ``0 * NaN = NaN`` would otherwise poison the
Jacobian -- safe because smootherstep is flat at ``s in {0, 1}``).

``R3`` closes the system differently per regime via
``w_pin = 1 - s*(1 - s_sub)`` (the weight on the outlet-pin form):

  * jet pump (s=0) and subcritical droop (s=1, s_sub=1): ``w_pin = 1`` ->
    ``Pt_outlet = recovery*discharge``. The outlet is what the ejector sets;
    the downstream boundary then forces ``P_py`` through this pin, and R0/R1
    set the primary Pt and entrainment. This pin is what excludes the spurious
    double-choked root (primary Pt below the back pressure) that a
    back-pressure-independent closure would otherwise admit.
  * critical (s=1, outlet below ``P_c*``, s_sub=0): ``w_pin = 0`` ->
    ``P_py = recovery*P_c*`` (a diagnostic) and the outlet floats FREE, set by
    the downstream network. Pinning it here would over-constrain the fixed
    outlet boundary (``Pt_outlet < P_c*`` legitimately).

``verify_solution_consistent`` now returns True unconditionally: every regime
is modeled, so there is no longer a subcritical design point to demote.

Jacobian: FULL analytic, all 4 rows, via the C++ (f, J) path
(``combaero.network._ejector_huang1999`` is now only the validation
reference; the hot path calls ``_solver_tools.ejector_*_and_jacobian``,
i.e. ``include/ejector.h``/``src/ejector.cpp``). The C++ Jacobians are
exact w.r.t. the four thermodynamic inputs (primary/secondary stagnation
pressure and temperature) holding ``gamma`` and ``r_gas`` FROZEN. Here
``gamma = choke_plane_gamma(secondary.Tt, ...)`` is a weak function of a
Newton unknown and ``r_gas`` depends only on composition (not a solved
unknown), so the element Jacobian is exact in the explicit (Pt, Tt)
dependence and drops only the tiny implicit ``d gamma / d Tt`` term -- a
standard frozen-property Jacobian. ``gamma`` is recomputed fresh every
residual call, so the converged root is still exact; only the Newton step
uses a slightly approximate slope. This matches the reference data's own
convention (``generate_reference.py`` treats gamma as a separate Jacobian
input, so its d/dt_e central-difference targets also hold gamma fixed).
``choke_plane_gamma`` therefore stays in Python.
"""

from __future__ import annotations

import math
from typing import Any

import combaero as cb
from combaero import _solver_tools
from combaero.network._ejector_huang1999 import ETA_P, ETA_S, EjectorGeometry
from combaero.network.components import MultiPortChamberElement, NetworkMixtureState


def _real_gamma(t: float, x: Any) -> float:
    """Ideal-gas-like gamma = cp/cv from combaero's own EOS at temperature
    t and composition x (mole fractions). cb.cp_mass/cv_mass are functions
    of (T, X) only, matching the pressure-independent NASA-polynomial-style
    cp0(T) convention the reference model's CoolProp-baked gamma also uses.
    """
    cp = cb.cp_mass(t, x)
    cv = cb.cv_mass(t, x)
    return cp / cv if cv > 0.0 else 1.4


def choke_plane_gamma(t_e: float, x: Any, *, iterations: int = 12) -> float:
    """gamma at the entrained-flow choking plane y-y (T_sy = T_e/((g+1)/2)),
    a self-consistent fixed point since gamma appears on both sides through
    T_sy -- same convention validated in the reference model (see
    ``validation/ejector/README.md``'s "gamma provenance" section), here
    evaluated live from combaero's EOS instead of baked from CoolProp.
    """
    gamma = _real_gamma(t_e, x)
    for _ in range(iterations):
        t_sy = t_e / ((gamma + 1.0) / 2.0)
        gamma = _real_gamma(t_sy, x)
    return gamma


class _D:
    """Minimal forward-mode dual: a value plus a gradient dict over the element's
    primitive unknowns (mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out). Lets the
    operating-regime residuals be assembled by chaining the C++ closures' analytic
    partials, exactly mirroring the FD-validated prototype (see
    validation/ejector/OPERATING_REGIMES_DESIGN.md sec 6b). Only the operators
    the assembly uses are defined."""

    __slots__ = ("v", "d")

    def __init__(self, v: float, d: dict[str, float] | None = None) -> None:
        self.v = float(v)
        self.d = dict(d) if d else {}

    @staticmethod
    def seed(v: float, name: str) -> _D:
        return _D(v, {name: 1.0})

    def __add__(self, b: _D | float) -> _D:
        if isinstance(b, _D):
            d = dict(self.d)
            for k, x in b.d.items():
                d[k] = d.get(k, 0.0) + x
            return _D(self.v + b.v, d)
        return _D(self.v + b, self.d)

    __radd__ = __add__

    def __neg__(self) -> _D:
        return _D(-self.v, {k: -x for k, x in self.d.items()})

    def __sub__(self, b: _D | float) -> _D:
        return self + (-b if isinstance(b, _D) else -b)

    def __rsub__(self, b: float) -> _D:
        return (-self) + b

    def __mul__(self, b: _D | float) -> _D:
        if isinstance(b, _D):
            d = {k: x * b.v for k, x in self.d.items()}
            for k, x in b.d.items():
                d[k] = d.get(k, 0.0) + self.v * x
            return _D(self.v * b.v, d)
        return _D(self.v * b, {k: x * b for k, x in self.d.items()})

    __rmul__ = __mul__

    def __truediv__(self, b: _D | float) -> _D:
        if isinstance(b, _D):
            inv = 1.0 / b.v
            d = {k: x * inv for k, x in self.d.items()}
            for k, x in b.d.items():
                d[k] = d.get(k, 0.0) - self.v * x * inv * inv
            return _D(self.v * inv, d)
        return _D(self.v / b, {k: x / b for k, x in self.d.items()})


def _chain(value: float, pairs: list[tuple[float, _D]]) -> _D:
    """Wrap a C++ closure output: value + sum(partial * grad(input dual))."""
    d: dict[str, float] = {}
    for p, inp in pairs:
        for k, x in inp.d.items():
            d[k] = d.get(k, 0.0) + p * x
    return _D(value, d)


def _smootherstep(t: _D, lo: float, hi: float) -> _D:
    """C1 clamped smootherstep on t, saturating to EXACTLY 0/1 outside [lo, hi]
    (so s_choke = 1 exactly at the primary-choke threshold and critical mode
    reduces exactly)."""
    x = (t - lo) / (hi - lo)
    if x.v <= 0.0:
        return _D(0.0)
    if x.v >= 1.0:
        return _D(1.0)
    val = x.v * x.v * x.v * (x.v * (x.v * 6.0 - 15.0) + 10.0)
    dval = 30.0 * x.v * x.v * (x.v - 1.0) * (x.v - 1.0)
    return _chain(val, [(dval, x)])


class EjectorElement(MultiPortChamberElement):
    """Supersonic ejector across operating regimes: 2 inlets (primary, secondary),
    1 outlet. Critical (double-choked), subcritical, and unchoked-primary
    (subsonic jet-pump) in one C1 residual system."""

    def __init__(
        self,
        id: str,
        primary_node: str,
        secondary_node: str,
        outlet_node: str,
        throat_area: float,
        nozzle_exit_area: float,
        mixing_area: float,
        recovery_efficiency: float = 1.0,
    ):
        """Build an ejector.

        Args:
            id: Element id.
            primary_node: Motive/high-pressure inlet port-MCN id.
            secondary_node: Suction/entrained inlet port-MCN id.
            outlet_node: Outlet port-MCN id.
            throat_area: Primary nozzle throat area A_t [m^2].
            nozzle_exit_area: Primary nozzle exit area A_p1 [m^2]. Must
                exceed throat_area (the nozzle must diverge to expand
                supersonically).
            mixing_area: Constant-area mixing-section area A_3 [m^2]. Must
                exceed nozzle_exit_area (room for the entrained flow).
            recovery_efficiency: Multiplies the lossless mixed stagnation
                pressure to give the critical back pressure P_c*. Defaults
                to 1.0 (no artificial loss); see this module's and
                ``_ejector_huang1999``'s docstrings for why that default is
                not fitted to any reference dataset.
        """
        # Validate geometry BEFORE super().__init__ so the port-face areas
        # passed below are guaranteed positive.
        if throat_area <= 0.0:
            raise ValueError(f"EjectorElement '{id}': throat_area must be positive.")
        if nozzle_exit_area <= throat_area:
            raise ValueError(
                f"EjectorElement '{id}': nozzle_exit_area ({nozzle_exit_area}) must "
                f"exceed throat_area ({throat_area}) -- the primary nozzle must "
                f"diverge to expand supersonically."
            )
        if mixing_area <= nozzle_exit_area:
            raise ValueError(
                f"EjectorElement '{id}': mixing_area ({mixing_area}) must exceed "
                f"nozzle_exit_area ({nozzle_exit_area}) -- the constant-area section "
                f"must leave room for the entrained flow."
            )
        if recovery_efficiency <= 0.0:
            raise ValueError(
                f"EjectorElement '{id}': recovery_efficiency must be positive, "
                f"got {recovery_efficiency}."
            )
        # Port-face areas from the ejector's OWN geometry (ordered
        # inlet_nodes + outlet_nodes = [primary, secondary, outlet]). Unlike a
        # generic tee, the ejector always knows these, so it supplies them
        # rather than relying on the base class inferring them from a
        # connecting ChannelElement -- an ejector can be fed by a bare
        # boundary + connection with no channel/diameter to infer from.
        # These set the port-MCN dynamic-head (Pt = P + 0.5*rho*v^2) only;
        # the ejector's own residuals use throat/nozzle_exit/mixing directly.
        #   primary  -> nozzle exit A_p1 (motive jet plane)
        #   secondary-> annular suction A_3 - A_p1
        #   outlet   -> mixing/diffuser duct A_3
        port_areas = [
            float(nozzle_exit_area),
            float(mixing_area - nozzle_exit_area),
            float(mixing_area),
        ]
        super().__init__(
            id,
            inlet_nodes=[primary_node, secondary_node],
            outlet_nodes=[outlet_node],
            port_areas=port_areas,
        )
        self.primary_node = primary_node
        self.secondary_node = secondary_node
        self.outlet_node = outlet_node
        self.throat_area = float(throat_area)
        self.nozzle_exit_area = float(nozzle_exit_area)
        self.mixing_area = float(mixing_area)
        self.recovery_efficiency = float(recovery_efficiency)

    @property
    def area_ratio_nozzle(self) -> float:
        """A_p1 / A_t."""
        return self.nozzle_exit_area / self.throat_area

    @property
    def area_ratio_mix(self) -> float:
        """A_3 / A_t."""
        return self.mixing_area / self.throat_area

    @property
    def geometry(self) -> EjectorGeometry:
        return EjectorGeometry(self.area_ratio_nozzle, self.area_ratio_mix)

    def row_scale_kinds(self) -> list[str]:
        """Rows 0-2 (choked flow, entrainment closure, mass conservation)
        are mass-flow-valued; row 3 (P_jct - P_c*) is pressure-valued --
        the OPPOSITE pattern from the base class's own impulse-CV rows.
        See the module docstring's residual design section and the
        solver.py scaling block this overrides."""
        return ["mdot", "mdot", "mdot", "p"]

    # smootherstep window for s_choke, on mp / choked_mass_flow(Pt_p): s = 0
    # when clearly unchoked (jet pump), 1 exactly at the primary-choke threshold.
    _S_CHOKE_LO = 0.90
    _S_CHOKE_HI = 0.999
    # s_sub window on outlet.Pt / (recovery * P_c*): 0 (critical, outlet free)
    # below the recovered critical back pressure, 1 (subcritical, outlet pinned)
    # above it. Centred on 1 with a narrow band for a sharp-but-smooth handover.
    _S_SUB_LO = 0.98
    _S_SUB_HI = 1.02
    _EPS_FRAC = 1.0e-3  # C-D nozzle / secondary-flux smooth-min rounding

    def residuals(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        primary, secondary, outlet = states

        # gamma/r_gas: frozen coefficients (see module docstring). Recomputed
        # each call so the converged root is exact; the closure Jacobians hold
        # both fixed (drops only the tiny implicit d(gamma)/d(Tt) term).
        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        arn = self.area_ratio_nozzle
        arm = self.area_ratio_mix
        A_t = self.throat_area
        A_e = self.nozzle_exit_area
        A_s = self.mixing_area - self.nozzle_exit_area
        rec = self.recovery_efficiency

        # Primitive unknowns as duals. port_mdots are junction-sign convention
        # (positive = out); primary/secondary are inlets (sign -1), so their
        # physical flow is the negation.
        mp_v = -float(port_mdots[0])
        ms_v = -float(port_mdots[1])
        mo_v = float(port_mdots[2])
        mp = _D.seed(mp_v, "mp")
        ms = _D.seed(ms_v, "ms")
        mdot_out = _D.seed(mo_v, "mdot_out")
        p_g = _D.seed(primary.Pt, "p_g")
        t_g = _D.seed(primary.Tt, "t_g")
        p_e = _D.seed(secondary.Pt, "p_e")
        t_e = _D.seed(secondary.Tt, "t_e")
        p_out = _D.seed(outlet.Pt, "p_out")

        # The owned unknown (base name "{id}.P_jct") is repurposed as the
        # PHYSICAL mixing-plane static pressure P_py -- the pressure common to
        # the primary jet and the entrained stream at the mixing-chamber
        # entrance. This is the quantity the two nozzle relations R0/R1 share,
        # so making it a genuine Newton unknown (instead of deriving it from
        # the primary nozzle inverse, which makes R0 vacuous in the unchoked
        # branch) keeps R0 AND R1 live and lets the downstream back pressure
        # select the operating regime. See OPERATING_REGIMES_DESIGN.md sec 6c.
        p_py = _D.seed(P_jct, "P_jct")

        # s_choke = smootherstep(mp / cap), cap = choked_mass_flow(Pt_p): 0 when
        # the primary is clearly unchoked (jet pump), 1 at the choke threshold.
        cap_v, dcap_dp_g, dcap_dt_g = _solver_tools.ejector_choked_mass_flow_and_jacobian(
            primary.Pt, primary.Tt, A_t, gamma, r_gas, ETA_P
        )
        cap = _chain(cap_v, [(dcap_dp_g, p_g), (dcap_dt_g, t_g)])
        s = _smootherstep(mp / cap, self._S_CHOKE_LO, self._S_CHOKE_HI)

        # R0: primary C-D nozzle.
        cd = _solver_tools.ejector_cd_nozzle_mass_flow_and_jacobian(
            primary.Pt, primary.Tt, p_py.v, A_t, A_e, gamma, r_gas, ETA_P, self._EPS_FRAC
        )
        cd_d = _chain(cd.mdot, [(cd.dmdot_dp0, p_g), (cd.dmdot_dt0, t_g), (cd.dmdot_dp_py, p_py)])
        r0 = mp - cd_d

        # R1/R3 blend the critical (choked-primary) and jet-pump (unchoked)
        # closures by s = s_choke. The critical closures
        # (``ejector_entrainment_ratio``, ``ejector_critical_back_pressure``)
        # assume a choked primary and return NaN once it is not (e.g. Pt_primary
        # ~ back pressure in the jet-pump root); the jet-pump closure is
        # symmetric at the choked end. Since smootherstep is FLAT at s in
        # {0, 1} (both the blend weight AND its gradient vanish there), the
        # inactive branch contributes exactly zero -- but ``0 * NaN = NaN``
        # would still poison the value and Jacobian. So evaluate each branch
        # ONLY where its weight is nonzero; the skipped branch is a hard zero
        # dual, which is C1-consistent with the flat smootherstep endpoint.
        omega_eff = ms / mp  # ACTUAL ratio for the discharge (R1 pins ms itself)
        crit_active = s.v > 0.0
        jet_active = s.v < 1.0

        if crit_active:
            entr = _solver_tools.ejector_entrainment_ratio_and_jacobian(
                primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, arn, arm, gamma
            )
            omega_crit = _chain(
                entr.value.omega,
                [
                    (entr.domega_dp_g, p_g),
                    (entr.domega_dt_g, t_g),
                    (entr.domega_dp_e, p_e),
                    (entr.domega_dt_e, t_e),
                ],
            )
            pc = _solver_tools.ejector_critical_back_pressure_and_jacobian(
                primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, arn, arm, gamma, r_gas, 1.0
            )
            p03_crit = _chain(
                pc.value.p_mixed_stagnation_pa,
                [(pc.dpc_dp_g, p_g), (pc.dpc_dt_g, t_g), (pc.dpc_dp_e, p_e), (pc.dpc_dt_e, t_e)],
            )
            # s_sub: 0 while the outlet stays below the recovered critical back
            # pressure (genuine critical mode -- outlet floats free), ramping to
            # 1 as it rises above it (subcritical droop -- outlet pinned).
            s_sub = _smootherstep(p_out / (rec * p03_crit), self._S_SUB_LO, self._S_SUB_HI)
        else:
            omega_crit = _D(0.0)
            p03_crit = _D(0.0)
            s_sub = _D(1.0)

        if jet_active:
            sec = _solver_tools.ejector_cd_nozzle_mass_flow_and_jacobian(
                secondary.Pt, secondary.Tt, p_py.v, A_s, A_s, gamma, r_gas, ETA_S, self._EPS_FRAC
            )
            ms_sec = _chain(
                sec.mdot, [(sec.dmdot_dp0, p_e), (sec.dmdot_dt0, t_e), (sec.dmdot_dp_py, p_py)]
            )
            jp = _solver_tools.ejector_jetpump_discharge_and_jacobian(
                primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, p_py.v, omega_eff.v, gamma, 1.0
            )
            p03_jet = _chain(
                jp.p03,
                [
                    (jp.dp03_dp_g, p_g),
                    (jp.dp03_dt_g, t_g),
                    (jp.dp03_dp_e, p_e),
                    (jp.dp03_dt_e, t_e),
                    (jp.dp03_dp_py, p_py),
                    (jp.dp03_domega, omega_eff),
                ],
            )
        else:
            ms_sec = _D(0.0)
            p03_jet = _D(0.0)

        # R1: blended entrainment  ms = s*omega_crit*mp + (1-s)*ms_sec.
        ms_model = s * omega_crit * mp + (1.0 - s) * ms_sec
        r1 = ms - ms_model

        # R2: mass conservation.
        r2 = mdot_out - mp - ms

        # R3: regime-dependent closure of the owned unknown.
        #   discharge = recovery * blend(P_c*, P03_jet) by s_choke.
        #   w_pin = 1 - s*(1 - s_sub): the weight on the OUTLET-PIN form.
        #     - jet-pump (s=0):            w_pin = 1  -> pin the outlet
        #     - subcritical droop (s=1,    w_pin = 1  -> pin the outlet
        #        outlet > P_c*, s_sub=1)
        #     - critical (s=1, outlet <    w_pin = 0  -> P_jct diagnostic, outlet
        #        P_c*, s_sub=0)                          free (set by the network)
        # In the critical branch the outlet is legitimately below the ejector's
        # critical back pressure, so pinning it would over-constrain the fixed
        # downstream boundary; there the owned unknown is instead the diagnostic
        # recovered stagnation. Everywhere the primary is unchoked, or choked
        # but back-pressure-loaded above P_c*, the outlet IS what the ejector
        # sets, so it is pinned and P_py floats to satisfy the two nozzles.
        discharge = rec * (s * p03_crit + (1.0 - s) * p03_jet)
        w_pin = 1.0 - s * (1.0 - s_sub)
        lhs = w_pin * p_out + (1.0 - w_pin) * p_py
        r3 = lhs - discharge

        # Map each residual's gradient (over the primitives) to solver columns.
        m_p = f"{self._port_element_ids[0]}.m_dot"
        m_s = f"{self._port_element_ids[1]}.m_dot"
        m_out = f"{self._port_element_ids[2]}.m_dot"
        pn_pt, pn_t = f"{self.primary_node}.Pt", f"{self.primary_node}.T"
        sn_pt, sn_t = f"{self.secondary_node}.Pt", f"{self.secondary_node}.T"
        on_pt = f"{self.outlet_node}.Pt"
        d_mp = -self._port_signs[0]  # d(mp_phys)/d(m_p unknown)
        d_ms = -self._port_signs[1]
        d_mo = self._port_signs[2]  # d(mdot_out_phys)/d(m_out unknown)
        col_of = {
            "mp": (m_p, d_mp),
            "ms": (m_s, d_ms),
            "mdot_out": (m_out, d_mo),
            "p_g": (pn_pt, 1.0),
            "t_g": (pn_t, 1.0),
            "p_e": (sn_pt, 1.0),
            "t_e": (sn_t, 1.0),
            "p_out": (on_pt, 1.0),
            "P_jct": (f"{self.id}.P_jct", 1.0),
        }

        def to_columns(res: _D) -> dict[str, float]:
            cols: dict[str, float] = {}
            for var, coeff in res.d.items():
                name, chain = col_of[var]
                if coeff != 0.0:
                    cols[name] = cols.get(name, 0.0) + coeff * chain
            return cols

        return (
            [r0.v, r1.v, r2.v, r3.v],
            {0: to_columns(r0), 1: to_columns(r1), 2: to_columns(r2), 3: to_columns(r3)},
        )

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        rel_tol: float = 1e-3,
    ) -> bool:
        """No post-solve demotion. The element now models the critical,
        subcritical, and unchoked-primary (jet-pump) regimes in one C1 residual
        system, so a converged root is physical by construction (the old
        ``outlet.Pt <= P_c*`` check demoted valid subcritical/jet-pump points
        and is dropped -- see OPERATING_REGIMES_DESIGN.md sec 6b)."""
        return True

    def diagnostics(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Parent per-port fields plus ejector-specific outputs: entrainment
        ratio, live gamma/R, and the critical-vs-actual back-pressure margin.
        """
        diag = super().diagnostics(states, P_jct, port_mdots)
        primary, secondary, outlet = states

        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        diag["gamma"] = gamma
        diag["r_gas"] = r_gas
        diag["outlet_pt_pa"] = float(outlet.Pt)
        diag["p_py_pa"] = float(P_jct)  # owned unknown = mixing-plane static

        # Operating regime: s_choke = smootherstep(mp/cap), reported so the GUI
        # can label critical (1) / jet-pump (0) / transition, plus the actual
        # entrainment ratio and the choked-mode reference quantities. The
        # critical closures (omega_critical, P_c*) assume a choked primary and
        # return NaN in the jet-pump regime, so they are reported only when
        # meaningful (s_choke > 0).
        if port_mdots is not None and len(port_mdots) == self.N:
            mp = -float(port_mdots[0])
            ms = -float(port_mdots[1])
            cap = _solver_tools.ejector_choked_mass_flow_and_jacobian(
                primary.Pt, primary.Tt, self.throat_area, gamma, r_gas, ETA_P
            )[0]
            t = (mp / cap - self._S_CHOKE_LO) / (self._S_CHOKE_HI - self._S_CHOKE_LO)
            t = min(max(t, 0.0), 1.0)
            s_choke = t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
            diag["s_choke"] = s_choke
            diag["m_dot_primary"] = mp
            diag["m_dot_secondary"] = ms
            diag["m_dot_outlet"] = float(port_mdots[2])
            # ACTUAL entrainment ratio -- regime-independent, always valid.
            diag["omega"] = ms / mp if mp != 0.0 else float("nan")
            # Regime label. critical_mode = 1 only when the primary is choked
            # AND the outlet sits below the ejector's critical back pressure
            # (double-choked plateau); 0 in the subcritical-droop and unchoked
            # jet-pump regimes, where the outlet pressure sets the operating
            # point.
            critical_mode = 0.0
            if s_choke > 0.0:
                omega_c = self._critical_omega(primary, secondary, gamma)
                p_c_star = self._critical_back_pressure(primary, secondary, gamma, r_gas)
                diag["omega_critical"] = omega_c
                diag["p_c_star_pa"] = p_c_star
                if s_choke >= 0.5 and math.isfinite(p_c_star) and float(outlet.Pt) <= p_c_star:
                    critical_mode = 1.0
            diag["critical_mode"] = critical_mode
        else:
            diag["omega"] = self._critical_omega(primary, secondary, gamma)
            diag["p_c_star_pa"] = self._critical_back_pressure(primary, secondary, gamma, r_gas)
        return diag

    def _critical_omega(
        self, primary: NetworkMixtureState, secondary: NetworkMixtureState, gamma: float
    ) -> float:
        return float(
            _solver_tools.ejector_entrainment_ratio_and_jacobian(
                primary.Pt,
                primary.Tt,
                secondary.Pt,
                secondary.Tt,
                self.area_ratio_nozzle,
                self.area_ratio_mix,
                gamma,
            ).value.omega
        )

    def _critical_back_pressure(
        self,
        primary: NetworkMixtureState,
        secondary: NetworkMixtureState,
        gamma: float,
        r_gas: float,
    ) -> float:
        return float(
            _solver_tools.ejector_critical_back_pressure_and_jacobian(
                primary.Pt,
                primary.Tt,
                secondary.Pt,
                secondary.Tt,
                self.area_ratio_nozzle,
                self.area_ratio_mix,
                gamma,
                r_gas,
                self.recovery_efficiency,
            ).value.p_c_pa
        )
