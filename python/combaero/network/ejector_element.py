"""EjectorElement: 3-port supersonic ejector on the momentum-CV junction
topology (``MultiPortChamberElement``).

Ports: ``primary_node`` (motive, high-pressure inlet), ``secondary_node``
(suction/entrained inlet), ``outlet_node``. Critical (double-choked) mode
only, reusing the validated closed-form physics in
``combaero.network._ejector_huang1999`` (Huang 1999's entrainment ratio,
Eqs. 1-8; Kracik & Dvorak 2016's mixing closure for the critical back
pressure, their Eqs. 7-13) -- see that module's docstring and
``validation/ejector/README.md`` for the full derivation, paper references,
and the accuracy comparison against three alternative closures that were
also evaluated.

Working fluid: combaero's real-fluid EOS covers combustion/air-relevant
species only (N2, O2, Ar, CO2, H2O, light hydrocarbons, CO, H2, NH3) -- no
halocarbon refrigerants. gamma and R are evaluated LIVE from that EOS at
each residual call, at the entrained-flow choking plane (same convention
validated in the reference model, there baked from CoolProp for R141b; here
computed on the fly via ``cb.cp_mass``/``cb.cv_mass``), using the secondary
stream's own composition -- Huang's model assumes a single working fluid, so
this is the representative gamma/R for the whole closed form.

Residual design (see the design discussion that produced this; not
re-derived here). Ports are ``[primary, secondary, outlet]``, so
``_port_signs = [-1, -1, +1]`` and the base class's ``unknowns()``/
``n_equations()`` (``{id}.P_jct``, N+1 = 4 rows) are inherited unchanged --
only the row CONTENT is ejector physics instead of impulse-CV physics:

    R0 = mp - mdot_choked(Pt_primary, Tt_primary, throat_area, gamma, eta_p)
    R1 = ms - omega(Pt_primary, Tt_primary, Pt_secondary, Tt_secondary) * mp
    R2 = mdot_out - mp - ms                        (mass conservation)
    R3 = P_jct - recovery_efficiency * P_c*(...)   (P_jct repurposed as P_c*)

``P_jct`` is NOT tied to the outlet port's actual pressure: in critical mode
the entrainment ratio is by definition independent of back pressure, so the
outlet's actual Pt is legitimately external (set by the downstream network,
same as any ordinary node) -- asserting ``Pt_outlet = P_c*`` would
over-constrain it. ``P_jct`` is instead a diagnostic upper bound, checked
post-solve by ``verify_solution_consistent`` (Pt_outlet must not exceed it --
violation means the design point is subcritical, which this element does not
model).

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

from typing import Any

import combaero as cb
from combaero import _solver_tools
from combaero.network._ejector_huang1999 import ETA_P, EjectorGeometry
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


class EjectorElement(MultiPortChamberElement):
    """Critical-mode supersonic ejector: 2 inlets (primary, secondary), 1 outlet."""

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
        super().__init__(
            id,
            inlet_nodes=[primary_node, secondary_node],
            outlet_nodes=[outlet_node],
        )
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

    def residuals(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        primary, secondary, _outlet = states

        # gamma/r_gas: frozen coefficients (see module docstring). gamma is a
        # weak function of secondary.Tt via the Python fixed point; r_gas is a
        # pure function of composition (not a Newton unknown). Recomputed each
        # call, so the converged root is exact; the C++ Jacobians below hold
        # both fixed, dropping only the tiny implicit d(gamma)/d(Tt) term.
        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        arn = self.area_ratio_nozzle
        arm = self.area_ratio_mix

        mp_choked, dmp_choked_dp_g, dmp_choked_dt_g = (
            _solver_tools.ejector_choked_mass_flow_and_jacobian(
                primary.Pt, primary.Tt, self.throat_area, gamma, r_gas, ETA_P
            )
        )
        entr = _solver_tools.ejector_entrainment_ratio_and_jacobian(
            primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, arn, arm, gamma
        )
        omega = entr.value.omega
        pc = _solver_tools.ejector_critical_back_pressure_and_jacobian(
            primary.Pt,
            primary.Tt,
            secondary.Pt,
            secondary.Tt,
            arn,
            arm,
            gamma,
            r_gas,
            self.recovery_efficiency,
        )
        pc_star = pc.value.p_c_pa

        # port_mdots are junction-sign convention (positive = out of
        # junction); primary/secondary are inlets (sign -1), so their
        # physical mass flow is the negation.
        mp = -float(port_mdots[0])
        ms = -float(port_mdots[1])
        mdot_out = float(port_mdots[2])

        r0 = mp - mp_choked
        r1 = ms - omega * mp
        r2 = mdot_out - mp - ms  # == sum(port_mdots); mass conservation
        r3 = P_jct - pc_star

        # Full analytic Jacobian. Column-name conventions (see solver.py's
        # element jac-relay loop): port mass flows are direct unknowns keyed
        # "{outer_id}.m_dot"; a physical port flow d(mp) = -_port_signs[i]
        # because port_mdots[i] = _port_signs[i] * (that unknown). Port
        # stagnation pressures "{node}.Pt" are direct unknowns; port
        # stagnation temperatures are emitted as "{node}.T" and chained
        # through the solver's temperature relay (same convention the
        # compressible orifice uses: pass Tt, key the column .T). "{id}.P_jct"
        # is this element's own direct unknown.
        m_p = f"{self._port_element_ids[0]}.m_dot"
        m_s = f"{self._port_element_ids[1]}.m_dot"
        pn_pt = f"{self.primary_node}.Pt"
        pn_t = f"{self.primary_node}.T"
        sn_pt = f"{self.secondary_node}.Pt"
        sn_t = f"{self.secondary_node}.T"
        # d(mp)/d(m_p unknown) and d(ms)/d(m_s unknown).
        d_mp = -self._port_signs[0]
        d_ms = -self._port_signs[1]

        row0 = {
            m_p: d_mp,
            pn_pt: -dmp_choked_dp_g,
            pn_t: -dmp_choked_dt_g,
        }
        row1 = {
            m_s: d_ms,
            m_p: -omega * d_mp,
            pn_pt: -mp * entr.domega_dp_g,
            pn_t: -mp * entr.domega_dt_g,
            sn_pt: -mp * entr.domega_dp_e,
            sn_t: -mp * entr.domega_dt_e,
        }
        mass_row: dict[str, float] = {}
        for i in range(self.N):
            mdot_var = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mdot_var] = mass_row.get(mdot_var, 0.0) + self._port_signs[i]
        row3 = {
            f"{self.id}.P_jct": 1.0,
            pn_pt: -pc.dpc_dp_g,
            pn_t: -pc.dpc_dt_g,
            sn_pt: -pc.dpc_dp_e,
            sn_t: -pc.dpc_dt_e,
        }

        return [r0, r1, r2, r3], {0: row0, 1: row1, 2: mass_row, 3: row3}

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        rel_tol: float = 1e-3,
    ) -> bool:
        """Post-solve check: the outlet's actual back pressure must not
        exceed the computed critical back pressure P_c* (stored in
        ``{id}.P_jct``). A violation means the converged operating point is
        subcritical, which this element does not model (critical-mode
        omega, used in the residuals, is only valid up to P_c*). Missing
        keys => True (do not police incomplete dicts, matching the
        convention used elsewhere in this module)."""
        pc_key = f"{self.id}.P_jct"
        pt_key = f"{self.outlet_node}.Pt"
        if pc_key not in sol or pt_key not in sol:
            return True
        return float(sol[pt_key]) <= float(sol[pc_key]) * (1.0 + rel_tol)

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
        entr = _solver_tools.ejector_entrainment_ratio_and_jacobian(
            primary.Pt,
            primary.Tt,
            secondary.Pt,
            secondary.Tt,
            self.area_ratio_nozzle,
            self.area_ratio_mix,
            gamma,
        )

        diag["gamma"] = gamma
        diag["r_gas"] = r_gas
        diag["omega"] = entr.value.omega
        diag["p_c_star_pa"] = float(P_jct)
        diag["outlet_pt_pa"] = float(outlet.Pt)
        diag["critical_mode"] = 1.0 if float(outlet.Pt) <= float(P_jct) else 0.0
        if port_mdots is not None and len(port_mdots) == self.N:
            diag["m_dot_primary"] = -float(port_mdots[0])
            diag["m_dot_secondary"] = -float(port_mdots[1])
            diag["m_dot_outlet"] = float(port_mdots[2])
        return diag
