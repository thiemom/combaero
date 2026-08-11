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

Jacobian: FD fallback (``jac = {}`` for rows 0, 1, 3; row 2's mass
conservation is linear so it gets an analytic entry). This mirrors
``MPCEv2Element``'s own documented prototype pattern -- analytic
(sympy-derived) Jacobians are planned for the C++ port, the next phase of
this work, not this element.
"""

from __future__ import annotations

from typing import Any

import combaero as cb
from combaero.network._ejector_huang1999 import (
    ETA_P,
    EjectorGeometry,
    choked_mass_flow,
    critical_back_pressure,
    entrainment_ratio,
)
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

        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        geom = self.geometry

        mp_choked = choked_mass_flow(
            primary.Pt, primary.Tt, self.throat_area, gamma, r_gas, eta=ETA_P
        )
        omega = entrainment_ratio(
            primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma
        ).omega
        pc_star = critical_back_pressure(
            primary.Pt,
            primary.Tt,
            secondary.Pt,
            secondary.Tt,
            geom,
            gamma,
            r_gas,
            self.recovery_efficiency,
        ).p_c_pa

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

        # Analytic Jacobian only for the linear mass row; rows 0/1/3 are FD
        # fallback -- documented prototype pattern (see module docstring),
        # analytic Jacobians land in the C++ port.
        mass_row: dict[str, float] = {}
        for i in range(self.N):
            mdot_var = f"{self._port_element_ids[i]}.m_dot"
            mass_row[mdot_var] = mass_row.get(mdot_var, 0.0) + self._port_signs[i]

        return [r0, r1, r2, r3], {2: mass_row}

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
        geom = self.geometry
        entr = entrainment_ratio(primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma)

        diag["gamma"] = gamma
        diag["r_gas"] = r_gas
        diag["omega"] = entr.omega
        diag["p_c_star_pa"] = float(P_jct)
        diag["outlet_pt_pa"] = float(outlet.Pt)
        diag["critical_mode"] = 1.0 if float(outlet.Pt) <= float(P_jct) else 0.0
        if port_mdots is not None and len(port_mdots) == self.N:
            diag["m_dot_primary"] = -float(port_mdots[0])
            diag["m_dot_secondary"] = -float(port_mdots[1])
            diag["m_dot_outlet"] = float(port_mdots[2])
        return diag
