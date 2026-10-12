"""Ideal momentum junction: a parameter-free closure beside the shipped tee (#517).

One junction static pressure ``p_jn`` (the element's ``P_jct`` unknown, which
for this closure is a STATIC pressure, not the stagnation pressure the Mynard
and constant-K tees carry in the same slot). Port ``i`` has an outward axis
``n_i`` and a signed flow ``m_i`` (positive out of the junction):

    supplier  k:  P_k - p_jn = 0
    receiver  j:  P_j - p_jn + (m_j / A_j) (u_j - V . n_j) = 0
    mass:         sum_i m_i = 0

with ``V`` the mass-weighted velocity vector of the supplied streams. The rows
are an exact momentum balance on each receiver in static pressure, with each
port's own density. Merge reduces to Blaisdell & Manson's theory; split gives
every receiver the vector Borda-Carnot loss ``rho/2 |v_j - v_s|^2``
(Idelchik Eq 7-6 with A' = 1). Nothing is fitted, and the model is not meant to
beat a tuned closure: it is a baseline, and a fallback where no tuned loss
applies.

Roles are not declared. Every port carries both, blended smoothly over 1e-4 of
the throughput (``include/ideal_junction.h``), so the residual is smooth through
flow reversal and at zero flow: no flow ratio, no direction classification, no
degenerate fallback, no soft barrier. The inlet/outlet declaration fixes the
canonical graph direction and the geometry only: a converged answer may run any
port either way (the solver feeds a reversed junction port by the sign of its
flow).
"""

from __future__ import annotations

import math
from functools import cache

import combaero as cb
from combaero import _core
from combaero.network.components import (
    MultiPortChamberBase,
    NetworkMixtureState,
)

#: Mirrors ``kIdealMaxPorts`` in include/ideal_junction.h; the kernel refuses
#: more ports, so the element refuses them at construction instead of mid-solve.
IDEAL_MAX_PORTS: int = 6


@cache
def _species_molar_mass(k: int) -> float:
    return float(cb.species_molar_mass(k))


class IdealMomentumTeeElement(MultiPortChamberBase):
    """Parameter-free momentum junction, parallel to ``MultiPortChamberElement``.

    Angles follow the tee convention of the other junction elements: each port's
    geometric angle from the main axis, inlets and outlets alike, so a straight
    run is 0 on both sides and a normal branch is 90. The outward axis the
    kernel needs is that angle for an outlet and ``pi - angle`` for an inlet --
    an inlet's flow arrives travelling along its angle, so it points away from
    the junction the other way. This reproduces the deflection convention of
    Blaisdell & Manson and of the validation tee (``joining_variants.tee``).

    Residual and Jacobian come whole from ``_core.ideal_junction_residuals_and_jacobian``.
    The rows depend on port static pressures, port flows and ``p_jn`` only; the
    port Pt does not enter. Density derivatives are ``rho / P`` (ideal gas, the
    same rule ``MultiPortChamberElement`` uses).

    RANGE: subsonic supplying ports. Every supplier sits at the junction
    static pressure, so a supplier that would really choke -- under-expanded,
    its static above the junction's -- has no state here; the solve drives it
    supersonic and stalls at Mach 1, where the port's mass flux peaks. Measured
    on Wang 2014 a = 2.44, q >= 0.8, M3 >= 0.4 (lateral driven to Mach ~1.55;
    Wang's own rig choked there, p. 1606). Not refused or flagged; a choked-
    supplier branch would be a model extension (#517).

    LOSS FACTOR (#526). ``loss_factor`` (one value for every port, or one per
    port in inlets-then-outlets order) scales each RECEIVER's loss against the
    mixed supplied stream, ``Pt_j = Pt_mix - f_j rho_j/2 |v_j - V|^2``, while
    the suppliers' mixing loss stays at its ideal value. The default 1.0 is the
    parameter-free closure exactly; any f >= 0 keeps the junction's mechanical
    dissipation non-negative, so no guard is needed. A port keeps its factor
    when it reverses, and the factor acts only while the port receives. In a
    split (one supplier) there is no mixing, so f scales the whole loss --
    and f -> 0 meets the energy check: the loss is written at low Mach
    (rho u^2/2) while the port Pt is isentropic, ~M^2/4 of the head larger,
    so a split below f ~ M^2/4 creates flow work and is refused. The
    factor is the user's rig-matching knob, like ``f_multiplier``: the scored
    closure is f = 1 (dividing bias on Bassett 2001, K straight +0.26, K branch
    +0.23), and the harness reports, never recommends, other values.

    Do not add a ``BorderCarnotLossElement`` on a port: the closure already
    carries the turning loss.
    """

    def __init__(
        self,
        id: str,
        inlet_nodes: list[str],
        outlet_nodes: list[str],
        inlet_angles_deg: list[float] | None = None,
        outlet_angles_deg: list[float] | None = None,
        port_areas: list[float] | None = None,
        loss_factor: float | list[float] = 1.0,
    ):
        super().__init__(
            id=id,
            inlet_nodes=inlet_nodes,
            outlet_nodes=outlet_nodes,
            inlet_angles_deg=inlet_angles_deg,
            outlet_angles_deg=outlet_angles_deg,
            port_areas=port_areas,
        )
        if self.N > IDEAL_MAX_PORTS:
            raise ValueError(
                f"IdealMomentumTeeElement '{id}': {self.N} ports, the kernel takes at most "
                f"{IDEAL_MAX_PORTS}."
            )
        factors = (
            [float(loss_factor)] * self.N
            if isinstance(loss_factor, int | float)
            else [float(f) for f in loss_factor]
        )
        if len(factors) != self.N:
            raise ValueError(
                f"IdealMomentumTeeElement '{id}': loss_factor has {len(factors)} entries "
                f"for {self.N} ports."
            )
        if not all(math.isfinite(f) and f >= 0.0 for f in factors):
            raise ValueError(
                f"IdealMomentumTeeElement '{id}': loss_factor must be finite and >= 0, "
                f"got {factors}."
            )
        #: Per-port receiver loss factor (inlets first, then outlets).
        self.loss_factor: list[float] = factors

    def axis_rad(self) -> list[float]:
        """Outward port axes, in the kernel's frame (see the class docstring)."""
        return [
            math.radians(float(a)) if s > 0.0 else math.pi - math.radians(float(a))
            for a, s in zip(self.port_angles_deg, self._port_signs, strict=True)
        ]

    def _kernel(
        self,
        states: list[NetworkMixtureState],
        p_jn: float,
        port_mdots: list[float],
        with_density_derivative: bool = True,
    ):
        geom = _core.IdealJunctionGeometry()
        geom.area = [float(a) for a in self.port_areas]
        geom.axis_rad = self.axis_rad()
        geom.port_sign = [float(s) for s in self._port_signs]
        if any(f != 1.0 for f in self.loss_factor):
            geom.loss_factor = list(self.loss_factor)
        p_static = [float(s.P) for s in states]
        rho = [float(s.density()) for s in states]
        # outer_mdot: the connecting elements' own flows; port flow = sign * outer.
        outer = [float(m) * float(s) for m, s in zip(port_mdots, self._port_signs, strict=True)]
        drho_dp = [
            r / p if p > 0.0 and with_density_derivative else 0.0
            for r, p in zip(rho, p_static, strict=True)
        ]
        return _core.ideal_junction_residuals_and_jacobian(
            p_static, rho, drho_dp, outer, float(p_jn), geom
        )

    def residuals(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        N = self.N
        kernel = self._kernel(states, P_jct, port_mdots)
        residuals = [float(v) for v in kernel.residual]
        # Column order fixed by include/ideal_junction.h: P[0..N-1],
        # outer_mdot[0..N-1], p_jn.
        seed_names = (
            [f"{self.port_nodes[i]}.P" for i in range(N)]
            + [f"{self._port_element_ids[i]}.m_dot" for i in range(N)]
            + [f"{self.id}.P_jct"]
        )
        jac: dict[int, dict[str, float]] = {}
        for row in range(N + 1):
            entries: dict[str, float] = {}
            for col, name in enumerate(seed_names):
                value = float(kernel.jacobian[row][col])
                if value != 0.0:
                    entries[name] = entries.get(name, 0.0) + value
            jac[row] = entries
        self._add_density_columns(states, P_jct, port_mdots, kernel, jac)
        return residuals, jac

    def _add_density_columns(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
        kernel: object,
        jac: dict[int, dict[str, float]],
    ) -> None:
        """d(row)/d(port T) and d(row)/d(port Y[k]) through the port density.

        A port node's T and Y are not unknowns: the solver derives them by
        mixing and relays a ``{port}.T`` / ``{port}.Y[k]`` entry through that
        chain. Without these the junction rows missed how a port's density
        follows its mixed temperature, measured at 1.2e-4 (scaled) on a hot
        forward merge and 8.7e-4 on a hot reversed one (#517).

        d(row)/d(rho_i) comes from the kernel without assuming its row form:
        a second call with d(rho)/dP = 0 leaves only the explicit P column, and
        the difference is d(row)/d(rho_i) * rho_i/P_i. Ideal gas then gives
        d(rho)/dT = -rho/T and d(rho)/dY_k = -rho W/W_k (W the mixture molar
        mass).
        """
        N = self.N
        explicit = self._kernel(states, P_jct, port_mdots, with_density_derivative=False)
        for i in range(N):
            P_i, T_i = float(states[i].P), float(getattr(states[i], "T", 0.0) or 0.0)
            if not (P_i > 0.0 and T_i > 0.0):
                continue
            rho_i = float(states[i].density())
            drho_dp = rho_i / P_i
            W = float(cb.mwmix(states[i].X))
            Y = list(states[i].Y)
            for row in range(N + 1):
                d_rho = (
                    float(kernel.jacobian[row][i]) - float(explicit.jacobian[row][i])
                ) / drho_dp
                if d_rho == 0.0:
                    continue
                entries = jac[row]
                name_T = f"{self.port_nodes[i]}.T"
                entries[name_T] = entries.get(name_T, 0.0) - d_rho * rho_i / T_i
                for k, y in enumerate(Y):
                    if y <= 0.0:
                        continue
                    name_Y = f"{self.port_nodes[i]}.Y[{k}]"
                    W_k = _species_molar_mass(k)
                    entries[name_Y] = entries.get(name_Y, 0.0) - d_rho * rho_i * W / W_k

    def seed_P_jct(self, pt_common: float, m_common: float, rho: float, common_i: int) -> float:
        """Initial ``p_jn`` from the common port's guessed Pt and flow.

        ``P_jct`` is a static pressure here. Where the common port supplies
        (a dividing tee) it sits AT ``p_jn``, so ``p_jn = Pt - m^2/(2 rho A^2)``;
        where it collects, the suppliers sit at ``p_jn`` and their heads are
        smaller than the common one by roughly the same amount, so the same
        subtraction is the better first guess there too. The other junction
        elements carry Pt in this slot and keep the solver's default.
        """
        A = float(self.port_areas[common_i] or 0.0)
        if not (A > 0.0 and rho > 0.0):
            return pt_common
        return pt_common - m_common * m_common / (2.0 * rho * A * A)

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        energy_rel_tol: float = 1e-9,
    ) -> bool:
        """Post-solve check that the junction did not create flow work.

        Unlike ``MultiPortChamberElement`` there is no direction check: a
        reversed port is a state this closure models, and the solver feeds a
        reversed junction port by the sign of its flow.

        The closure cannot create flow work: measured
        volumetrically relative to ``p_jn``, every exact root dissipates
        ``sum_in m_k |v_k|^2/2 - M |V|^2/2 + sum_out m_j |v_j - V|^2/2 >= 0``
        (Jensen on the mass-weighted mean, plus a sum of squares). What this
        guards is the SOLVE: a stalled state reported as converged, or a
        future change to the rows. Pt is used mass-weighted, as the other
        junctions do, because the solution dict carries no density.
        """
        flows: list[float] = []
        pts: list[float] = []
        for i in range(self.N):
            m_key = f"{self._port_element_ids[i]}.m_dot"
            pt_key = f"{self.port_nodes[i]}.Pt"
            if not self._port_element_ids[i] or m_key not in sol or pt_key not in sol:
                return True
            flows.append(float(self._port_signs[i]) * float(sol[m_key]))
            pts.append(float(sol[pt_key]))
        net_outflow_of_energy = sum(f * pt for f, pt in zip(flows, pts, strict=True))
        scale = sum(abs(f) for f in flows) * max(abs(pt) for pt in pts)
        return net_outflow_of_energy <= energy_rel_tol * scale

    def diagnostics(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Port bookkeeping, ``p_jn``, and the loss coefficients the closure produced.

        ``port_i_K`` is referenced to the common port's dynamic head (its own
        ``Pt - P``, what the port node calls dynamic head): ``(Pt_c - Pt_i)/q_c``
        on a dividing junction, ``(Pt_i - Pt_c)/q_c`` on a joining one -- the
        conventions of Bassett's K5/K6 and K11/K12. The common port of THAT
        definition is read off the CONVERGED flows (the lone port flowing the
        other way); where there is no lone port the K fields are omitted.

        For a three-port tee the named aliases are GEOMETRIC, as the GUI wires
        them: ``m_dot_com`` / ``m_dot_straight`` / ``m_dot_branch`` are the
        declared ports' flows, signed positive in the declared direction (a
        negative one ran reversed; ``n_reversed_ports`` counts them). The K
        aliases -- ``K_straight``/``K_branch`` dividing, ``K11``/``K12``
        joining -- and ``mass_flow_ratio`` are emitted only when the converged
        regime IS the declared one, so a tee that ran as a merge although wired
        as a split never reports a merge coefficient under a split's name.
        """
        diag = super().diagnostics(states, P_jct, port_mdots)
        diag["p_jn"] = float(P_jct)
        if port_mdots is None or len(port_mdots) != self.N:
            return diag
        outer = [float(s) * float(m) for s, m in zip(self._port_signs, port_mdots, strict=True)]
        diag["n_reversed_ports"] = float(sum(m < 0.0 for m in outer))
        declared_common = None
        if self.N == 3:
            declared_common = 0 if len(self.inlet_nodes) == 1 else self.N - 1
            d_straight, d_branch = (i for i in range(self.N) if i != declared_common)
            diag["m_dot_com"] = outer[declared_common]
            diag["m_dot_straight"] = outer[d_straight]
            diag["m_dot_branch"] = outer[d_branch]
        out = [i for i in range(self.N) if port_mdots[i] > 0.0]
        inn = [i for i in range(self.N) if port_mdots[i] < 0.0]
        if len(inn) == 1 and len(out) == self.N - 1:
            common, dividing = inn[0], True
        elif len(out) == 1 and len(inn) == self.N - 1:
            common, dividing = out[0], False
        else:
            return diag
        q_c = float(states[common].Pt) - float(states[common].P)
        if not q_c > 0.0:
            return diag
        m_c = abs(float(port_mdots[common]))
        others = [i for i in range(self.N) if i != common]
        for i in others:
            dpt = float(states[common].Pt) - float(states[i].Pt)
            diag[f"port_{i}_K"] = (dpt if dividing else -dpt) / q_c
            diag[f"port_{i}_flow_fraction"] = abs(float(port_mdots[i])) / m_c
        if declared_common is not None and common == declared_common:
            if dividing:
                diag["K_straight"] = diag[f"port_{d_straight}_K"]
                diag["K_branch"] = diag[f"port_{d_branch}_K"]
            else:
                diag["K11"] = diag[f"port_{d_straight}_K"]
                diag["K12"] = diag[f"port_{d_branch}_K"]
            diag["mass_flow_ratio"] = diag[f"port_{d_branch}_flow_fraction"]
        return diag
