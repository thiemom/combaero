"""Second-law check per element and a control-volume ledger (#545).

Every converged solve is checked twice after the fact:

* **Per element: no flow work manufactured.** Over each element's control
  volume the entropy generation must not be negative::

      sigma = sum_out m s - sum_in m s - integral dQ / T  >=  0

  with ``s = cb.s_mass(Tt, X, Pt)`` (stagnation state; it includes the ideal
  mixing term) and the heat taken at the fluid's static temperature, which
  makes ``sigma`` the INTERNAL entropy generation -- a necessary condition
  that needs no wall temperature. Both sums and the heat integral are C++
  (``cb.entropy_generation``, ``cb.heat_entropy_flux``). For an adiabatic
  two-port of fixed composition sigma is exactly ``-m R ln(Pt_out / Pt_in)``,
  compressible flow included: "total pressure does not rise". Heating at Mach
  M costs total pressure and cooling may raise it (Rayleigh); the heat term
  allows exactly that, through the static temperature.

  Inside an element only its end Mach numbers are known (its own ``mach_in``
  / ``mach_out`` where it reports them), so the heat is taken along a linear
  Mach path and bracketed by constant-Mach paths at either end:
  ``sigma_norm_low <= sigma_norm <= sigma_norm_high``. ``sigma_norm = sigma /
  (R sum_in m)`` is dimensionless; for a two-port it is ``-ln(Pt_out /
  Pt_in)``. A verdict is ``"pass"``, ``"fail"`` (``sigma_norm_high <
  -SIGMA_TOL``: even the most favourable path destroys entropy) or
  ``"not_checked"`` with the reason -- never a silent pass.

* **The network: mass and energy close over its control volume.** The
  streams crossing at the flow boundaries (``NetworkNode.is_flow_boundary``),
  the heat given to the network, and the heat leaving it (``Q_wall_out``, a
  wall's heat landing on a boundary) or withheld (``Q_withheld``, heat a
  stagnant node cannot take). Enthalpy is absolute (formation included), so a
  combustor closes the balance as it is.

Both are REPORTS in this release: nothing is demoted. See #545 for the
rollout to rejection once the check has proven itself.

Control volumes
---------------
* A two-port element. The stream enters at its upstream node's state and
  leaves at the downstream node's total pressure (a merge chamber's main
  inlet: at the chamber's main face, which is what it delivers), with the upstream
  composition and the enthalpy its walls added (the solver deposits wall
  heat on the downstream node, so the node's own temperature is a mix).
* A multi-port junction. Each port's stream at its port node's state.
* A chain of crossflow segments (a liner's backside duct). A segment's
  momentum spans half of each station, so one segment is not a control
  volume; the chain is. It is built from the segments' ``prev_seg`` /
  ``next_seg`` links; its stations are the nodes they declare
  (``from_kappa`` / ``to_kappa``), less any boundary. Streams crossing it:
  the chain's two ends, and every other element at a station (the holes), at
  the station's total pressure. Heat on a station node is taken as raising
  the node's mixed stream to the node's temperature.

Not checked: elements without a solved flow, the legacy
``TeeJunctionElement``, and a stagnant control volume (no flow through it).
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import combaero as cb

#: Below this, sigma_norm counts as negative. Chosen from a sweep of 238 element
#: instances across the suite's networks plus 50 random junctions: the most
#: negative best-case sigma_norm_high seen was 0 (a lossless connection sits
#: exactly at 0) and the smallest real loss was 1.5e-5 (an iterative junction),
#: so nothing came within three orders of magnitude of this bound from below.
#: 1e-8 is 100x tighter than the former 1e-6 -- it now catches a sub-ppm
#: manufactured total pressure -- while keeping a ~10x margin over the ~1e-9
#: solver round-off. Global, not per element: no element class showed a noise
#: floor that would need a looser per-element bound (#545).
SIGMA_TOL = 1e-8


@dataclass
class _Stream:
    m: float  # > 0
    Tt: float
    Pt: float
    Y: Any


@dataclass
class _Volume:
    inflows: list[_Stream] = field(default_factory=list)
    outflows: list[_Stream] = field(default_factory=list)
    #: Entropy the heat brought in [W/K], one entry per heat input, as
    #: (Mach linear along the path, lowest, highest): inside an element only
    #: the end Mach numbers are known, so the path is bracketed.
    heat: list[tuple[float, float, float]] = field(default_factory=list)


def _X(Y: Any) -> np.ndarray:
    return np.asarray(cb.mass_to_mole(list(Y)), dtype=float)


def _sigma(cv: _Volume) -> dict[str, float]:
    """sigma [W/K], sigma_norm [-] and its band over the heat's Mach path, in C++."""

    def streams(sts: list[_Stream]) -> list[Any]:
        return [cb.MassStream(st.m, st.Tt, st.Pt, list(st.Y)) for st in sts]

    s_in, s_out = streams(cv.inflows), streams(cv.outflows)
    S_lin = sum(h[0] for h in cv.heat)
    S_low = sum(min(h) for h in cv.heat)
    S_high = sum(max(h) for h in cv.heat)
    central = cb.entropy_generation(s_in, s_out, S_lin)
    # Less heat entropy, more generation: the least favourable path is S_high.
    best = cb.entropy_generation(s_in, s_out, S_low) if cv.heat else central
    worst = cb.entropy_generation(s_in, s_out, S_high) if cv.heat else central
    return {
        "sigma": float(central.sigma),
        "sigma_norm": float(central.sigma_norm),
        "sigma_norm_low": float(worst.sigma_norm),
        "sigma_norm_high": float(best.sigma_norm),
    }


class SecondLawReport:
    """Evaluates both checks on a solver's converged vector ``x``."""

    def __init__(self, solver: Any, x: np.ndarray, sol: dict[str, Any]) -> None:
        self.solver = solver
        self.net = solver.network
        self.x = x
        self.sol = sol
        self._states: dict[str, Any] = {}

    # -- states and flows -------------------------------------------------

    def _state(self, nid: str) -> Any:
        if nid not in self._states:
            self._states[nid] = self.solver._get_node_state(self.net.nodes[nid], self.x)
        return self._states[nid]

    def _flow(self, eid: str) -> float | None:
        idx = self.solver._unknown_indices.get(eid, [])
        if not idx or not self.solver.unknown_names[idx[0]].endswith(".m_dot"):
            return None
        return float(self.x[idx[0]])

    def _mach(self, nid: str) -> float:
        return float(self.sol.get(f"{nid}.mach", 0.0) or 0.0)

    def _wall_heat_into(self, eid: str) -> float:
        """Heat the thermal walls give element ``eid`` [W] (Q > 0 runs A -> B)."""
        q = 0.0
        for wid, wall in self.net.walls.items():
            Q = self.sol.get(f"{wid}.Q")
            if Q is None:
                continue
            if wall.element_a == eid:
                q -= float(Q)
            if wall.element_b == eid:
                q += float(Q)
        return q

    def _element_mach(self, eid: str, up: str, dn: str, forward: bool) -> tuple[float, float]:
        """The duct's Mach at the element's inlet and outlet, in flow order.

        An element that reports ``mach_in`` / ``mach_out`` knows its duct's;
        the end nodes may be plenums at Mach 0.
        """
        diag = (self.sol.get("__element_diag__") or {}).get(eid) or {}
        m_from = diag.get("mach_in")
        m_to = diag.get("mach_out")
        if m_from is None or m_to is None:
            return self._mach(up), self._mach(dn)
        a, b = abs(float(m_from)), abs(float(m_to))
        return (a, b) if forward else (b, a)

    def _heat_entropy(
        self, m: float, T0_in: float, T0_out: float, Y: Any, M_in: float, M_out: float
    ) -> tuple[float, float, float]:
        Y = list(Y)
        lin = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_in, M_out))
        if M_in == M_out:
            return lin, lin, lin
        a = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_in, M_in))
        b = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_out, M_out))
        return lin, min(a, b, lin), max(a, b, lin)

    def _stream_at(self, nid: str, m: float, Pt: float | None = None) -> _Stream:
        st = self._state(nid)
        return _Stream(m, float(st.Tt), float(st.Pt if Pt is None else Pt), st.Y)

    def _heated(self, inflow: _Stream, Q: float) -> float:
        """Total temperature after ``Q`` [W] is added to ``inflow``."""
        if not Q:
            return inflow.Tt
        X = _X(inflow.Y)
        h = float(cb.h_mass(inflow.Tt, X)) + Q / inflow.m
        T = float(cb.calc_T_from_h_mass(h, X, inflow.Tt))
        # The inverse does not refuse an enthalpy below its range: it returned
        # -2.4e25 K for a wall taking 3.5 MJ/kg from a stream holding 1.45.
        if not (math.isfinite(T) and T > 0.0):
            raise ValueError(
                f"the heat ({Q:.4g} W) takes more than the stream's enthalpy "
                f"({inflow.m:.4g} kg/s at {inflow.Tt:.1f} K)"
            )
        return T

    def _outlet_Pt(self, eid: str, dn: str, m: float, forward: bool) -> float:
        """Total pressure the element delivers at its outlet node.

        A merge chamber's declared main inlet arrives at the chamber's main
        FACE, not at its mixed outlet state (#471): the side jets can lift the
        chamber's Pt above the main duct's own inlet Pt without the duct
        manufacturing anything, so closing the duct's control volume at the
        chamber Pt would charge it with the jets' momentum. The face is the
        one the solver hands the element (NetworkSolver._main_face).
        """
        node = self.net.nodes[dn]
        st = self._state(dn)
        if (
            forward
            and getattr(node, "main_inlet", None) == eid
            and getattr(node, "has_side_streams", lambda: False)()
        ):
            face_in = copy.copy(st)
            face_in.m_dot = m
            return float(node.main_face_state(face_in, m_main=m)[1])
        return float(st.Pt)

    # -- control volumes ---------------------------------------------------

    def _two_port(self, eid: str, element: Any) -> _Volume | str:
        m = self._flow(eid)
        if m is None:
            return "no solved flow"
        up, dn = (
            (element.from_node, element.to_node)
            if m >= 0.0
            else (
                element.to_node,
                element.from_node,
            )
        )
        m = abs(m)
        if m < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        forward = up == element.from_node
        inflow = self._stream_at(up, m)
        Q = self._wall_heat_into(eid)
        Tt_out = self._heated(inflow, Q)
        outflow = _Stream(m, Tt_out, self._outlet_Pt(eid, dn, m, forward), inflow.Y)
        cv = _Volume([inflow], [outflow])
        if Q:
            M_in, M_out = self._element_mach(eid, up, dn, forward)
            cv.heat.append(self._heat_entropy(m, inflow.Tt, Tt_out, inflow.Y, M_in, M_out))
        return cv

    def _junction(self, element: Any) -> _Volume | str:
        cv = _Volume()
        for i, pid in enumerate(element.port_nodes):
            outer = element._port_element_ids[i]
            m = self._flow(outer) if outer else None
            if m is None:
                return f"port {pid} has no solved flow"
            m_out = float(element._port_signs[i]) * m  # > 0 leaves the junction
            if m_out == 0.0:
                continue
            st = self._stream_at(pid, abs(m_out))
            (cv.outflows if m_out > 0.0 else cv.inflows).append(st)
        if sum(st.m for st in cv.inflows) < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        return cv

    def _chain(self, members: list[str]) -> _Volume | str:
        """A crossflow chain: its stations are the nodes its members declare."""
        els = self.net.elements
        stations = self._stations(members)
        cv = _Volume()
        for eid, element in els.items():
            touches = [n for n in (element.from_node, element.to_node) if n in stations]
            if not touches:
                continue
            m = self._flow(eid)
            if m is None:
                return f"element {eid} at a station has no solved flow"
            if eid in members:
                # a chain end: the member's far node is outside the chain
                outside = [n for n in (element.from_node, element.to_node) if n not in stations]
                if not outside:
                    continue
                far = outside[0]
                into_chain = (m >= 0.0) == (element.from_node == far)
                st_node = far if into_chain else touches[0]
                stream = self._stream_at(st_node, abs(m), Pt=self._state(far).Pt)
                if not into_chain:
                    # its wall heat lands on the far node: add it on the way out
                    Q = self._wall_heat_into(eid)
                    if Q:
                        T0_in = stream.Tt
                        stream.Tt = self._heated(stream, Q)
                        forward = element.to_node == far
                        M_in, M_out = self._element_mach(eid, touches[0], far, forward)
                        cv.heat.append(
                            self._heat_entropy(stream.m, T0_in, stream.Tt, stream.Y, M_in, M_out)
                        )
            else:
                # a hole (or any other element) at a station
                station = touches[0]
                into_chain = (m >= 0.0) == (element.to_node == station)
                src = element.from_node if m >= 0.0 else element.to_node
                stream = self._stream_at(
                    src if into_chain else station, abs(m), Pt=self._state(station).Pt
                )
            if stream.m:
                (cv.inflows if into_chain else cv.outflows).append(stream)
        for nid in stations:
            node = self.net.nodes[nid]
            q = sum(float(eb.Q) for eb in getattr(node, "energy_boundaries", None) or [])
            q -= float(self.sol.get(f"{nid}.Q_withheld", 0.0))
            m_thru = sum(m for m, _, _ in self.solver._node_inflows(nid, self.x))
            if q and m_thru > 0.0:
                # the node's heat raised its mixed stream to the node's Tt
                st = self._state(nid)
                X = _X(st.Y)
                T0_out = float(st.Tt)
                h_before = float(cb.h_mass(T0_out, X)) - q / m_thru
                T0_in = float(cb.calc_T_from_h_mass(h_before, X, T0_out))
                M = self._mach(nid)
                cv.heat.append(self._heat_entropy(m_thru, T0_in, T0_out, st.Y, M, M))
        if sum(st.m for st in cv.inflows) < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        return cv

    def _stations(self, members: list[str]) -> set[str]:
        """Nodes a member declares a station (``from_kappa`` / ``to_kappa``).

        A boundary is never inside a control volume: an impingement array's
        last segment declares its exit (a PressureBoundary) a station whose
        onward flow is zero, and the stream into it leaves the chain.
        """
        out: set[str] = set()
        for eid in members:
            e = self.net.elements[eid]
            if e.from_kappa is not None:
                out.add(e.from_node)
            if e.to_kappa is not None:
                out.add(e.to_node)
        return {nid for nid in out if not self.net.nodes[nid].is_boundary}

    def _chains(self) -> list[list[str]]:
        from .components import CrossflowSegmentElement

        segs = {
            eid: e for eid, e in self.net.elements.items() if isinstance(e, CrossflowSegmentElement)
        }
        seen: set[str] = set()
        chains: list[list[str]] = []
        for eid in segs:
            if eid in seen:
                continue
            chain, todo = [], [eid]
            while todo:
                cur = todo.pop()
                if cur in seen or cur not in segs:
                    continue
                seen.add(cur)
                chain.append(cur)
                todo += [n for n in (segs[cur].prev_seg, segs[cur].next_seg) if n]
            chains.append(sorted(chain))
        return chains

    # -- reports -----------------------------------------------------------

    def element_verdicts(self) -> dict[str, dict[str, Any]]:
        from .components import MultiPortChamberBase, TeeJunctionElement

        out: dict[str, dict[str, Any]] = {}

        def record(eids: list[str], build: Any, group: str | None) -> None:
            # A REPORT must never take the solve down: a control volume that
            # cannot be evaluated (a converged but clamped state, heat beyond
            # what a stream holds) is not checked, and says why.
            try:
                cv = build()
                verdict_or_reason: dict[str, Any] | str = cv if isinstance(cv, str) else _sigma(cv)
            except (ValueError, RuntimeError) as exc:
                verdict_or_reason = f"cannot be evaluated: {exc}"
            if isinstance(verdict_or_reason, str):
                verdict: dict[str, Any] = {"verdict": "not_checked", "reason": verdict_or_reason}
            else:
                verdict = verdict_or_reason
                # Fail only where even the most favourable Mach path does.
                ok = verdict["sigma_norm_high"] >= -SIGMA_TOL
                verdict["verdict"] = "pass" if ok else "fail"
            if group:
                verdict["group"] = group
            for eid in eids:
                out[eid] = dict(verdict)

        chained: set[str] = set()
        for chain in self._chains():
            if not self._stations(chain):
                continue  # no side streams: each segment is a plain two-port
            chained.update(chain)
            record(chain, lambda chain=chain: self._chain(chain), group="+".join(chain))
        for eid, element in self.net.elements.items():
            if eid in chained:
                continue
            if isinstance(element, TeeJunctionElement):
                record([eid], lambda: "legacy TeeJunctionElement", None)
            elif isinstance(element, MultiPortChamberBase):
                record([eid], lambda element=element: self._junction(element), None)
            else:
                record([eid], lambda eid=eid, element=element: self._two_port(eid, element), None)
        return out

    def _wall_side_entropy(self, eid: str, Q_in: float) -> tuple[float, float, float]:
        """Entropy [W/K] the heat ``Q_in`` brings into element ``eid``'s stream,
        integrated along the stream's own temperature path (lin, low, high)."""
        element = self.net.elements.get(eid)
        if element is None:
            raise ValueError(f"'{eid}' is a node, not a flowing element")
        n_walls = sum(1 for w in self.net.walls.values() if eid in (w.element_a, w.element_b))
        if n_walls > 1:
            raise ValueError(f"'{eid}' carries {n_walls} walls; their heat paths are shared")
        m = self._flow(eid)
        if m is None or abs(m) < cb.MIXER_HEAT_MDOT_FLOOR:
            raise ValueError(f"'{eid}' has no solved flow")
        up, dn = (
            (element.from_node, element.to_node)
            if m >= 0.0
            else (element.to_node, element.from_node)
        )
        inflow = self._stream_at(up, abs(m))
        Tt_out = self._heated(inflow, Q_in)
        M_in, M_out = self._element_mach(eid, up, dn, up == element.from_node)
        return self._heat_entropy(abs(m), inflow.Tt, Tt_out, inflow.Y, M_in, M_out)

    def wall_verdicts(self) -> dict[str, dict[str, Any]]:
        """Each ThermalWall as its own control volume (#551).

        The heat leaves side A's stream and enters side B's, each at that
        stream's own temperature along its path, so the wall generates
        ``sigma_wall = S_A + S_B`` [W/K] with ``S = integral dQ / T`` (negative
        on the side giving heat). It is negative when heat runs from the
        stream that is colder on average to the hotter one -- a wall model
        with no capacity limit can do that (a hot stream cooled below the cold
        stream's inlet). Evaluating it at the temperatures the wall model sees
        instead would be vacuous: ``Q = U (T_A - T_B)`` makes ``Q (1/T_B -
        1/T_A)`` non-negative by construction.

        A necessary condition, not a sufficient one: the path means can stay
        ordered while the end temperatures cross. ``sigma_wall_norm`` is
        ``sigma_wall`` over the mean ``|S|`` of the two sides. A REPORT only.
        """
        out: dict[str, dict[str, Any]] = {}
        for wid, wall in self.net.walls.items():
            Q = self.sol.get(f"{wid}.Q")
            if Q is None:
                out[wid] = {"verdict": "not_checked", "reason": "no heat solved for this wall"}
                continue
            try:
                # Q > 0 runs A -> B.
                S_a = self._wall_side_entropy(wall.element_a, -float(Q))
                S_b = self._wall_side_entropy(wall.element_b, float(Q))
            except (ValueError, RuntimeError) as exc:
                out[wid] = {"verdict": "not_checked", "reason": f"cannot be evaluated: {exc}"}
                continue
            scale = 0.5 * (abs(S_a[0]) + abs(S_b[0]))
            if scale <= 0.0:
                out[wid] = {"verdict": "pass", "sigma_wall": 0.0, "sigma_wall_norm": 0.0}
                continue
            sigma = S_a[0] + S_b[0]
            # The most favourable paths: the largest entropy on both sides.
            sigma_high = S_a[2] + S_b[2]
            out[wid] = {
                "verdict": "pass" if sigma_high / scale >= -SIGMA_TOL else "fail",
                "sigma_wall": sigma,
                "sigma_wall_norm": sigma / scale,
                "sigma_wall_norm_low": (S_a[1] + S_b[1]) / scale,
                "sigma_wall_norm_high": sigma_high / scale,
            }
        return out

    def ledger(self) -> dict[str, Any]:
        """Mass and energy over the network's control volume."""
        rows: list[dict[str, Any]] = []

        def add(node: str, via: str, m_in: float, Tt: float, Y: Any) -> None:
            H = m_in * float(cb.h_mass(Tt, _X(Y)))
            rows.append({"node": node, "via": via, "m_dot_in": m_in, "H_in": H})

        for nid, node in self.net.nodes.items():
            if not node.is_flow_boundary:
                continue
            ups = self.net.get_upstream_elements(nid)
            downs = self.net.get_downstream_elements(nid)
            if ups and downs:
                # an injection: only its own stream crosses here
                Y_inj = node.Y if node.Y is not None else self.solver._default_Y
                add(nid, nid, float(node.m_dot), float(node.Tt), Y_inj)
                continue
            for element in {e.id: e for e in [*ups, *downs]}.values():
                ind = self.solver._unknown_indices.get(element.id, [])
                if not ind:
                    continue
                # signed source -> sink, as the element itself reports it here
                m = float(element.flow_at_node(nid, self.x, ind))
                at_source = nid in element.all_source_nodes()
                if at_source == (m >= 0.0):
                    st = self._state(nid)
                    add(nid, element.id, abs(m), float(st.Tt), st.Y)
                    continue
                # leaving, at the mixed specific enthalpy of what feeds it
                far = element.all_sink_nodes() if at_source else element.all_source_nodes()
                m_far = H_far = 0.0
                for src in far:
                    if src == nid:
                        continue
                    m_src = abs(float(element.flow_at_node(src, self.x, ind)))
                    st = self._state(src)
                    m_far += m_src
                    H_far += m_src * float(cb.h_mass(float(st.Tt), _X(st.Y)))
                H = abs(m) * H_far / m_far if m_far > 0.0 else 0.0
                rows.append({"node": nid, "via": element.id, "m_dot_in": -abs(m), "H_in": -H})

        q_given = 0.0
        for nid, node in self.net.nodes.items():
            for eb in getattr(node, "energy_boundaries", None) or []:
                if eb.id != f"_wall_{nid}":
                    q_given += float(eb.Q)
            q_given += float(self.sol.get(f"{nid}.Q_fraction", 0.0))
        q_out = sum(float(v) for k, v in self.sol.items() if k.endswith(".Q_wall_out"))
        q_withheld = sum(float(v) for k, v in self.sol.items() if k.endswith(".Q_withheld"))

        m_in = sum(r["m_dot_in"] for r in rows if r["m_dot_in"] > 0.0)
        m_out = -sum(r["m_dot_in"] for r in rows if r["m_dot_in"] < 0.0)
        H_in = sum(r["H_in"] for r in rows if r["m_dot_in"] > 0.0)
        H_out = -sum(r["H_in"] for r in rows if r["m_dot_in"] < 0.0)
        energy = H_in + q_given - H_out - q_out - q_withheld
        h_scale = max(abs(H_in), abs(H_out), sum(abs(r["H_in"]) for r in rows), 1.0)
        return {
            "streams": rows,
            "m_dot_in": m_in,
            "m_dot_out": m_out,
            "mass_closure": m_in - m_out,
            "mass_closure_rel": (m_in - m_out) / max(m_in, m_out, cb.MIXER_HEAT_MDOT_FLOOR),
            "H_in": H_in,
            "H_out": H_out,
            "Q_given": q_given,
            "Q_wall_out": q_out,
            "Q_withheld": q_withheld,
            "energy_closure": energy,
            "energy_closure_rel": energy / h_scale,
        }
