"""Flow splitter: one inlet to N outlets in given fractions (#903).

A rig divides one supply into metered fractions. ``FlowSplitter`` is that
division as one configuration, a mass-flow controller per outlet: outlet i
carries f_i of the splitter's throughput, with sum f_i = 1. One outlet may be
declared the remainder (``None``), 1 - sum of the others, so the sum cannot
drift; a set that does not sum to 1 is refused, never renormalised.

``FlowSplitter`` is the configuration, not a network element. It adds one
``SplitterOutletElement`` per outlet, each a two-port from the inlet node to
its outlet node, so the solver sees nothing new::

    inlet --{id}__{o1}--> o1      m_1 = f_1 M     (valve: dPt free)
          --{id}__{o2}--> o2      m_2 = f_2 M     (valve: dPt free)
          --{id}__{oN}--> oN      Pt_oN = Pt_in   (the reference: open)

with M = sum_j m_j. Each metered outlet's row is linear in the flows and
holds no pressure, so its total-pressure drop is free: the drop its valve
would have to take, reported with the effective area that meters it. Of N
fractions only N - 1 are independent, so one outlet, the **reference**, has
no valve. It has to be the path that needs the most pressure: a metered
outlet needing more than the reference gets would need a pump, a root the
second-law check refuses.

**The reference is found, not picked** (``auto_reference``, default True).
The solve starts from the user's ``reference``, else the remainder, else a
first guess (the path with the largest seed resistance at its fraction). If
the root is refused because a metered outlet would need a pump,
``NetworkSolver.solve`` moves the reference to the outlet that needs the most
pressure at that root and solves again from it (``solve_with_reference_swaps``):
an outer loop around the smooth Newton problem, at most N - 1 swaps per
splitter and never back to a reference already tried. At a fixed supply flow
one swap is exact (the outlet pressures do not depend on the reference); at a
fixed supply pressure the flows move and a second swap may follow. The swaps
are reported (``__reference_swaps__``, the message, each outlet's
``reference_swaps``). With no feasible reference, or ``auto_reference=False``,
the refusal stands and names the outlet to make the reference.

Each outlet carries the inlet's Tt and composition; Pt follows each outlet's
own path.
"""

from __future__ import annotations

import math
import time
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .components import (
    SPLITTER_FRACTION_SUM_TOL,
    NetworkNode,
    PressureBoundary,
    SplitterOutletElement,
)

if TYPE_CHECKING:
    from .graph import FlowNetwork
    from .solver import NetworkSolver

#: Density [kg/m^3], pressure [Pa] and flow [kg/s] the first reference guess
#: evaluates the outlet paths' seed resistances at. Only their ORDER is used,
#: and a guess only: a wrong one costs one swap.
_GUESS_RHO, _GUESS_P, _GUESS_MDOT = 1.2, 1.0e5, 0.1

#: The per-outlet diagnostics ``summarize`` collects, when present.
_OUTLET_KEYS = (
    "m_dot",
    "fraction",
    "fraction_achieved",
    "is_reference",
    "Pt_in",
    "Pt_out",
    "valve_dPt",
    "valve_CdA",
    "valve_Cd",
    "valve_zeta",
    "valve_choked",
    "valve_needs_pump",
)


@dataclass
class FlowSplitter:
    """One inlet node split to N outlet nodes in given mass-flow fractions.

    Parameters
    ----------
    id : str
        The splitter's id; its outlet elements are ``{id}__{outlet node}``.
    inlet : str
        The node the splitter draws from (a plenum, or a mass-flow or
        pressure boundary supplying it). It must be in the network.
    outlets : dict[str, float | None]
        Outlet node id -> fraction of the throughput, each in (0, 1). At most
        one ``None``: the remainder, 1 - the sum of the others. Without a
        remainder the fractions must sum to 1 within
        ``SPLITTER_FRACTION_SUM_TOL``; they are never renormalised.
    reference : str or None
        The outlet node whose valve is open (``Pt_out = Pt_in``) to start
        from. Default: the remainder, else the first guess (``add_to``).
        With ``auto_reference`` the solve moves it when it is infeasible.
    areas : dict[str, float] or None
        Valve reference areas [m^2] per outlet node, for the reported
        ``valve_Cd`` and ``valve_zeta``. Outlets without one report the
        effective area ``valve_CdA`` only.
    auto_reference : bool
        Move the reference to the outlet needing the most pressure when a
        solve is refused for a pump at a metered outlet (default True).
        False keeps the reference as given, and such a solve is refused.
    """

    id: str
    inlet: str
    outlets: dict[str, float | None]
    reference: str | None = None
    areas: dict[str, float] | None = None
    auto_reference: bool = True

    def __post_init__(self) -> None:
        where = f"FlowSplitter {self.id!r}"
        self.outlets = dict(self.outlets)
        if len(self.outlets) < 2:
            raise ValueError(f"{where}: needs at least two outlets, got {len(self.outlets)}")
        if self.inlet in self.outlets:
            raise ValueError(f"{where}: the inlet {self.inlet!r} is also an outlet")
        rest = [o for o, f in self.outlets.items() if f is None]
        if len(rest) > 1:
            raise ValueError(f"{where}: at most one remainder outlet (None), got {rest}")
        for o, f in self.outlets.items():
            if f is None:
                continue
            if not (isinstance(f, int | float) and math.isfinite(f) and 0.0 < f < 1.0):
                raise ValueError(
                    f"{where}: outlet {o!r} fraction must be in (0, 1), got {f!r}. An "
                    "outlet that carries nothing is no outlet: remove it."
                )
        given = math.fsum(float(f) for f in self.outlets.values() if f is not None)
        if rest:
            left = 1.0 - given
            if not left > SPLITTER_FRACTION_SUM_TOL:
                raise ValueError(
                    f"{where}: the given fractions sum to {given!r}, leaving the "
                    f"remainder {rest[0]!r} nothing; they are not renormalised."
                )
        elif not abs(given - 1.0) <= SPLITTER_FRACTION_SUM_TOL:
            raise ValueError(
                f"{where}: the fractions sum to {given!r}, not 1 (tolerance "
                f"{SPLITTER_FRACTION_SUM_TOL:g}); they are not renormalised. Declare "
                "one outlet the remainder (None) to take what the others leave."
            )
        if self.reference is None:
            # Without a remainder the first guess is made in add_to, which
            # sees the paths.
            self.reference = rest[0] if rest else None
        elif self.reference not in self.outlets:
            raise ValueError(f"{where}: reference {self.reference!r} is not one of its outlets")
        self.areas = dict(self.areas or {})
        for o, a in self.areas.items():
            if o not in self.outlets:
                raise ValueError(f"{where}: area given for {o!r}, which is not an outlet")
            if not (math.isfinite(a) and a > 0.0):
                raise ValueError(f"{where}: outlet {o!r} area must be positive, got {a!r}")

    @property
    def fractions(self) -> dict[str, float]:
        """Outlet node -> its fraction, the remainder resolved."""
        given = math.fsum(float(f) for f in self.outlets.values() if f is not None)
        return {o: (1.0 - given if f is None else float(f)) for o, f in self.outlets.items()}

    def element_id(self, outlet: str) -> str:
        """The id of the element feeding ``outlet``."""
        return f"{self.id}__{outlet}"

    def first_guess(self, net: FlowNetwork) -> str:
        """The outlet to start from when none is given and there is no
        remainder: the one whose path needs the most pressure by its seed
        resistance K at its fraction, K f^2 (each outlet node feeding one
        element, the sinks at one pressure). Else the first outlet."""
        first = next(iter(self.outlets))
        need: dict[str, float] = {}
        sinks: set[float] = set()
        for o, f in self.fractions.items():
            downs = net.get_downstream_elements(o) if o in net.nodes else []
            if len(downs) != 1:
                return first
            try:
                k = downs[0]._seed_resistance(_GUESS_RHO, _GUESS_P, _GUESS_MDOT)
            except (ValueError, RuntimeError, TypeError, AttributeError):
                return first
            if k is None or not math.isfinite(k):
                return first
            sink: NetworkNode | None = net.nodes.get(downs[0].to_node)
            if isinstance(sink, PressureBoundary):
                sinks.add(float(sink.Pt))
            need[o] = float(k) * f * f
        if len(sinks) > 1:
            return first
        return max(need, key=lambda o: need[o])

    def add_to(self, net: FlowNetwork) -> list[SplitterOutletElement]:
        """Add the outlet elements to ``net`` (its inlet and outlet nodes must
        be in it already) and return them. Without a given reference or a
        remainder, the start is ``first_guess``."""
        if self.reference is None:
            self.reference = self.first_guess(net)
        ids = [self.element_id(o) for o in self.outlets]
        elements = [
            SplitterOutletElement(
                self.element_id(o),
                from_node=self.inlet,
                to_node=o,
                fraction=f,
                outlets=ids,
                reference=(o == self.reference),
                area=self.areas.get(o) if self.areas else None,
                splitter=self.id,
                auto_reference=self.auto_reference,
            )
            for o, f in self.fractions.items()
        ]
        for e in elements:
            net.add_element(e)
        return elements

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The splitter's result from a solve's ``__element_diag__``: the
        throughput, the reference the solve ended with (``reference``; the
        one it started from is ``reference_start``), and per outlet node its
        flow, fraction set and reached, and the valve's drop and metering."""
        outlets: dict[str, dict[str, Any]] = {}
        for o in self.outlets:
            d = element_diag.get(self.element_id(o), {})
            outlets[o] = {k: d[k] for k in _OUTLET_KEYS if k in d}
        m_in = math.fsum(float(v.get("m_dot", 0.0)) for v in outlets.values())
        used = [o for o, v in outlets.items() if v.get("is_reference")]
        return {
            "m_dot": m_in,
            "reference": used[0] if len(used) == 1 else self.reference,
            "reference_start": self.reference,
            "needs_pump": sorted(o for o, v in outlets.items() if v.get("valve_needs_pump")),
            "outlets": outlets,
        }


# ---------------------------------------------------------------------------
# The reference, found by the solve (#903)
# ---------------------------------------------------------------------------


def auto_reference_groups(net: FlowNetwork) -> dict[str, list[SplitterOutletElement]]:
    """{splitter id: its outlet elements} for the splitters whose reference
    the solve may move (``auto_reference``)."""
    groups: dict[str, list[SplitterOutletElement]] = {}
    for e in net.elements.values():
        if isinstance(e, SplitterOutletElement) and e.auto_reference:
            groups.setdefault(e.splitter, []).append(e)
    return groups


def _pumping(
    group: list[SplitterOutletElement], sol: dict[str, Any], at_root: bool = True
) -> list[str]:
    """The group's metered outlets whose stream would have to gain total
    pressure: at a root, those the second-law check refused; at the best
    iterate of a solve that did not converge (a wrong reference can stall a
    cold start before it reaches the refused root), every one flowing
    against its drop there."""
    failed = set(sol.get("__second_law_failures__") or [])
    out = []
    for e in group:
        if e.reference or (at_root and e.id not in failed):
            continue
        m = float(sol.get(f"{e.id}.m_dot", 0.0))
        drop = float(sol.get(f"{e.from_node}.Pt", 0.0)) - float(sol.get(f"{e.to_node}.Pt", 0.0))
        if m * drop < 0.0:  # flow against the drop: a pump
            out.append(e.id)
    return out


def _set_reference(group: list[SplitterOutletElement], eid: str) -> None:
    for e in group:
        e.reference = e.id == eid


def _attempt(solver: NetworkSolver, kw: dict[str, Any]) -> tuple[dict[str, Any], list[Any]]:
    """One solve, its warnings held back: only the attempt returned speaks."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sol = solver.solve(**kw)
    return sol, list(caught)


def _raised_start(
    solver: NetworkSolver, sol: dict[str, Any], raise_to: dict[str, float]
) -> list[float] | None:
    """The refused root as a start, each swapped splitter's inlet node at the
    total pressure its new reference needs (its P too: an inlet node is at
    rest)."""
    x = sol.get("__x_solution__")
    if x is None:
        return None
    x = list(x)
    index = {name: i for i, name in enumerate(solver.unknown_names)}
    for node, pt in raise_to.items():
        for key in (f"{node}.P", f"{node}.Pt"):
            if key in index and math.isfinite(pt):
                x[index[key]] = pt
    return x


def solve_with_reference_swaps(solver: NetworkSolver, kw: dict[str, Any]) -> dict[str, Any]:
    """Solve, moving each auto-reference splitter's reference until no
    metered outlet needs a pump (``NetworkSolver.solve`` calls this when the
    network holds such a splitter).

    An outer loop around the smooth problem, never a switch inside Newton.
    After a root refused only for pumps at metered outlets, each such
    splitter's reference moves to the outlet whose node needs the most total
    pressure at that root, and the network is solved again from that root
    (cold if the warm start fails). Per splitter at most N - 1 swaps, never
    back to a reference already tried. A solve that then succeeds is
    returned with the swaps reported; otherwise the first refusal is
    returned, with the references restored and the swaps tried listed, so
    its message still names the outlet to make the reference.

    Reported: ``__reference_swaps__`` = {splitter: {"reference", "start",
    "swapped_from" (the references left, in order)}}, ``__reference_swap_solves__``
    (solves run, the first included), ``__reference_swap_time__`` [s], and
    on each outlet's diagnostics ``reference_swaps`` (swaps of its splitter).
    """
    groups = auto_reference_groups(solver.network)
    start = {sid: next(e.id for e in g if e.reference) for sid, g in groups.items()}
    t0 = time.time()
    first, first_warn = _attempt(solver, kw)
    sol, caught = first, first_warn
    n_solves = 1
    tried = {sid: {ref} for sid, ref in start.items()}
    left: dict[str, list[str]] = {sid: [] for sid in groups}
    while not sol.get("__success__", False):
        converged = bool(sol.get("__converged__", False))
        pumps = {sid: _pumping(g, sol, at_root=converged) for sid, g in groups.items()}
        failed = set(sol.get("__second_law_failures__") or []) if converged else set()
        others = failed - {eid for p in pumps.values() for eid in p}
        if not any(pumps.values()) or others:
            break  # refused for something a reference cannot fix
        moved = False
        raise_to: dict[str, float] = {}
        for sid, p in pumps.items():
            g = groups[sid]
            if not p or len(left[sid]) >= len(g) - 1:
                continue
            need = {e.id: float(sol.get(f"{e.to_node}.Pt", -math.inf)) for e in g}
            new = max(need, key=lambda k: need[k])
            if new in tried[sid]:
                continue
            old = next(e.id for e in g if e.reference)
            _set_reference(g, new)
            # The start: the refused root with the inlet raised to the new
            # reference's pressure, which is its row (at a fixed supply flow,
            # the root itself).
            raise_to[g[0].from_node] = need[new]
            tried[sid].add(new)
            left[sid].append(old)
            moved = True
        if not moved:
            break
        timeout = kw.get("timeout")
        rest = None if timeout is None else max(timeout - (time.time() - t0), 1.0)
        warm = {**kw, "timeout": rest, "x0": _raised_start(solver, sol, raise_to)}
        sol, caught = _attempt(solver, warm)
        n_solves += 1
        if not sol.get("__converged__", False):
            sol, caught = _attempt(solver, {**kw, "timeout": rest})
            n_solves += 1
    names = {e.id: e.to_node for g in groups.values() for e in g}
    swapped = {sid: v for sid, v in left.items() if v}
    if swapped and not sol.get("__success__", False):
        # No feasible reference: the first refusal, its references restored.
        tried_text = "; ".join(
            f"'{sid}': "
            + " -> ".join(names[r] for r in [*v, next(e.id for e in groups[sid] if e.reference)])
            for sid, v in swapped.items()
        )
        for sid, g in groups.items():
            _set_reference(g, start[sid])
        sol, caught = first, first_warn
        what = (
            "found no feasible reference"
            if sol.get("__converged__", False)
            else "did not reach a solution"
        )
        sol["__message__"] = (
            f"{sol.get('__message__', '')} Automatic reference swaps {what} (tried {tried_text})."
        )
        swapped = {}
    elif swapped:
        sol["__message__"] = f"{sol.get('__message__', '')} " + " ".join(
            f"Flow splitter '{sid}': reference moved from '{names[v[0]]}' to "
            f"'{names[next(e.id for e in groups[sid] if e.reference)]}', the outlet "
            "needing the most pressure."
            for sid, v in swapped.items()
        )
    for w in caught:
        warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
    sol["__reference_swaps__"] = {
        sid: {
            "reference": names[next(e.id for e in g if e.reference)],
            "start": names[start[sid]],
            "swapped_from": [names[r] for r in swapped.get(sid, [])],
        }
        for sid, g in groups.items()
    }
    sol["__reference_swap_solves__"] = n_solves
    sol["__reference_swap_time__"] = time.time() - t0
    diag = sol.get("__element_diag__") or {}
    for sid, g in groups.items():
        for e in g:
            if e.id in diag:
                diag[e.id]["reference_swaps"] = float(len(swapped.get(sid, [])))
    return sol
