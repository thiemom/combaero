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
no valve; it is the remainder when there is one. It should be the path that
needs the most pressure: a metered outlet needing more than the reference
gets would need a pump, and the solve is refused by the second-law check,
naming the outlet to make the reference.

Each outlet carries the inlet's Tt and composition; Pt follows each outlet's
own path.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .components import SPLITTER_FRACTION_SUM_TOL, SplitterOutletElement

if TYPE_CHECKING:
    from .graph import FlowNetwork

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
        The outlet node whose valve is open (``Pt_out = Pt_in``). Default: the
        remainder; required when there is none.
    areas : dict[str, float] or None
        Valve reference areas [m^2] per outlet node, for the reported
        ``valve_Cd`` and ``valve_zeta``. Outlets without one report the
        effective area ``valve_CdA`` only.
    """

    id: str
    inlet: str
    outlets: dict[str, float | None]
    reference: str | None = None
    areas: dict[str, float] | None = None

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
            if not rest:
                raise ValueError(
                    f"{where}: name the reference outlet, the one whose valve is open "
                    "(the path needing the most pressure); without a remainder there "
                    "is no default."
                )
            self.reference = rest[0]
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

    def add_to(self, net: FlowNetwork) -> list[SplitterOutletElement]:
        """Add the outlet elements to ``net`` (its inlet and outlet nodes must
        be in it already) and return them."""
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
            )
            for o, f in self.fractions.items()
        ]
        for e in elements:
            net.add_element(e)
        return elements

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The splitter's result from a solve's ``__element_diag__``: the
        throughput, and per outlet node its flow, fraction set and reached,
        and the valve's drop and metering."""
        outlets: dict[str, dict[str, Any]] = {}
        for o in self.outlets:
            d = element_diag.get(self.element_id(o), {})
            outlets[o] = {k: d[k] for k in _OUTLET_KEYS if k in d}
        m_in = math.fsum(float(v.get("m_dot", 0.0)) for v in outlets.values())
        return {
            "m_dot": m_in,
            "reference": self.reference,
            "needs_pump": sorted(o for o, v in outlets.items() if v.get("valve_needs_pump")),
            "outlets": outlets,
        }
