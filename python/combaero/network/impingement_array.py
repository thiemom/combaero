"""Impingement array, channel exit: one configuration, built from rows (#465).

The configuration Florschuetz, Truman and Metzger (1981) measured: a jet
plate fed from a supply plenum, impinging on a target, the spent air leaving
in ONE direction down the gap between them. That is all this models. A bypass
crossflow (#467) and spent air leaving through the target (#468) are
different configurations and get their own elements.

``ImpingementArray`` is the configuration, not a network element: it adds the
rows to a network, one ``ImpingementPlateElement`` and one
``ImpingementCrossflowElement`` per row, joined by internal crossflow nodes,
so the solver sees nothing new. Ids are namespaced under the array's own id::

    supply --+------------------+------------------+
             |                  |                  |
         {id}__p1           {id}__p2           {id}__pN
             |                  |                  |
         {id}__c1 --{id}__x1-- {id}__c2 -- ... -- {id}__cN --{id}__xN-- exit

``add_wall`` distributes one thermal wall to the hot side across the rows,
and ``summarize`` collects the per-row diagnostics under the array's id.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb

from .components import (
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    PlenumNode,
    ThermalWall,
    WallLayer,
)

if TYPE_CHECKING:
    from .graph import FlowNetwork


@dataclass
class ImpingementArray:
    """A jet-plate impingement array whose spent air exits down the channel.

    Parameters
    ----------
    id : str
        The array's id; every row element and node is namespaced under it.
    n_rows : int
        Spanwise rows, counted from the closed (upstream) end.
    d_jet, xn_d, yn_d, z_d, span, plate_thickness, pattern, correlation, Cd,
    Nu_multiplier :
        As ``ImpingementPlateElement``; every row shares them.
    roughness : float
        Of the target and jet plate in the gap, for the crossflow friction.
    """

    id: str
    n_rows: int
    d_jet: float
    xn_d: float
    yn_d: float
    z_d: float
    span: float
    plate_thickness: float
    pattern: Literal["inline", "staggered"] = "inline"
    correlation: str = "fixed"
    Cd: float = cb.FLORSCHUETZ_1981_DEFAULT_CD
    Nu_multiplier: float = 1.0
    roughness: float = 0.0
    _rows: list[tuple[str, str, str]] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        if int(self.n_rows) != self.n_rows or self.n_rows < 1:
            raise ValueError(f"ImpingementArray {self.id!r}: n_rows must be a positive integer")

    def plate_id(self, row: int) -> str:
        return f"{self.id}__p{row}"

    def segment_id(self, row: int) -> str:
        return f"{self.id}__x{row}"

    def node_id(self, row: int) -> str:
        return f"{self.id}__c{row}"

    def add_to(self, net: FlowNetwork, supply: str, exit: str) -> None:
        """Add the rows between ``supply`` (the jet plenum) and ``exit``."""
        n = int(self.n_rows)
        for i in range(1, n + 1):
            net.add_node(PlenumNode(self.node_id(i)))
        self._rows = []
        for i in range(1, n + 1):
            to = self.node_id(i + 1) if i < n else exit
            net.add_element(
                ImpingementPlateElement(
                    self.plate_id(i),
                    supply,
                    self.node_id(i),
                    d_jet=self.d_jet,
                    xn_d=self.xn_d,
                    yn_d=self.yn_d,
                    z_d=self.z_d,
                    span=self.span,
                    plate_thickness=self.plate_thickness,
                    pattern=self.pattern,
                    correlation=self.correlation,
                    Cd=self.Cd,
                    row=i,
                    Nu_multiplier=self.Nu_multiplier,
                )
            )
            net.add_element(
                ImpingementCrossflowElement(
                    self.segment_id(i),
                    self.node_id(i),
                    to,
                    length=self.xn_d * self.d_jet,
                    height=self.z_d * self.d_jet,
                    span=self.span,
                    roughness=self.roughness,
                )
            )
            self._rows.append((self.plate_id(i), self.segment_id(i), self.node_id(i)))

    def add_wall(
        self,
        net: FlowNetwork,
        wall_id: str,
        hot_element: str,
        layers: list[WallLayer],
        R_fouling: float = 0.0,
        array_is_side_a: bool = False,
    ) -> list[str]:
        """One wall between ``hot_element`` and the array's target, as one
        ``ThermalWall`` per row, each on that row's target footprint.

        The hot side is evaluated once per row at its own state, so every row
        sees the same hot-gas h and T_aw: a hot side that changes along the
        array needs its own segments. Returns the per-row wall ids.
        """
        ids = []
        for i in range(1, int(self.n_rows) + 1):
            plate = self.plate_id(i)
            a, b = (plate, hot_element) if array_is_side_a else (hot_element, plate)
            wid = f"{wall_id}__r{i}"
            net.add_wall(
                ThermalWall(
                    id=wid,
                    element_a=a,
                    element_b=b,
                    layers=list(layers),
                    contact_area=net.elements[plate].surface.area,
                    R_fouling=R_fouling,
                )
            )
            ids.append(wid)
        net.thermal_coupling_enabled = True
        return ids

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The array's own diagnostics from its rows' (``__element_diag__``).

        Scalars for the array as a whole, plus a per-row list for each row
        quantity. The pressure block is supply-to-exit.
        """
        n = int(self.n_rows)
        plates = [element_diag[self.plate_id(i)] for i in range(1, n + 1)]
        segments = [element_diag[self.segment_id(i)] for i in range(1, n + 1)]
        m_total = sum(p["m_dot"] for p in plates)
        out: dict[str, Any] = {
            "m_dot": float(m_total),
            "n_rows": float(n),
            "n_holes": float(sum(p["n_holes"] for p in plates)),
            "Pt_in": float(plates[0]["Pt_in"]),
            "P_out": float(segments[-1]["P_out"]),
            "dP": float(plates[0]["Pt_in"] - segments[-1]["P_out"]),
            "Re_j_min": float(min(p["Re_j"] for p in plates)),
            "Re_j_max": float(max(p["Re_j"] for p in plates)),
            "Gc_Gj_max": float(max(p["Gc_Gj"] for p in plates)),
            "Nu_mean": float(sum(p["Nu"] for p in plates) / n),
            "htc_mean": float(sum(p["htc"] for p in plates) / n),
            "Cd": float(plates[0]["Cd"]),
            "Cd_in_range": float(all(p["Cd_in_range"] for p in plates)),
            "surface_extrapolated": float(any(p["surface_extrapolated"] for p in plates)),
            "jet_flow_nonuniformity": float(
                max(p["m_dot"] for p in plates) / min(p["m_dot"] for p in plates)
                if min(p["m_dot"] for p in plates) > 0.0
                else float("inf")
            ),
        }
        for key in ("m_dot", "Re_j", "Gc_Gj", "Nu", "htc"):
            out[f"rows_{key}"] = [float(p[key]) for p in plates]
        return out
