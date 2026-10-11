"""Effusion liner: a duct-fed effusion wall as one configuration (#471).

A backside duct (an annulus between liner and casing) carries coolant from
``coolant_in`` to ``coolant_out``; its inner wall is an effusion plate that
bleeds coolant into the hot gas. Both sides have crossflow: the coolant runs
past the holes, the gas past their exits.

``EffusionLiner`` is the configuration, not a network element. It adds N
stations along the duct, each a ``PlenumNode`` (Pt = P: the node pressure is
the duct's static pressure, the drive a hole sees), joined by
``CrossflowSegmentElement`` s, and one ``EffusionPlateElement`` panel per
station into the gas node::

    coolant_in --s0-- b1 --s1-- b2 -- ... -- bN --sN-- coolant_out
                      |         |               |
                    panel1    panel2    ...   panelN     (into gas)

* Each station is a BLEED: the segments carry Bassett's separating
  straight-run momentum (kappa = 0.75, half per side; ``station_half_drop``).
* ``s0`` enters the duct from ``coolant_in`` as a reservoir (``entry_K``);
  ``sN`` dumps into ``coolant_out`` (static continuity).
* Each panel is DUCT-FED: its Cd is McGreehan-Schotsch with the supply-side
  U1/Vi from the duct's mean velocity at its station. Its wall is the #475
  plate wall; the coolant-side heat transfer is still Andrews' still-plenum
  correlation, flagged ``coolant_ht_plenum_assumed`` (crossflow-supply heat
  transfer is the next step, pending sources).

Ids are namespaced under the liner's id: ``{id}__b{i}``, ``{id}__s{i}``,
``{id}__p{i}``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb

from .components import CrossflowSegmentElement, EffusionPlateElement, PlenumNode

if TYPE_CHECKING:
    from .graph import FlowNetwork

# Bassett et al. (2001) measured the separating straight run at main/branch
# area ratios psi = F_C/F_B of 1 to 3. An effusion station bleeds through
# holes far smaller than the duct; K5 is psi-independent by derivation, so
# the station is used there, but flagged.
BASSETT_PSI_MEASURED_MAX = 3.0


@dataclass
class EffusionLiner:
    """A duct-fed effusion liner, N stations along the backside duct.

    Parameters
    ----------
    id : str
    n_segments : int
        Stations (panels) along the duct. More resolve the bleed's decline
        along the liner -- the first panel passes more than the last.
    length, width : float
        Liner panel length (along the backside flow) and width [m]. Each
        station's panel is ``length/n_segments`` by ``width``.
    duct_height : float
        Backside duct height [m]; its cross-section is ``duct_height*width``.
    hole_diameter, wall_thickness, pitch_x, pitch_y, angle_deg :
        As ``EffusionPlateElement``; every panel shares them.
    entry_K : float
        Loss on entering the duct from ``coolant_in`` (0.5 sharp).
    wall_conductivity, gas_film, gas_augmentation, turbulence_intensity,
    gas_heat_flux, internal_Nu_multiplier, edge_radius :
        As ``EffusionPlateElement``.
    """

    id: str
    n_segments: int
    length: float
    width: float
    duct_height: float
    hole_diameter: float
    wall_thickness: float
    pitch_x: float
    pitch_y: float
    angle_deg: float = 90.0
    entry_K: float = 0.5
    roughness: float = 0.0
    edge_radius: float = 0.0
    wall_conductivity: float = 20.0
    gas_film: Literal["none", "baldauf_sellers"] = "none"
    gas_augmentation: float = 1.0
    turbulence_intensity: float = 0.05
    gas_heat_flux: float = 0.0
    internal_Nu_multiplier: float = 1.0

    def __post_init__(self) -> None:
        if int(self.n_segments) != self.n_segments or self.n_segments < 1:
            raise ValueError(f"EffusionLiner {self.id!r}: n_segments must be a positive integer")
        for name in ("length", "width", "duct_height"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"EffusionLiner {self.id!r}: {name} must be positive")

    @property
    def duct_area(self) -> float:
        return self.duct_height * self.width

    @property
    def duct_Dh(self) -> float:
        return 2.0 * self.duct_height * self.width / (self.duct_height + self.width)

    def node_id(self, i: int) -> str:
        return f"{self.id}__b{i}"

    def segment_id(self, i: int) -> str:
        return f"{self.id}__s{i}"

    def panel_id(self, i: int) -> str:
        return f"{self.id}__p{i}"

    def add_to(
        self, net: FlowNetwork, coolant_in: str, coolant_out: str, gas: str | list[str]
    ) -> None:
        """Add the stations between ``coolant_in`` and ``coolant_out``, each
        panel discharging into ``gas`` (one node, or one per station for an
        axial gas-side gradient)."""
        n = int(self.n_segments)
        gas_nodes = [gas] * n if isinstance(gas, str) else list(gas)
        if len(gas_nodes) != n:
            raise ValueError(f"EffusionLiner {self.id!r}: give one gas node or {n}")
        dx = self.length / n
        bleed = cb.STATION_KAPPA_BLEED_BASSETT
        A, Dh = self.duct_area, self.duct_Dh
        for i in range(1, n + 1):
            net.add_node(PlenumNode(self.node_id(i)))
        # s0: reservoir entry into station 1.
        net.add_element(
            CrossflowSegmentElement(
                self.segment_id(0),
                coolant_in,
                self.node_id(1),
                length=0.5 * dx,
                area=A,
                Dh=Dh,
                roughness=self.roughness,
                entry_K=self.entry_K,
                to_kappa=bleed,
                next_seg=self.segment_id(1),
            )
        )
        for i in range(1, n + 1):
            last = i == n
            net.add_element(
                CrossflowSegmentElement(
                    self.segment_id(i),
                    self.node_id(i),
                    coolant_out if last else self.node_id(i + 1),
                    length=(0.5 if last else 1.0) * dx,
                    area=A,
                    Dh=Dh,
                    roughness=self.roughness,
                    from_kappa=bleed,
                    prev_seg=self.segment_id(i - 1),
                    to_kappa=None if last else bleed,
                    next_seg=None if last else self.segment_id(i + 1),
                )
            )
            net.add_element(
                EffusionPlateElement(
                    self.panel_id(i),
                    self.node_id(i),
                    gas_nodes[i - 1],
                    hole_diameter=self.hole_diameter,
                    wall_thickness=self.wall_thickness,
                    pitch_x=self.pitch_x,
                    pitch_y=self.pitch_y,
                    panel_length=dx,
                    panel_width=self.width,
                    angle_deg=self.angle_deg,
                    correlation="McGreehanSchotsch",
                    edge_radius=self.edge_radius,
                    internal_Nu_multiplier=self.internal_Nu_multiplier,
                    wall_conductivity=self.wall_conductivity,
                    gas_film=self.gas_film,
                    gas_augmentation=self.gas_augmentation,
                    turbulence_intensity=self.turbulence_intensity,
                    gas_heat_flux=self.gas_heat_flux,
                    crossflow_segments=(self.segment_id(i - 1), self.segment_id(i)),
                    crossflow_area=A,
                )
            )

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The liner's own diagnostics from its stations' (``__element_diag__``)."""
        n = int(self.n_segments)
        panels = [element_diag[self.panel_id(i)] for i in range(1, n + 1)]
        m_in = float(element_diag[self.segment_id(0)]["m_dot"])
        m_out = float(element_diag[self.segment_id(n)]["m_dot"])
        bleed = [float(p["m_dot"]) for p in panels]
        hole_area = math.pi * self.hole_diameter**2 / 4.0
        psi = [self.duct_area / (p["n_holes"] * hole_area) for p in panels]
        out: dict[str, Any] = {
            "m_dot": m_in,
            "m_dot_out": m_out,
            "m_dot_bleed": float(sum(bleed)),
            "bleed_fraction": float(sum(bleed) / m_in) if m_in else 0.0,
            "n_segments": float(n),
            "n_holes": float(sum(p["n_holes"] for p in panels)),
            "U1_over_Vi_max": float(max(p.get("U1_over_Vi", 0.0) for p in panels)),
            "crossflow_cd_beyond_8pct": float(
                any(p.get("crossflow_cd_beyond_8pct", 0.0) for p in panels)
            ),
            "crossflow_cd_degraded": float(
                any(p.get("crossflow_cd_degraded", 0.0) for p in panels)
            ),
            "bleed_nonuniformity": float(max(bleed) / min(bleed)) if min(bleed) > 0.0 else math.inf,
            "station_psi_min": float(min(psi)),
            "station_psi_beyond_bassett": float(min(psi) > BASSETT_PSI_MEASURED_MAX),
            "coolant_ht_plenum_assumed": 1.0,
            "is_ingesting": float(any(p.get("is_ingesting", 0.0) for p in panels)),
        }
        if all("T_wall_hot" in p for p in panels):
            out["T_wall_hot_max"] = float(max(p["T_wall_hot"] for p in panels))
            out["rows_T_wall_hot"] = [float(p["T_wall_hot"]) for p in panels]
        out["rows_m_dot_bleed"] = bleed
        out["rows_U1_over_Vi"] = [float(p.get("U1_over_Vi", 0.0)) for p in panels]
        out["rows_Cd"] = [float(p["Cd"]) for p in panels]
        out["rows_P_backside"] = [float(p["P_in"]) for p in panels]
        return out
