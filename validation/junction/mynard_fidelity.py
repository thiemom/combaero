"""Mynard & Valen-Sendstad 2015 against its own figures: fidelity and Ref3D.

This is the one place the junction closure meets its own paper. Two questions,
kept apart (docs/VALIDATION_POLICY.md):

**Fidelity -- does our port reproduce Mynard's model?** Judged on his red
"Unified0D" curves (Figs 4, 6-11), with the port in his configuration
(``FAITHFUL``: eta on, no joining alpha, no dividing-streamline recovery). A
miss is our bug. The red lines are splines through the model evaluated at
the CFD sample points -- on the steep Fig 6d curve the spline leaves the model
between samples -- so the check is made AT the samples: the model's value at
each Ref3D marker's abscissa must lie on the drawn red curve, read in the
marker's own pixel column. That reading needs no x calibration; the
abscissa itself is the design value the marker snaps to (``_DESIGN``), with
the snap residual reported.

**Ref3D -- what does each configuration do on Mynard's own reference?**
Ref3D is laminar CFD, Re 350-2400, flat inlet profile, junction corners
smoothed. It is the author's own data, but our PRODUCTION closure deliberately
departs from his model (eta off, Hager/Bassett recovery on, joining alpha),
for reasons measured on turbulent data. The score shows what those departures
cost or gain in his regime; it is not a target, and nothing here is tuned on
it -- his envelope is far from the full range.

Geometry conventions, resolved against the red curves and stated per flow
type in ``_junction`` (the panel registry says where a caption disagrees):
branch axes point away from the junction; areas are relative to the main
duct; the common branch carries unit flow.

Run:

    uv run python -m validation.junction.mynard_fidelity
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from combaero.network._mynard2010 import junction_loss_coefficient
from combaero.network.mpce_element import MultiPortChamberElement

_DATA = Path(__file__).resolve().parent / "data" / "mynard2015"

#: Mynard's published model (the Matlab JunctionLossCoefficient.m).
FAITHFUL = dict(joining_etransfer_alpha=0.0, eta_scale=1.0, dividing_streamline_recovery=0.0)

#: What MultiPortChamberElement evaluates by default.
PRODUCTION = dict(
    joining_etransfer_alpha=MultiPortChamberElement.DEFAULT_JOINING_ETRANSFER_ALPHA,
    eta_scale=MultiPortChamberElement.DEFAULT_ETA_SCALE,
)

# Design abscissae the CFD markers sit on. Flow fractions and angles are round
# numbers; area ratios are squares of diameter ratios k/6 (the insets print
# 1/psi = 0.25 and 0.44, i.e. (3/6)^2 and (4/6)^2).
_DESIGN = {
    "lam": np.round(np.arange(0.05, 1.0, 0.05), 6),
    "phi": np.arange(0.0, 181.0, 5.0),
    "area": np.array([(k / 6.0) ** 2 for k in range(1, 7)]),
}


@dataclass(frozen=True)
class PanelSpec:
    key: str
    flow_type: int  # Mynard's Fig 1 numbering; 34 = his Fig 4 Y-junction
    sweep: str  # "Re", "lam", "phi", "area"
    branch: str  # which coefficient, in our geometry: "side", "straight", "a"
    note: str = ""


PANELS: tuple[PanelSpec, ...] = (
    PanelSpec("fig04", 34, "phi", "a", "collector at the plotted angle; the other 90 deg from it"),
    *(
        PanelSpec(f"fig06{c}", 1, s, b)
        for c, s, b in zip(
            "abcdefgh",
            ["Re", "lam", "phi", "area"] * 2,
            ["side"] * 4 + ["straight"] * 4,
            strict=True,
        )
    ),
    *(PanelSpec(f"fig07{c}", 2, s, "a") for c, s in zip("abcd", ["Re", "lam", "phi", "area"], strict=True)),
    *(PanelSpec(f"fig08{c}", 3, s, "a") for c, s in zip("abc", ["Re", "lam", "phi"], strict=True)),
    *(
        PanelSpec(
            f"fig09{c}",
            4,
            s,
            b,
            "caption puts the side inlet on top; Mynard's own Unified0D, and the "
            "dead-branch limit K -> -1, put it on the bottom row",
        )
        for c, s, b in zip(
            "abcdefgh",
            ["Re", "lam", "phi", "area"] * 2,
            ["straight"] * 4 + ["side"] * 4,
            strict=True,
        )
    ),
    *(PanelSpec(f"fig10{c}", 5, s, "a") for c, s in zip("abcd", ["Re", "lam", "phi", "area"], strict=True)),
    *(PanelSpec(f"fig11{c}", 6, s, "a") for c, s in zip("abc", ["Re", "lam", "phi"], strict=True)),
)


def _junction(flow_type: int, lam: float, phi_deg: float | None, area: float) -> tuple[np.ndarray, ...]:
    """(Q, A, theta, slots) for one of Mynard's six flow types.

    ``slots`` maps a branch name to its index in the kernel's K output (one K
    per collector when diverging, per supplier when converging, in array
    order). Defaults are Mynard's: lam 0.5, 90 deg for a T and 45 for a Y,
    equal areas.
    """
    pi = math.pi
    d = math.radians(phi_deg if phi_deg is not None else (45.0 if flow_type in (3, 6) else 90.0))
    if flow_type == 1:
        # Diverging T. lam = side fraction; phi = angle between the side axis
        # and the straight outlet's; area = side / main.
        A = np.array([1.0, 1.0, area])
        Q = np.array([1.0, -(1.0 - lam), -lam])
        th = np.array([0.0, pi, pi - d])
        slots = {"straight": 0, "side": 1}
    elif flow_type == 2:
        # Side inlet feeding both straight outlets. lam = fraction to outlet
        # "a"; phi = flow deflection into it; area = inlet / outlet.
        A = np.array([area, 1.0, 1.0])
        Q = np.array([1.0, -lam, -(1.0 - lam)])
        th = np.array([pi - d, 0.0, pi])
        slots = {"a": 0}
    elif flow_type == 3:
        # Diverging Y. lam = fraction to outlet "a"; phi = each outlet's
        # deflection.
        A = np.ones(3)
        Q = np.array([1.0, -lam, -(1.0 - lam)])
        th = np.array([0.0, pi - d, pi + d])
        slots = {"a": 0}
    elif flow_type == 34:
        # Fig 4: Y with 90 deg between the outlets, turned; phi = outlet "a"'s
        # deflection, the other is deflected 90 - phi the other way.
        A = np.ones(3)
        Q = np.array([1.0, -0.5, -0.5])
        th = np.array([0.0, pi - d, pi + (pi / 2.0 - d)])
        slots = {"a": 0}
    elif flow_type == 4:
        # Converging T. lam = side fraction; phi = angle between the side
        # axis and the OUTLET's (so a small phi opposes the main flow -- the
        # same physical angle as type 1's, measured from the downstream
        # branch). The area sweep holds both inlet velocities equal, so the
        # side fraction follows its area: lam = area / (1 + area).
        if area != 1.0:
            lam = area / (1.0 + area)
        A = np.array([1.0, 1.0, area])
        Q = np.array([-1.0, 1.0 - lam, lam])
        th = np.array([0.0, pi, -d])
        slots = {"straight": 0, "side": 1}
    elif flow_type == 5:
        # Both straight inlets into the side. lam = inlet "a"'s fraction;
        # area = outlet / inlet; Fig 10d holds both inlet Re equal (lam 0.5).
        A = np.array([1.0, 1.0, area])
        Q = np.array([lam, 1.0 - lam, -1.0])
        th = np.array([pi, 0.0, -d])
        slots = {"a": 0}
    elif flow_type == 6:
        # Converging Y. lam = inlet "a"'s fraction; phi = each inlet's
        # deflection.
        A = np.ones(3)
        Q = np.array([-1.0, lam, 1.0 - lam])
        th = np.array([0.0, pi + d, pi - d])
        slots = {"a": 0}
    else:
        raise ValueError(flow_type)
    return Q, A, th, slots


def model_K(spec: PanelSpec, x: float, config: dict | None = None) -> float:
    """The kernel's K for one panel at abscissa x (ignored on a Re panel)."""
    kw = {"lam": 0.5, "phi_deg": None, "area": 1.0}
    if spec.sweep == "lam":
        kw["lam"] = x
    elif spec.sweep == "phi":
        kw["phi_deg"] = x
    elif spec.sweep == "area":
        kw["area"] = x
    Q, A, th, slots = _junction(spec.flow_type, **kw)
    r = junction_loss_coefficient(Q / A, A, th, **(FAITHFUL if config is None else config))
    return float(r.K[slots[spec.branch]])


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def _csv(name: str) -> np.ndarray:
    a = np.loadtxt(_DATA / name, delimiter=",", comments="#", skiprows=2, ndmin=2)
    return a.reshape(-1, 2)


def calibration() -> dict:
    return yaml.safe_load((_DATA / "calibration.yaml").read_text())


def snap(spec: PanelSpec, x: float) -> float:
    if spec.sweep == "Re":
        return x
    grid = _DESIGN[spec.sweep]
    return float(grid[np.argmin(np.abs(grid - x))])


@dataclass(frozen=True)
class SamplePoint:
    panel: str
    x: float  # design abscissa the marker snaps to
    snap_px: float  # marker's distance from it, after the marker recalibration
    K_ref3d: float
    K_faithful: float
    K_production: float


@dataclass(frozen=True)
class PanelFidelity:
    panel: str
    n_red: int  # red points scored
    max_px: float
    rms_px: float
    max_K: float  # max_px in units of K


def _x_map(spec: PanelSpec, x_markers: np.ndarray) -> tuple[float, float]:
    """Refit x from the markers: each sits on a design value, which pins the
    axis better than the tick-label centroids (up to 2 px off on some
    panels). Identity on a Reynolds panel, where x does not enter."""
    if spec.sweep == "Re":
        return 1.0, 0.0
    design = np.array([snap(spec, x) for x in x_markers])
    if len(np.unique(design)) < 2:
        return 1.0, 0.0
    b, a = np.polyfit(x_markers, design, 1)
    return float(b), float(a)


def _spline(xk: np.ndarray, yk: np.ndarray, kind: str):
    from scipy.interpolate import CubicSpline, PchipInterpolator

    if len(xk) < 4 or kind == "linear":
        return lambda x: np.interp(x, xk, yk)
    if kind == "pchip":
        return PchipInterpolator(xk, yk)
    return CubicSpline(xk, yk, bc_type="not-a-knot")  # Matlab's spline()


def panel_fidelity(
    spec: PanelSpec, curve: str = "unified0d", config: dict | None = None, kind: str = "spline"
) -> PanelFidelity:
    """Distance of every visible red point from the model drawn Mynard's way:
    evaluated at his samples, splined between them."""
    cal = calibration()[spec.key]
    markers = _csv(f"mynard_{spec.key}_ref3d.csv")
    red = _csv(f"mynard_{spec.key}_{curve}.csv")
    sx, sk = cal["px_per_x"], cal["px_per_K"]
    if spec.sweep == "Re":
        d = np.abs(red[:, 1] - model_K(spec, 0.0, config)) * sk
    else:
        b, a = _x_map(spec, markers[:, 0])
        xr = a + b * red[:, 0]
        xk = np.unique([snap(spec, a + b * x) for x in markers[:, 0]])
        f = _spline(xk, np.array([model_K(spec, x, config) for x in xk]), kind)
        inside = (xr >= xk[0]) & (xr <= xk[-1])
        xr, kr = xr[inside], red[inside, 1]
        grid = np.linspace(xk[0], xk[-1], 2000)
        gk = f(grid)
        d = np.hypot((xr[:, None] - grid[None, :]) * sx, (kr[:, None] - gk[None, :]) * sk).min(axis=1)
    return PanelFidelity(spec.key, len(d), float(d.max()), float(np.sqrt((d**2).mean())), float(d.max() / sk))


def samples() -> list[SamplePoint]:
    cal_all = calibration()
    out = []
    for spec in PANELS:
        markers = _csv(f"mynard_{spec.key}_ref3d.csv")
        b, a = _x_map(spec, markers[:, 0])
        for x_m, k_ref in markers:
            xm = a + b * float(x_m)
            x = snap(spec, xm)
            out.append(
                SamplePoint(
                    panel=spec.key,
                    x=x,
                    snap_px=abs(x - xm) * cal_all[spec.key]["px_per_x"] if spec.sweep != "Re" else 0.0,
                    K_ref3d=float(k_ref),
                    K_faithful=model_K(spec, x),
                    K_production=model_K(spec, x, PRODUCTION),
                )
            )
    return out


def _main() -> None:
    pts = samples()
    print("FIDELITY -- port in Mynard's configuration vs his drawn Unified0D (model at his samples, splined)")
    print(f"{'panel':8s} {'n':>4s} {'max px':>7s} {'rms px':>7s} {'max K':>7s} {'pchip max':>9s} {'snap px':>8s}")
    for spec in PANELS:
        f = panel_fidelity(spec)
        g = panel_fidelity(spec, kind="pchip")
        snap_px = max(p.snap_px for p in pts if p.panel == spec.key)
        print(f"{spec.key:8s} {f.n_red:4d} {f.max_px:7.2f} {f.rms_px:7.2f} {f.max_K:7.3f} {g.max_px:9.2f} {snap_px:8.2f}")
    f4 = panel_fidelity(PANELS[0], "unified0d_eta0", dict(FAITHFUL, eta_scale=0.0))
    print(f"fig04 eta_j=0 dashed: n {f4.n_red} max {f4.max_px:.2f} px rms {f4.rms_px:.2f}")
    print()
    print("REF3D -- Mynard's own CFD (laminar, Re 350-2400); production departs on purpose; NOT tuned on")
    print(f"{'':26s} {'faithful':>23s} {'production':>23s}")
    print(f"{'group':22s} {'n':>3s} {'MAE':>7s} {'median':>7s} {'bias':>7s} {'MAE':>7s} {'median':>7s} {'bias':>7s}")
    for label, sel in ref3d_groups(pts):
        ef = np.array([p.K_faithful - p.K_ref3d for p in sel])
        ep = np.array([p.K_production - p.K_ref3d for p in sel])
        print(
            f"{label:22s} {len(sel):3d}"
            f" {np.abs(ef).mean():7.3f} {np.median(np.abs(ef)):7.3f} {ef.mean():+7.3f}"
            f" {np.abs(ep).mean():7.3f} {np.median(np.abs(ep)):7.3f} {ep.mean():+7.3f}"
        )


_TYPE_NAME = {
    1: "T diverging",
    2: "T from side",
    3: "Y diverging",
    34: "Y diverging, Fig 4",
    4: "T converging",
    5: "T into side",
    6: "Y converging",
}


def ref3d_groups(pts: list[SamplePoint]) -> list[tuple[str, list[SamplePoint]]]:
    """Ref3D samples grouped by flow type and branch, then the psi != 1
    sweeps on their own (the only place the joining alpha acts), then all."""
    spec_of = {s.key: s for s in PANELS}
    groups: dict[str, list[SamplePoint]] = {}
    for p in pts:
        s = spec_of[p.panel]
        name = _TYPE_NAME[s.flow_type] + ("" if s.branch == "a" else f" {s.branch}")
        groups.setdefault(name, []).append(p)
    out = list(groups.items())
    out.append(("joining, psi != 1", [p for p in pts if p.panel in ("fig09d", "fig09h") and p.x != 1.0]))
    out.append(("all", pts))
    return out

if __name__ == "__main__":
    _main()
