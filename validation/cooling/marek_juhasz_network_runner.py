"""NASA TN D-7156 (Marek and Juhasz 1973) as a combaero network (#764).

A 304 stainless plate in a combustor's exhaust, cooled by a 2-D film slot on
its hot face, a parallel-flow convection passage under it, or both. Ten
thermocouples brazed into grooves on the HOT face measure the wall
temperature at x = 1.02 ... 11.99 cm downstream of the slot (1.22 cm pitch).

The network, per run (Table II, SI, transcribed as printed):

    hot_in --[lead]--> g0 --[h1]--> g1 -- ... --[h10]--> g10 --[tail]--> out
                        ^
    film_in --[SlotFilmElement]
    conv_in --[c1]--> k1 -- ... --[c10]--> k10 --[conv_exit]--> conv_out

with one ThermalWall per thermocouple, centred on it (1.22 cm long, 10.2 cm
wide), coupling h_i (film side, ``film=`` the slot) to c_i through 0.159 cm
of stainless. The walls span x = 0.41-12.60 cm behind the slot. The
backside's other face (the stagnant gap over the slave-air-cooled underside)
is adiabatic. In a run without convection air the passage is stagnant too:
it carries a token flow whose face is switched off (``STAGNANT_FLOW``), so
the wall sits at the film's adiabatic wall temperature. The predicted wall
temperature is each wall's hot-face temperature, where the thermocouples
sit. Every solve starts from the outlet-referenced warm start: the default
cold start fails on these networks and then retries the same way, at
about four times the cost.

The slot discharges into g0, at the head of the walled duct, as a user would
build it. Segments 2-10 therefore read gas that already holds the film's
flow (#694 un-dilutes only the first), and M is formed on that mixed flow.
``Options(slot_into="tail")`` is a diagnostic in which every segment sees
undiluted gas; it moves the scores by under 1 K.

WHAT A MISS MEANS. Nothing here is fitted to D-7156. The film is D-7958's
kernel (Marek and Tacina 1975) at the SlotFilmElement default Tu 0.14
(C_m 0.035), not the C_m 0.03 the authors fitted to these data; the hot-side
h without a film is combaero's own duct h (Gnielinski), not the authors'
tuned 0.0365 Re^0.8; the film-side h is Seban (1960). Same lab, different
study: every number this module reports is CROSS-SOURCE ACCURACY
(docs/VALIDATION_POLICY.md), never fidelity, and nothing is tuned to it.

LABELLED ASSUMPTIONS (the user's sign-off, 2026-11-03):

- Hot gas is vitiated air. The preheater burns ASTM A-1 (Jet A-1) to lift
  ambient air (290 K, dry) to 589 K; the test combustor fires the same fuel
  for the 1033 K condition (assumed, labelled). combaero carries no kerosene
  species; the heaviest alkane it has, n-heptane (NC7H16), stands in, burnt
  gaseous at 290 K to CO2 and H2O, with f/a solved for each temperature rise
  (290 K -> 589 K, then 589 K -> the run's T_h; condition-1 runs fire only
  the preheater, straight to T_h). Coolant flows stay plain dry air.
- The hot-gas flow area 0.0370 m^2 (38.1 x 10.6 cm duct less the 3.3 x
  10.2 cm test section) returns Table I's velocities within 1 %. The
  hydraulic diameter is 4 A / P over that section, both perimeters
  included: 11.9 cm.
- T_h is the rake average (pattern factor 0.33), T_s as printed (the
  authors inferred it, their film thermocouple having failed).
- A slot lip of 1 mm: it sets only the slot's pressure drop, irrelevant with
  the film flow prescribed.
- The wall's k(T): combaero has no type-304 entry; ss316 (its nearest
  austenitic stainless) stands in. The printed SI k is a factor ~17 too low
  (the US value is right) and is not used. The thin wall makes the choice
  immaterial: a constant k of 13.8 or 22.5 W/m K (304 at 311 K and 811 K)
  moves the combined stage by 0.3 K.
- Thermocouples 1-3 are scored and reported apart from 4-10: combaero has no
  axial conduction and no slot-lip heat sink, and the authors saw conduction
  effects out to 3.8 cm.
"""

from __future__ import annotations

import csv
import math
import warnings
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SLOT_FILM_DEFAULT_TU
from validation.cooling.schema import DATA_ROOT

DATA = DATA_ROOT / "marek_juhasz1973"

# --- geometry (pp. 9-11, Fig. 4) ------------------------------------------
P_ATM = 101325.0
SLOT_HEIGHT = 0.0038  # s, 0.38 cm
SPAN = 0.102  # test surface width, 10.2 cm
WALL_THICKNESS = 0.00159  # 0.159 cm type-304
PASSAGE_HEIGHT = 0.003  # convection passage, 0.3 cm
X_LEAD = 0.060  # x_o, leading edge to slot
TC_X = tuple(0.0102 + 0.0122 * i for i in range(10))  # thermocouples [m]
TC_PITCH = 0.0122
LIP_THICKNESS = 0.001  # ASSUMED (irrelevant with the film flow prescribed)
# Hot gas: 38.1 x 10.6 cm duct less the 3.3 x 10.2 cm test section.
HOT_AREA = 0.381 * 0.106 - 0.033 * 0.102
HOT_PERIMETER = 2.0 * (0.381 + 0.106) + 2.0 * (0.033 + 0.102)
HOT_DH = 4.0 * HOT_AREA / HOT_PERIMETER
PASSAGE_AREA = PASSAGE_HEIGHT * SPAN
PASSAGE_DH = 2.0 * PASSAGE_HEIGHT * SPAN / (PASSAGE_HEIGHT + SPAN)
SEGMENT_AREA = TC_PITCH * SPAN
WALL_MATERIAL = "ss316"  # stands in for type 304 (see the module docstring)
#: Runs without convection air: the passage is stagnant and the plate's
#: underside adiabatic. A ThermalWall needs a stream on each side, so the
#: passage carries a token flow whose face is switched off (its h scaled by
#: 1e-6, under the wall's 0.01 W/m^2K floor): the wall then passes ~mW and
#: sits at its filmed side's adiabatic wall temperature. This models the
#: adiabatic face; it is not a rig-matching tuner and moves no cooled run.
STAGNANT_FLOW = 1e-3
ADIABATIC_FACE_NU_MULTIPLIER = 1e-6
T_BAND = 11.0  # K, the overall accuracy of the test-surface temperatures (p. 13)

# --- the vitiated hot gas -------------------------------------------------
T_AMBIENT = 290.0
T_PREHEAT = 589.0
SURROGATE = "NC7H16"
#: Condition 1 (589 K): the preheater alone; above this, the combustor fires.
PREHEAT_ONLY_BELOW = 700.0

# --- run classes (Table II; p. 13) ----------------------------------------
UNCOOLED = (1, 14)
FILM_ONLY = (*range(2, 14), *range(20, 26), 33, *range(38, 44), 55, 56, 57)
CONVECTION_ONLY = (*range(15, 20), 34, 35, 36, *range(44, 48), 58, 59, 60)
COMBINED = (*range(26, 33), 37, *range(48, 55), 61, 62)
#: Combined runs inside or above Seban's box (m >= 1); the rest are m < 1,
#: where the film-side h hands back to the duct's.
COMBINED_M_GE_1 = (29, 30, 31, 32, 37, 54, 62)
TC_NEAR = (0, 1, 2)  # thermocouples 1-3 (x 1.02-3.45 cm)
TC_FAR = tuple(range(3, 10))  # thermocouples 4-10

STAGES = {
    "uncooled": UNCOOLED,
    "convection": CONVECTION_ONLY,
    "film": FILM_ONLY,
    "combined": COMBINED,
}


@dataclass(frozen=True)
class Run:
    run: int
    m_hot: float
    T_hot: float
    m_film: float
    T_film: float
    m_conv: float
    T_conv_in: float
    T_conv_exit: float
    T_wall: tuple[float, ...]

    @property
    def stage(self) -> str:
        return next(k for k, v in STAGES.items() if self.run in v)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


@lru_cache(maxsize=1)
def runs() -> dict[int, Run]:
    out = {}
    for r in _rows(DATA / "table2a_experimental_data_SI.csv"):
        tw = tuple(float(v) for k, v in r.items() if k.startswith("Tw_"))
        out[int(r["run"])] = Run(
            run=int(r["run"]),
            m_hot=float(r["m_hot_kg_s"]),
            T_hot=float(r["T_hot_K"]),
            m_film=float(r["m_film_kg_s"]),
            T_film=float(r["T_film_in_K"]),
            m_conv=float(r["m_conv_kg_s"]),
            T_conv_in=float(r["T_conv_in_K"]),
            T_conv_exit=float(r["T_conv_exit_K"]),
            T_wall=tw,
        )
    return out


# --- composition ----------------------------------------------------------


def _stream(T: float, X: list[float], m: float) -> cb.Stream:
    s = cb.Stream()
    s.set_T(T)
    s.set_P(P_ATM)
    s.set_X(list(X))
    s.set_mdot(m)
    return s


def _fuel_X() -> list[float]:
    X = [0.0] * len(cb.species.names)
    X[cb.species.names.index(SURROGATE)] = 1.0
    return X


def _burn(T_in: float, X_in: list[float], T_out: float) -> tuple[float, list[float]]:
    """Fuel at T_AMBIENT into X_in at T_in, complete combustion to T_out:
    (fuel/oxidiser mass ratio, product mole fractions)."""
    fuel = cb.set_fuel_stream_for_Tad(
        T_out, _stream(T_AMBIENT, _fuel_X(), 0.0), _stream(T_in, X_in, 1.0), tol=1e-3
    )
    fa = float(fuel.mdot)
    Y = (np.asarray(cb.mole_to_mass(list(X_in))) + fa * np.asarray(cb.mole_to_mass(_fuel_X()))) / (
        1.0 + fa
    )
    X_r = cb.mass_to_mole(list(Y))
    return fa, list(cb.complete_combustion_to_CO2_H2O(list(X_r)))


@lru_cache(maxsize=None)
def vitiated(T_hot: float) -> tuple[tuple[float, ...], float]:
    """(mass fractions, overall fuel/air mass ratio) of the hot gas at T_hot."""
    X_air = list(cb.species.dry_air())
    if T_hot < PREHEAT_ONLY_BELOW:
        fa, X = _burn(T_AMBIENT, X_air, T_hot)
    else:
        fa1, X1 = _burn(T_AMBIENT, X_air, T_PREHEAT)
        fa2, X = _burn(T_PREHEAT, X1, T_hot)
        fa = fa1 + fa2 * (1.0 + fa1)
    return tuple(cb.mole_to_mass(X)), fa


# --- the network ----------------------------------------------------------


@dataclass(frozen=True)
class Options:
    film_htc: str = "seban1960"
    turbulence_intensity: float = SLOT_FILM_DEFAULT_TU
    wall_k: float | None = None  # a constant k instead of the material's k(T)
    hot_gas: str = "vitiated"  # or "air" (diagnostic)
    #: Where the slot's flow joins the hot stream: "head" (g0, the network a
    #: user builds: segments 2-10 then read the mixed gas, and M the mixed
    #: flow) or "tail" (g10, a diagnostic: every segment sees the undiluted
    #: hot gas and M is the report's rho_s u_s / rho_h u_h).
    slot_into: str = "head"


def build(run: Run, opt: Options = Options()) -> FlowNetwork:
    Y_air = cb.mole_to_mass(cb.species.dry_air())
    Y_hot = list(vitiated(run.T_hot)[0]) if opt.hot_gas == "vitiated" else Y_air
    D_hot = math.sqrt(4.0 * HOT_AREA / math.pi)
    D_pass = math.sqrt(4.0 * PASSAGE_AREA / math.pi)
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("hot_in", m_dot=run.m_hot, Tt=run.T_hot, Y=Y_hot))
    g.add_node(PressureBoundary("out", Pt=P_ATM, Tt=run.T_hot, Y=Y_hot))
    g.add_node(MassFlowBoundary("conv_in", m_dot=max(run.m_conv, STAGNANT_FLOW), Tt=run.T_conv_in, Y=Y_air))
    g.add_node(PressureBoundary("conv_out", Pt=P_ATM, Tt=run.T_conv_in, Y=Y_air))
    for i in range(11):
        g.add_node(PlenumNode(f"g{i}"))
    for i in range(1, 11):
        g.add_node(PlenumNode(f"k{i}"))

    def hot(eid: str, a: str, b: str, length: float, area: float = 0.0) -> ChannelElement:
        return ChannelElement(
            eid,
            a,
            b,
            length=length,
            diameter=D_hot,
            Dh=HOT_DH,
            roughness=0.0,
            surface=ConvectiveSurface(area=area) if area else None,
        )

    g.add_element(hot("lead", "hot_in", "g0", X_LEAD))
    if run.m_film > 0.0:
        g.add_node(MassFlowBoundary("film_in", m_dot=run.m_film, Tt=run.T_film, Y=Y_air))
        g.add_element(
            SlotFilmElement(
                "slot",
                "film_in",
                "g0" if opt.slot_into == "head" else "g10",
                slot_height=SLOT_HEIGHT,
                span=SPAN,
                lip_thickness=LIP_THICKNESS,
                turbulence_intensity=opt.turbulence_intensity,
                film_htc=opt.film_htc,
            )
        )
    for i in range(1, 11):
        g.add_element(hot(f"h{i}", f"g{i - 1}", f"g{i}", TC_PITCH, SEGMENT_AREA))
        g.add_element(
            ChannelElement(
                f"c{i}",
                "conv_in" if i == 1 else f"k{i - 1}",
                f"k{i}",
                length=TC_PITCH,
                diameter=D_pass,
                Dh=PASSAGE_DH,
                roughness=0.0,
                surface=ConvectiveSurface(
                    area=SEGMENT_AREA,
                    Nu_multiplier=1.0 if run.m_conv > 0.0 else ADIABATIC_FACE_NU_MULTIPLIER,
                ),
            )
        )
        x0 = TC_X[i - 1] - 0.5 * TC_PITCH
        layer = (
            WallLayer(WALL_THICKNESS, opt.wall_k)
            if opt.wall_k is not None
            else WallLayer(WALL_THICKNESS, 17.0, material=WALL_MATERIAL)
        )
        g.add_wall(
            ThermalWall(
                id=f"w{i}",
                element_a=f"h{i}",
                element_b=f"c{i}",
                layers=[layer],
                contact_area=SEGMENT_AREA,
                film="slot" if run.m_film > 0.0 else None,
                film_start=x0,
                film_end=x0 + TC_PITCH,
            )
        )
    g.add_element(hot("tail", "g10", "out", 0.01))
    g.add_element(
        ChannelElement(
            "conv_exit", "k10", "conv_out", length=0.01, diameter=D_pass, Dh=PASSAGE_DH, roughness=0.0
        )
    )
    return g


@dataclass
class Result:
    run: Run
    T_wall: tuple[float, ...]
    T_conv_exit: float | None
    extras: dict = field(default_factory=dict)

    @property
    def errors(self) -> np.ndarray:
        return np.asarray(self.T_wall) - np.asarray(self.run.T_wall)


def solve(run: Run, opt: Options = Options()) -> Result:
    """Solve one run; the predicted wall temperature is each wall's hot-face
    temperature (the thermocouples are brazed into the hot face)."""
    net = build(run, opt)
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(init_strategy="outletref_warmstart")
    if not r["__success__"]:
        raise RuntimeError(f"run {run.run}: {r.get('__message__')}")
    keys = ("eta_film", "film_blowing_ratio", "h_a", "h_b", "film_h_weight_M", "T_aw_a", "Q")
    extras = {k: [r.get(f"w{i}.{k}") for i in range(1, 11)] for k in keys}
    T_wall = tuple(float(r[f"w{i}.T_interface"][0]) for i in range(1, 11))
    T_exit = float(r["k10.Tt"]) if run.m_conv > 0.0 else None
    extras["material_in_range"] = [r.get(f"w{i}.material_in_range") for i in range(1, 11)]
    return Result(run, T_wall, T_exit, extras)


@dataclass(frozen=True)
class Score:
    n: int
    bias: float
    rms: float
    worst: float
    within_band: float  # share of points inside +-T_BAND

    def __str__(self) -> str:
        return (
            f"n {self.n:4d}  bias {self.bias:+7.1f} K  RMS {self.rms:6.1f} K  "
            f"max {self.worst:+7.1f} K  in band {100 * self.within_band:5.1f} %"
        )


def score(errors: np.ndarray) -> Score:
    e = np.asarray(errors, dtype=float).ravel()
    worst = float(e[np.argmax(np.abs(e))])
    return Score(
        int(e.size),
        float(e.mean()),
        float(np.sqrt((e**2).mean())),
        worst,
        float((np.abs(e) <= T_BAND).mean()),
    )


def _solve_star(args: tuple[int, Options]) -> Result:
    return solve(runs()[args[0]], args[1])


def solve_all(run_ids: tuple[int, ...], opt: Options = Options(), jobs: int = 1) -> list[Result]:
    if jobs <= 1:
        return [solve(runs()[i], opt) for i in run_ids]
    from concurrent.futures import ProcessPoolExecutor

    try:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            return list(ex.map(_solve_star, [(i, opt) for i in run_ids]))
    except PermissionError:  # no process pool here (a sandbox): serially
        return solve_all(run_ids, opt, jobs=1)


def stage_scores(results: list[Result]) -> dict[str, Score]:
    """Pooled, TC 1-3, TC 4-10."""
    E = np.array([r.errors for r in results])
    return {
        "all": score(E),
        "TC1-3": score(E[:, list(TC_NEAR)]),
        "TC4-10": score(E[:, list(TC_FAR)]),
    }


#: Diagnostics, never adopted: Tu 0.12 is C_m 0.03, the authors' fit to
#: these data; 0.07 and 0.23 are D-7958's lowest and next blockage steps.
TU_SWEEP = (0.07, 0.12, 0.14, 0.23)
#: Constant-k wall bounds, type 304 at 311 K and 811 K (the report's US
#: values, 8 and 13 Btu/hr-ft-F).
WALL_K_SWEEP = (13.8, 22.5)


def _lines(res: list[Result], indent: str = "  ") -> list[str]:
    return [f"{indent}{k:8s} {sc}" for k, sc in stage_scores(res).items()]


def _split(res: list[Result]) -> tuple[list[Result], list[Result]]:
    hi = [r for r in res if r.run.run in COMBINED_M_GE_1]
    lo = [r for r in res if r.run.run not in COMBINED_M_GE_1]
    return hi, lo


def report(jobs: int = 4, per_run: bool = False) -> str:
    lines = [
        "NASA TN D-7156 (Marek & Juhasz 1973) network case -- CROSS-SOURCE ACCURACY",
        "(same lab, different study; nothing tuned to these data). Wall T predicted - measured,",
        f"band +-{T_BAND:.0f} K (p. 13). TC 1-3 sit where conduction to the slot lip acts",
        "(the authors saw it to 3.8 cm), which combaero does not model.",
    ]
    base: dict[str, list[Result]] = {}
    for stage, ids in STAGES.items():
        res = solve_all(ids, jobs=jobs)
        base[stage] = res
        lines.append(f"\n[{stage}] {len(ids)} runs")
        lines += _lines(res)
        if stage in ("convection", "film"):
            for label, sel in (
                ("589 K (condition 1)", [r for r in res if r.run.T_hot < PREHEAT_ONLY_BELOW]),
                ("1033 K (conditions 2-4)", [r for r in res if r.run.T_hot >= PREHEAT_ONLY_BELOW]),
            ):
                lines.append(f"  {label}: {len(sel)} runs  {stage_scores(sel)['all']}")
        if stage == "combined":
            hi, lo = _split(res)
            lines.append(f"  m >= 1 (Seban's box or above): runs {[r.run.run for r in hi]}")
            lines += _lines(hi, "    ")
            lines.append(f"  m < 1 (film-side h hands back to the duct's): {len(lo)} runs")
            lines += _lines(lo, "    ")
        exits = [(r.T_conv_exit - r.run.T_conv_exit) for r in res if r.T_conv_exit is not None]
        if exits:
            lines.append(f"  convection exit T (energy balance)  {score(np.array(exits))}")
        if per_run:
            for r in res:
                M = r.extras["film_blowing_ratio"][0]
                lines.append(
                    f"    run {r.run.run:2d}  M {M if M is None else round(M, 2)!s:5s} "
                    f"err {np.round(r.errors, 0).astype(int).tolist()}"
                )
    lines.append("\nMETHOD COMPARISON (combined stage): the film-side h")
    for model in ("seban1960", "none", "ballal_lefebvre1973"):
        res = base["combined"] if model == "seban1960" else solve_all(
            COMBINED, Options(film_htc=model), jobs=jobs
        )
        hi, lo = _split(res)
        lines.append(f"  {model}")
        lines.append(f"    all      {stage_scores(res)['all']}")
        lines.append(f"    m >= 1   {stage_scores(hi)['all']}")
        lines.append(f"    m < 1    {stage_scores(lo)['all']}")
    lines.append("\nDIAGNOSTIC (not adopted): film-only stage against Tu (C_m = 0.25 Tu)")
    for tu in TU_SWEEP:
        res = base["film"] if tu == SLOT_FILM_DEFAULT_TU else solve_all(
            FILM_ONLY, Options(turbulence_intensity=tu), jobs=jobs
        )
        lines.append(f"  Tu {tu:4.2f}  {stage_scores(res)['all']}")
    lines.append("\nDIAGNOSTIC: combined stage, constant wall k [W/m K] instead of ss316 k(T)")
    for k in WALL_K_SWEEP:
        lines.append(f"  k {k:4.1f}  {stage_scores(solve_all(COMBINED, Options(wall_k=k), jobs=jobs))['all']}")
    lines.append("\nDIAGNOSTIC: the slot's flow joins the hot stream after the wall (undiluted)")
    for stage in ("film", "combined"):
        res = solve_all(STAGES[stage], Options(slot_into="tail"), jobs=jobs)
        lines.append(f"  {stage:10s} {stage_scores(res)['all']}")
    lines.append("\nDIAGNOSTIC: hot gas as plain air instead of vitiated")
    for stage in ("convection", "combined"):
        res = solve_all(STAGES[stage], Options(hot_gas="air"), jobs=jobs)
        lines.append(f"  {stage:10s} {stage_scores(res)['all']}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys

    print(report(per_run="--per-run" in sys.argv))
