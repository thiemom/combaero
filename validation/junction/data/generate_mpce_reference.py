"""Emit golden values for the C++ whole-element MPCE (f, J) port.

Residual values come from the shipping Python element, so the port's gate is
equivalence with what runs today.

The DERIVATIVES are central differences of that residual -- NOT the Python's
own analytic Jacobian. That is deliberate. The Python assembles its Jacobian
from an explicit dKQ/dmdot block plus a hand-derived dR/dP column for the
common port alone, and misses the other two ports' static-pressure columns
(K depends on every port's velocity, and every velocity on its own density).
Measured at a converged state, those columns are off by 4.4e-3 and 2.3e-3.
The C++ seeds the whole element, so it produces them; comparing against the
Python's analytic Jacobian would therefore mark the C++ wrong where it is
right. Finite differences of the residual are the truth for both.

Run:

    uv run python validation/junction/data/generate_mpce_reference.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import combaero as cb  # noqa: E402
from combaero.network.components import NetworkMixtureState  # noqa: E402
from combaero.network.mpce_element import MultiPortChamberElement  # noqa: E402

_OUT = Path(__file__).with_name("mpce_reference_data.h")
_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_FD_REL = 1e-6

# Seed layout, mirrored from include/mpce_junction.h.
_N = 3
_SEEDS = 3 * _N + 1


def _element(theta_deg: list[float], areas: list[float], joining: bool,
             alpha: float, eta: float) -> MultiPortChamberElement:
    ports = ["p0", "p1", "p2"]
    if joining:
        inlets, outlets = ports[:2], ports[2:]
        in_ang, out_ang = theta_deg[:2], theta_deg[2:]
        direction = "merge"
    else:
        inlets, outlets = ports[:1], ports[1:]
        in_ang, out_ang = theta_deg[:1], theta_deg[1:]
        direction = "branch"
    element = MultiPortChamberElement(
        id="jct", inlet_nodes=inlets, outlet_nodes=outlets,
        inlet_angles_deg=in_ang, outlet_angles_deg=out_ang,
        port_areas=areas, flow_direction=direction, strict=False,
        joining_etransfer_alpha=alpha, eta_scale=eta,
    )
    element._port_element_ids = ["e0", "e1", "e2"]
    return element


def _states(P: list[float], Pt: list[float], T: float) -> list[NetworkMixtureState]:
    # Tt and m_dot play no part in the junction residual; the element reads
    # Pt, and density() reads T, P and Y.
    return [
        NetworkMixtureState(P=p, Pt=pt, T=T, Tt=T, m_dot=0.0, Y=list(_Y))
        for p, pt in zip(P, Pt, strict=True)
    ]


def _residual(case: dict, x: np.ndarray) -> np.ndarray:
    """Residual as a function of the 10 seeded unknowns."""
    P = list(x[0:3])
    Pt = list(x[3:6])
    outer = list(x[6:9])
    pt_jct = float(x[9])
    element = case["element"]
    states = _states(P, Pt, case["T"])
    mdots = [s * m for s, m in zip(element._port_signs, outer, strict=True)]
    res, _ = element.residuals(states, pt_jct, mdots)
    return np.array(res, dtype=float)


def _cases() -> list[dict]:
    d45 = 45.0
    d90 = 90.0
    d156 = 156.26
    out = []

    def add(cid, theta, areas, joining, P, Pt, outer, pt_jct, T=320.0,
            alpha=0.0, eta=0.0):
        out.append(dict(
            id=cid, element=_element(theta, areas, joining, alpha, eta),
            theta=theta, areas=areas, joining=joining, alpha=alpha, eta=eta,
            P=P, Pt=Pt, outer=outer, pt_jct=pt_jct, T=T,
        ))

    # Separating: port 0 in, ports 1 and 2 out. port_signs = [-1, +1, +1],
    # so a positive outer mdot on port 0 flows INTO the junction.
    add("sep_equal_area_90", [0.0, 0.0, d90], [0.01, 0.01, 0.01], False,
        P=[2.00e5, 1.97e5, 1.95e5], Pt=[2.06e5, 2.02e5, 2.00e5],
        outer=[0.90, 0.55, 0.35], pt_jct=2.05e5)
    add("sep_equal_area_45", [0.0, 0.0, d45], [0.01, 0.01, 0.01], False,
        P=[2.00e5, 1.97e5, 1.96e5], Pt=[2.06e5, 2.02e5, 2.01e5],
        outer=[0.90, 0.65, 0.25], pt_jct=2.05e5)
    add("sep_area_ratio_4", [0.0, 0.0, d45], [0.01, 0.01, 0.0025], False,
        P=[2.00e5, 1.97e5, 1.80e5], Pt=[2.06e5, 2.02e5, 2.00e5],
        outer=[0.90, 0.72, 0.18], pt_jct=2.05e5)
    add("sep_small_branch", [0.0, 0.0, d90], [0.01, 0.01, 0.01], False,
        P=[2.00e5, 1.98e5, 1.99e5], Pt=[2.06e5, 2.03e5, 2.02e5],
        outer=[0.90, 0.86, 0.04], pt_jct=2.05e5)
    add("sep_eta_on", [0.0, 0.0, d90], [0.01, 0.01, 0.01], False,
        P=[2.00e5, 1.97e5, 1.95e5], Pt=[2.06e5, 2.02e5, 2.00e5],
        outer=[0.90, 0.55, 0.35], pt_jct=2.05e5, eta=1.0)
    add("sep_high_pressure", [0.0, 0.0, d45], [0.02, 0.02, 0.005], False,
        P=[9.0e5, 8.9e5, 8.5e5], Pt=[9.2e5, 9.05e5, 9.0e5],
        outer=[6.0, 4.5, 1.5], pt_jct=9.15e5, T=520.0)

    # Joining: ports 0 and 1 in, port 2 out. port_signs = [-1, -1, +1].
    add("join_equal_area_90", [0.0, d90, 0.0], [0.01, 0.01, 0.01], True,
        P=[2.02e5, 2.01e5, 1.98e5], Pt=[2.08e5, 2.06e5, 2.02e5],
        outer=[0.55, 0.35, 0.90], pt_jct=2.03e5)
    add("join_obtuse_lateral", [0.0, d156, 0.0], [0.01, 0.01, 0.01], True,
        P=[2.02e5, 2.01e5, 1.98e5], Pt=[2.08e5, 2.06e5, 2.02e5],
        outer=[0.55, 0.35, 0.90], pt_jct=2.03e5)
    add("join_area_ratio_2p5", [0.0, d45, 0.0], [0.01, 0.004, 0.01], True,
        P=[2.02e5, 1.95e5, 1.98e5], Pt=[2.08e5, 2.10e5, 2.02e5],
        outer=[0.55, 0.35, 0.90], pt_jct=2.03e5)
    add("join_alpha_on", [0.0, d45, 0.0], [0.01, 0.004, 0.01], True,
        P=[2.02e5, 1.95e5, 1.98e5], Pt=[2.08e5, 2.10e5, 2.02e5],
        outer=[0.55, 0.35, 0.90], pt_jct=2.03e5, alpha=0.2)
    add("join_low_pressure", [0.0, d90, 0.0], [0.005, 0.005, 0.005], True,
        P=[0.6e5, 0.59e5, 0.55e5], Pt=[0.64e5, 0.63e5, 0.58e5],
        outer=[0.10, 0.06, 0.16], pt_jct=0.60e5, T=280.0)
    return out


def _fmt(values) -> str:
    return ", ".join(f"{float(v):.17g}" for v in np.ravel(values))


def main() -> None:
    rows = []
    for case in _cases():
        x = np.array(case["P"] + case["Pt"] + case["outer"] + [case["pt_jct"]])
        res = _residual(case, x)

        jac_fd = np.zeros((len(res), _SEEDS))
        for j in range(_SEEDS):
            h = max(abs(x[j]) * _FD_REL, 1e-9)
            xp = x.copy(); xp[j] += h
            xm = x.copy(); xm[j] -= h
            jac_fd[:, j] = (_residual(case, xp) - _residual(case, xm)) / (2.0 * h)

        # rho and drho/dP at the evaluation point, from the same states the
        # element sees, so the C++ starts from identical thermodynamics.
        states = _states(case["P"], case["Pt"], case["T"])
        rho = np.array([float(s.density()) for s in states])
        drho = np.zeros(_N)
        for i, s in enumerate(states):
            h = max(abs(case["P"][i]) * _FD_REL, 1.0)
            hi = NetworkMixtureState(P=case["P"][i] + h, Pt=case["Pt"][i], T=case["T"],
                                     Tt=case["T"], m_dot=0.0, Y=list(_Y))
            lo = NetworkMixtureState(P=case["P"][i] - h, Pt=case["Pt"][i], T=case["T"],
                                     Tt=case["T"], m_dot=0.0, Y=list(_Y))
            drho[i] = (float(hi.density()) - float(lo.density())) / (2.0 * h)

        element = case["element"]
        theta_rad = [math.radians(t) for t in element.port_angles_deg]
        rows.append(
            "    {"
            f'"{case["id"]}",\n'
            f"     {{{_fmt(case['P'])}}}, {{{_fmt(case['Pt'])}}},\n"
            f"     {{{_fmt(rho)}}}, {{{_fmt(drho)}}},\n"
            f"     {{{_fmt(case['outer'])}}}, {case['pt_jct']:.17g},\n"
            f"     {{{_fmt(case['areas'])}}}, {{{_fmt(theta_rad)}}},\n"
            f"     {{{_fmt(element._port_signs)}}},\n"
            f"     {case['alpha']:.17g}, {case['eta']:.17g},\n"
            f"     {{{_fmt(res)}}},\n"
            f"     {{{_fmt(jac_fd)}}}}},"
        )

    body = "\n".join(rows)
    _OUT.write_text(
        "// AUTO-GENERATED by validation/junction/data/generate_mpce_reference.py\n"
        "// Do not edit by hand. Golden values for the C++ whole-element (f, J)\n"
        "// port of MultiPortChamberElement, produced by the Python that ships.\n"
        "//\n"
        "// The Jacobian rows are CENTRAL DIFFERENCES of the Python RESIDUAL, not\n"
        "// the Python's own analytic Jacobian: the Python's misses the two\n"
        "// non-common ports' static-pressure columns, which whole-element\n"
        "// seeding supplies. See the generator's docstring.\n"
        "#pragma once\n\n"
        "#include <array>\n\n"
        "namespace combaero::validation::junction {\n\n"
        "struct MpceCase {\n"
        "  const char *id;\n"
        "  std::array<double, 3> p_static;\n"
        "  std::array<double, 3> p_total;\n"
        "  std::array<double, 3> rho;\n"
        "  std::array<double, 3> drho_dp;\n"
        "  std::array<double, 3> outer_mdot;\n"
        "  double pt_jct;\n"
        "  std::array<double, 3> area;\n"
        "  std::array<double, 3> theta_rad;\n"
        "  std::array<double, 3> port_sign;\n"
        "  double joining_etransfer_alpha;\n"
        "  double eta_scale;\n"
        "  std::array<double, 4> residual;\n"
        "  std::array<double, 40> jacobian_fd; // [row][seed], row-major\n"
        "};\n\n"
        f"inline constexpr std::array<MpceCase, {len(rows)}> kMpceCases{{{{\n"
        f"{body}\n}}}};\n\n"
        "} // namespace combaero::validation::junction\n"
    )
    print(f"wrote {_OUT} ({len(rows)} cases)")


if __name__ == "__main__":
    main()
