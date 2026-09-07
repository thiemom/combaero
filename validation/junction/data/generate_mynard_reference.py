"""Emit golden values for the C++ Mynard closure port.

The C++ port's gate is equivalence with the Python that ships, so the
reference numbers are produced BY that Python rather than typed in. Run:

    uv run python validation/junction/data/generate_mynard_reference.py

Cases are chosen to reach every branch the port takes, not to look tidy:
both flow directions, both mask orientations, an equal-area case and an
extreme area ratio, angles either side of the second-quadrant flip, a
collinear lateral (the dividing-streamline term) and one that just misses it,
the energy-transfer term on and off, and the joining asymmetry term on and
off. Derivatives are central differences in each velocity, which is what the
C++ dual partials are checked against.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from combaero.network._mynard2010 import junction_loss_coefficient  # noqa: E402

_OUT = Path(__file__).with_name("mynard_reference_data.h")
_FD_STEP = 1e-7


def _cases() -> list[dict]:
    d45 = math.radians(45.0)
    d90 = math.radians(90.0)
    d156 = math.radians(156.26)
    pi = math.pi
    return [
        # --- dividing: one supplier, two collectors -------------------------
        dict(id="div_equal_area_90", U=[10.0, -6.0, -4.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, d90], alpha=0.0, eta=0.0),
        dict(id="div_equal_area_45", U=[10.0, -7.0, -3.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, d45], alpha=0.0, eta=0.0),
        dict(id="div_area_ratio_4", U=[10.0, -6.0, -16.0], A=[0.01, 0.01, 0.0025],
             theta=[0.0, pi, d45], alpha=0.0, eta=0.0),
        dict(id="div_small_split", U=[10.0, -9.8, -0.2], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, d90], alpha=0.0, eta=0.0),
        dict(id="div_eta_on", U=[10.0, -6.0, -4.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, d90], alpha=0.0, eta=1.0),
        # Collinear continuing collector: the dividing-streamline term fires.
        dict(id="div_collinear", U=[10.0, -6.0, -4.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, pi], alpha=0.0, eta=0.0),
        # And one that just misses the collinearity tolerance.
        dict(id="div_near_collinear", U=[10.0, -6.0, -4.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, pi - 1.0e-4], alpha=0.0, eta=0.0),
        # --- joining: two suppliers, one collector --------------------------
        dict(id="join_equal_area_90", U=[6.0, 4.0, -10.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, d90, pi], alpha=0.0, eta=0.0),
        dict(id="join_obtuse_lateral", U=[6.0, 4.0, -10.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, d156, pi], alpha=0.0, eta=0.0),
        dict(id="join_area_ratio_2p5", U=[6.0, 10.0, -10.0], A=[0.01, 0.004, 0.01],
             theta=[0.0, d45, pi], alpha=0.0, eta=0.0),
        dict(id="join_alpha_on", U=[6.0, 10.0, -10.0], A=[0.01, 0.004, 0.01],
             theta=[0.0, d45, pi], alpha=0.2, eta=0.0),
        dict(id="join_alpha_equal_area", U=[6.0, 4.0, -10.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, d45, pi], alpha=0.2, eta=0.0),
        dict(id="join_eta_on", U=[6.0, 4.0, -10.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, d90, pi], alpha=0.0, eta=1.0),
        # A supplier at negative angle, to flip the pseudosupplier direction.
        dict(id="join_negative_angle", U=[6.0, 4.0, -10.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, -d45, pi], alpha=0.0, eta=0.0),
        # --- branches the cases above leave unreached ------------------------
        # An angle below -pi, so the initial wrap's negative-fmod path is taken.
        # Nothing forbids a caller declaring -270 degrees.
        dict(id="div_angle_below_minus_pi", U=[10.0, -6.0, -4.0], A=[0.01, 0.01, 0.01],
             theta=[0.0, pi, -1.5 * pi], alpha=0.0, eta=0.0),
        # The joining asymmetry term with only ONE supplier: it must do
        # nothing, and "div_area_ratio_4" is the alpha=0 twin to compare with.
        dict(id="div_alpha_on_is_inert", U=[10.0, -6.0, -16.0], A=[0.01, 0.01, 0.0025],
             theta=[0.0, pi, d45], alpha=0.2, eta=0.0),
        # The common port at index 0 rather than in the middle or at the end,
        # so the common-port search is not always satisfied on its first look.
        dict(id="join_common_port_first", U=[-10.0, 6.0, 4.0], A=[0.01, 0.01, 0.01],
             theta=[pi, 0.0, d90], alpha=0.0, eta=0.0),
        dict(id="div_common_port_last", U=[-6.0, -4.0, 10.0], A=[0.01, 0.01, 0.01],
             theta=[pi, d90, 0.0], alpha=0.0, eta=0.0),
    ]


def _evaluate(case: dict) -> tuple[np.ndarray, np.ndarray]:
    r = junction_loss_coefficient(
        np.array(case["U"], dtype=float),
        np.array(case["A"], dtype=float),
        np.array(case["theta"], dtype=float),
        joining_etransfer_alpha=case["alpha"],
        eta_scale=case["eta"],
    )
    C = np.zeros(3)
    C[: len(r.C)] = r.C
    K = np.full(2, np.nan)
    if r.K is not None:
        K[: len(r.K)] = np.atleast_1d(r.K)
    return C, K


def _derivatives(case: dict) -> tuple[np.ndarray, np.ndarray]:
    """d(C)/d(U_j) and d(K)/d(U_j), central difference, [quantity][j]."""
    dC = np.zeros((3, 3))
    dK = np.zeros((2, 3))
    for j in range(3):
        hi = dict(case)
        lo = dict(case)
        step = _FD_STEP * max(1.0, abs(case["U"][j]))
        hi["U"] = list(case["U"])
        lo["U"] = list(case["U"])
        hi["U"][j] += step
        lo["U"][j] -= step
        C_hi, K_hi = _evaluate(hi)
        C_lo, K_lo = _evaluate(lo)
        dC[:, j] = (C_hi - C_lo) / (2.0 * step)
        dK[:, j] = (K_hi - K_lo) / (2.0 * step)
    return dC, dK


def _fmt(values) -> str:
    return ", ".join(f"{float(v):.17g}" for v in np.ravel(values))


def main() -> None:
    rows = []
    for case in _cases():
        C, K = _evaluate(case)
        dC, dK = _derivatives(case)
        rows.append(
            "    {"
            f'"{case["id"]}",\n'
            f"     {{{_fmt(case['U'])}}}, {{{_fmt(case['A'])}}}, {{{_fmt(case['theta'])}}},\n"
            f"     {case['alpha']:.17g}, {case['eta']:.17g},\n"
            f"     {{{_fmt(C)}}}, {{{_fmt(K)}}},\n"
            f"     {{{_fmt(dC)}}}, {{{_fmt(dK)}}}}},"
        )

    body = "\n".join(rows)
    _OUT.write_text(
        "// AUTO-GENERATED by validation/junction/data/generate_mynard_reference.py\n"
        "// Do not edit by hand. Golden values for the C++ port of the Mynard\n"
        "// Unified0D junction closure, produced by the Python that ships\n"
        "// (combaero.network._mynard2010) so the gate is equivalence with it.\n"
        "//\n"
        "// dC/dK derivatives are CENTRAL DIFFERENCES of the Python, which is an\n"
        "// independent check on the C++ dual partials: the two are computed by\n"
        "// different methods in different languages.\n"
        "#pragma once\n\n"
        "#include <array>\n\n"
        "namespace combaero::validation::junction {\n\n"
        "struct MynardCase {\n"
        "  const char *id;\n"
        "  std::array<double, 3> u;\n"
        "  std::array<double, 3> area;\n"
        "  std::array<double, 3> theta;\n"
        "  double joining_etransfer_alpha;\n"
        "  double eta_scale;\n"
        "  std::array<double, 3> c;       // per branch; supplier entries zero\n"
        "  std::array<double, 2> k;       // per non-common port; NaN if absent\n"
        "  std::array<double, 9> dc_du;   // [branch][velocity], row-major\n"
        "  std::array<double, 6> dk_du;   // [k entry][velocity], row-major\n"
        "};\n\n"
        f"inline constexpr std::array<MynardCase, {len(rows)}> kMynardCases{{{{\n"
        f"{body}\n}}}};\n\n"
        "} // namespace combaero::validation::junction\n"
    )
    print(f"wrote {_OUT} ({len(rows)} cases)")


if __name__ == "__main__":
    main()
