#!/usr/bin/env python3
"""Generate the C2 bridges of the NASA-9 thermo (#549).

The NASA-9 polynomials are C0 in cp at their joins (1000 K, and 6000 K for six
species) but not C1: dcp/dT jumps 1-4 % (H2 20 %). Every solver Jacobian that
uses dcp/dT steps there.

This replaces cp/R, inside a window [T_j - w, T_j + w] around each join, by a
degree-7 polynomial in t = (T - T_j) / w fixed by eight conditions:

* value, first and second derivative of the NASA polynomial at both window
  edges (Hermite): cp is C2 everywhere;
* the bridge's  integral cp dT  and  integral cp / T dT  across the window
  equal the NASA rise of H/R and S/R across it (including NASA's own tiny h and
  s steps at the join).

So h and s are IDENTICAL to the NASA-9 originals everywhere outside the
windows -- formation enthalpies, absolute entropies, the 298.15 K reference --
and differ only inside a window, by the bridge's local cp deviation.

Input: the NASA-9 coefficients already compiled into
``include/thermo_transport_data.h`` (the generated header; the generator's own
merged_species.json source is not in the repository). Output:
``include/thermo_c2_refit_data.h``. Run from the repository root:

    uv run python thermo_data_generator/generate_c2_refit.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "include" / "thermo_transport_data.h"
OUT = ROOT / "include" / "thermo_c2_refit_data.h"

#: Window half-width as a fraction of the join temperature: 50 K at 1000 K,
#: 300 K at 6000 K. Wide enough that the bridge's curvature stays physical,
#: narrow enough that the cp deviation stays far inside the data's accuracy.
HALF_WIDTH_FRACTION = 0.05

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def parse_nasa9(path: Path = SRC) -> tuple[list[str], list[list[tuple[float, float, list[float]]]]]:
    text = path.read_text()
    names = re.findall(r'"(\w+)"', re.search(r"species_names = \{([^}]*)\}", text).group(1))
    block = text[text.index("nasa_coeffs = {") : text.index("transport_props")]
    species = []
    for line in block.splitlines():
        if not line.startswith("{{{"):
            continue
        ivs = re.findall(r"\{(" + _NUM + r"), (" + _NUM + r"), \{([^}]*)\}\}", line)
        species.append([(float(a), float(b), [float(x) for x in c.split(",")]) for a, b, c in ivs])
    assert len(species) == len(names), (len(species), len(names))
    return names, species


# NASA-9 per interval: a[0..6], b1 = a[8], b2 = a[9]
def cp_R(a: list[float], T: float) -> float:
    return a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4


def dcp_R(a: list[float], T: float) -> float:
    return -2 * a[0] / T**3 - a[1] / T**2 + a[3] + 2 * a[4] * T + 3 * a[5] * T**2 + 4 * a[6] * T**3


def d2cp_R(a: list[float], T: float) -> float:
    return 6 * a[0] / T**4 + 2 * a[1] / T**3 + 2 * a[4] + 6 * a[5] * T + 12 * a[6] * T**2


def H_R(a: list[float], T: float) -> float:
    """H/R [K] = (h/RT) * T."""
    return (
        -a[0] / T
        + a[1] * np.log(T)
        + a[2] * T
        + a[3] * T**2 / 2
        + a[4] * T**3 / 3
        + a[5] * T**4 / 4
        + a[6] * T**5 / 5
        + a[8]
    )


def S_R(a: list[float], T: float) -> float:
    return (
        -a[0] / (2 * T**2)
        - a[1] / T
        + a[2] * np.log(T)
        + a[3] * T
        + a[4] * T**2 / 2
        + a[5] * T**3 / 3
        + a[6] * T**4 / 4
        + a[9]
    )


_GL_X, _GL_W = np.polynomial.legendre.leggauss(32)


def bridge(left: list[float], right: list[float], Tj: float, w: float) -> np.ndarray:
    """Coefficients c[0..7] of cp/R = sum c_k t^k, t = (T - Tj)/w in [-1, 1]."""
    TL, TR = Tj - w, Tj + w
    A = np.zeros((8, 8))
    b = np.zeros(8)
    k = np.arange(8)
    for row, (t, a, T) in enumerate(((-1.0, left, TL), (1.0, right, TR))):
        A[3 * row + 0] = t**k
        A[3 * row + 1] = np.where(k >= 1, k * t ** np.maximum(k - 1, 0), 0.0)
        A[3 * row + 2] = np.where(k >= 2, k * (k - 1) * t ** np.maximum(k - 2, 0), 0.0)
        b[3 * row + 0] = cp_R(a, T)
        b[3 * row + 1] = dcp_R(a, T) * w
        b[3 * row + 2] = d2cp_R(a, T) * w * w
    # integral cp dT = w * integral p dt ;  integral cp/T dT = w * integral p/(Tj + w t) dt
    A[6] = w * np.array([(1 - (-1) ** (j + 1)) / (j + 1) for j in k])
    A[7] = w * np.array([np.sum(_GL_W * _GL_X**j / (Tj + w * _GL_X)) for j in k])
    b[6] = H_R(right, TR) - H_R(left, TL)
    b[7] = S_R(right, TR) - S_R(left, TL)
    return np.linalg.solve(A, b)


def bridges_for(
    intervals: list[tuple[float, float, list[float]]],
) -> list[tuple[float, float, np.ndarray]]:
    out = []
    for (_, tmax, left), (tmin, _, right) in zip(intervals, intervals[1:], strict=False):
        assert tmax == tmin, (tmax, tmin)
        w = HALF_WIDTH_FRACTION * tmax
        out.append((tmax, w, bridge(left, right, tmax, w)))
    return out


def write_header(names: list[str], species: list) -> None:
    lines = [
        "#ifndef THERMO_C2_REFIT_DATA_H",
        "#define THERMO_C2_REFIT_DATA_H",
        "",
        "// !! DO NOT EDIT THIS FILE MANUALLY !!",
        "// Generated by thermo_data_generator/generate_c2_refit.py from the NASA-9",
        "// coefficients in include/thermo_transport_data.h (#549). Regenerate with:",
        "//   uv run python thermo_data_generator/generate_c2_refit.py",
        "//",
        "// Per species, one bridge per NASA-9 join: cp/R = sum c[k] t^k with",
        "// t = (T - T_join) / half_width on |t| < 1, C2 with the NASA polynomials",
        "// at both window edges, and h, s matched to NASA-9 outside the window.",
        "",
        "#include <array>",
        "#include <vector>",
        "",
        "namespace combaero {",
        "",
        f"constexpr double C2_BRIDGE_HALF_WIDTH_FRACTION = {HALF_WIDTH_FRACTION!r};",
        "",
        "struct C2Bridge {",
        "    double T_join;      // [K]",
        "    double half_width;  // [K]",
        "    std::array<double, 8> c;  // cp/R = sum c[k] t^k",
        "};",
        "",
        "const std::vector<std::vector<C2Bridge>> c2_bridges = {",
    ]
    for name, ivs in zip(names, species, strict=True):
        parts = []
        for tj, w, c in bridges_for(ivs):
            coeffs = ", ".join(f"{x:.17g}" for x in c)
            parts.append(f"{{{tj!r}, {w!r}, {{{coeffs}}}}}")
        lines.append(f"    {{{', '.join(parts)}}},  // {name}")
    lines[-1] = lines[-1].replace("},  //", "}  //")
    lines += ["};", "", "}  // namespace combaero", "", "#endif  // THERMO_C2_REFIT_DATA_H", ""]
    OUT.write_text("\n".join(lines))


def main() -> None:
    names, species = parse_nasa9()
    write_header(names, species)
    print(
        f"wrote {OUT.relative_to(ROOT)}: {sum(len(s) - 1 for s in species)} bridges, {len(names)} species"
    )


if __name__ == "__main__":
    main()
