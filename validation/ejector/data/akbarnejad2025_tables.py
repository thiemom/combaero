"""Tables 1-3 of Akbarnejad & Ziabasharhagh (2025), Thermal Engineering 72(1).

S. Akbarnejad, M. Ziabasharhagh, "Development of a Theoretical Model for
Predicting Performance of a Gas Ejector", Thermal Engineering 72(1) (2025).
Source PDF: docs/ejector/Theoretical Model for Predicting Performance of a Gas
Ejector.pdf (gitignored).

What these numbers are
----------------------
CFD (ANSYS Fluent 19.2, 2-D axisymmetric, standard k-epsilon, wall functions).
The working fluid is water vapour on both streams, as an IDEAL GAS with
constant specific heat, thermal conductivity and viscosity. They are not
measurements: in the cross-source scorecard they are labelled CFD.

The paper's own model is a scaling law. It takes one reference
characteristic of an ejector (CFD or test) and rescales its critical and
breakdown pressures to new supply conditions. It does not predict either
pressure from geometry, so it supplies no geometry-only P_b0 to adopt.
Its "theoretical model" rows are therefore not used here; only the CFD rows
are.

Not printed in the paper: the stream temperatures. Steam at each stated
pressure is taken as saturated (``T_SAT_C``, steam-table values) by the
scorecard. That is an assumption, and omega is sensitive to it (it follows
the square root of the temperature ratio). The ratio P_b0 / P_c* is far less
sensitive to it.

Transcribed by hand from the PDF text layer (2026-10-19) and checked row by
row against the printed tables.
"""

from __future__ import annotations

#: Table 1, diameters [mm]: primary nozzle throat, nozzle outlet, constant-area
#: section.
GEOMETRY_MM: dict[str, tuple[float, float, float]] = {
    "E1": (3.0, 5.4, 7.0),
    "E2": (26.0, 52.0, 150.0),
}

#: Table 2 (E1) and Table 3 (E2), CFD rows only:
#: (ejector, primary kPa, secondary kPa, entrainment ratio, critical back
#: pressure kPa, breakdown back pressure kPa).
CFD_ROWS: list[tuple[str, float, float, float, float, float]] = [
    ("E1", 400.0, 60.0, 0.365, 83.0, 96.5),
    ("E1", 300.0, 45.0, 0.364, 60.0, 72.5),
    ("E1", 600.0, 45.0, 0.189, 116.0, 123.5),
    ("E1", 800.0, 60.0, 0.188, 157.0, 167.0),
    ("E1", 600.0, 60.0, 0.245, 116.0, 130.0),
    ("E1", 500.0, 60.0, 0.284, 99.0, 113.0),
    ("E1", 600.0, 90.0, 0.369, 128.0, 146.0),
    ("E2", 700.0, 22.85, 1.234, 42.4, 45.3),
    ("E2", 600.0, 22.85, 1.417, 39.2, 41.7),
    ("E2", 800.0, 22.85, 1.080, 45.5, 48.8),
]

#: Saturation temperature [C] of water at each pressure [kPa] in the tables
#: (steam tables); the assumed stream temperatures, see the module docstring.
T_SAT_C: dict[float, float] = {
    22.85: 63.0,
    45.0: 78.7,
    60.0: 85.9,
    90.0: 96.7,
    300.0: 133.5,
    400.0: 143.6,
    500.0: 151.8,
    600.0: 158.8,
    700.0: 165.0,
    800.0: 170.4,
}
