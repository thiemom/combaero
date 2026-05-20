"""
Validation of the Vatistas n-vortex model against digitised experimental data
from the original paper (Vatistas et al., Exp. Fluids 11, 73-76, 1991).

Data files live in python/tests/data/vatistas/.  Axes:
  fig2a_exp_data.csv : x = r_bar, y = V0_bar (normalised tangential velocity)
  fig2b_exp_data.csv : x = r_bar, y = delta_P_bar (normalised pressure, 0->1)
  fig_1b_n=1.csv     : x = r_bar, y = |Vr_bar| (magnitude of normalised radial velocity)
  fig1b_n=2.csv      : same, n=2
  fig1b_n=3.csv      : same, n=3

Tolerance strategy: 90th-percentile of absolute errors rather than max, to avoid
failures from isolated digitisation noise.
"""

import math
import pathlib

import numpy as np
import pandas as pd
import pytest

import combaero as cb

DATA_DIR = pathlib.Path(__file__).parent / "data" / "vatistas"


def load(filename: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(DATA_DIR / filename)
    return df["x"].to_numpy(), df[" y"].to_numpy()


# ---------------------------------------------------------------------------
# Fig 2a -- tangential velocity profile (n=2, experimental scatter)
# ---------------------------------------------------------------------------


def test_fig2a_tangential_velocity():
    r_bar, v0_bar_data = load("fig2a_exp_data.csv")
    v0_bar_model = np.array([cb.vatistas_v0_bar(rb, 2.0) for rb in r_bar])
    abs_err = np.abs(v0_bar_model - v0_bar_data)
    p90 = np.percentile(abs_err, 90)
    assert p90 < 0.15, (
        f"Fig 2a: 90th-percentile absolute error {p90:.4f} exceeds 0.15 "
        f"(mean={abs_err.mean():.4f}, max={abs_err.max():.4f})"
    )


# ---------------------------------------------------------------------------
# Fig 2b -- pressure profile (n=2, experimental scatter)
# ---------------------------------------------------------------------------


def test_fig2b_pressure_profile():
    r_bar, dp_bar_data = load("fig2b_exp_data.csv")
    # Normalised: delta_P_bar = I(r_bar) / (pi/4) = (2/pi)*arctan(r_bar^2) for n=2
    dp_bar_model = np.array([(2.0 / math.pi) * math.atan(rb**2) for rb in r_bar])
    abs_err = np.abs(dp_bar_model - dp_bar_data)
    p90 = np.percentile(abs_err, 90)
    assert p90 < 0.15, (
        f"Fig 2b: 90th-percentile absolute error {p90:.4f} exceeds 0.15 "
        f"(mean={abs_err.mean():.4f}, max={abs_err.max():.4f})"
    )


# ---------------------------------------------------------------------------
# Fig 1b -- radial velocity profiles (n=1, 2, 3)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "filename, n",
    [
        ("fig_1b_n=1.csv", 1.0),
        ("fig1b_n=2.csv", 2.0),
        ("fig1b_n=3.csv", 3.0),
    ],
)
def test_fig1b_radial_velocity(filename: str, n: float):
    r_bar, vr_bar_data = load(filename)
    # Inside the core (r_bar < 1) all n curves overlap and digitisation is unreliable;
    # restrict comparison to r_bar >= 1 where the curves are well separated.
    mask = r_bar >= 1.0
    r_bar, vr_bar_data = r_bar[mask], vr_bar_data[mask]
    # CSV stores |Vr_bar| (positive magnitude); formula gives negative values
    vr_bar_model = np.array([abs(cb.vatistas_vr_bar(rb, n)) for rb in r_bar])
    abs_err = np.abs(vr_bar_model - vr_bar_data)
    p90 = np.percentile(abs_err, 90)
    assert p90 < 0.30, (
        f"Fig 1b (n={n}): 90th-percentile absolute error {p90:.4f} exceeds 0.30 "
        f"(mean={abs_err.mean():.4f}, max={abs_err.max():.4f})"
    )
