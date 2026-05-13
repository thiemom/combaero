"""Combustor equivalence-ratio sweep using NetworkRunner.

Loads the GUI-saved ``combustor.json`` network (compressible, CH4/air),
injects node labels, and sweeps fuel mass flow to span equivalence ratios
from 0.3 to 0.9.  Produces a T vs. phi plot.

Run::

    uv run python/examples/network_runner_combustion_sweep.py

Set ``COMBAERO_SAVE_PLOTS=1`` to save the figure instead of showing it.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt  # noqa: E402
from plot_utils import show_or_save  # noqa: E402

from gui.backend.runner import NetworkRunner  # noqa: E402

HERE = Path(__file__).parent

# Node IDs in combustor.json (no labels set in the GUI export)
_AIR_INLET = "node_1778698958490"
_FUEL_INLET = "node_1778698961265"
_COMBUSTOR = "node_1778699014636"
_OUTLET = "node_1778699154284"

# Stoichiometric FAR for CH4/air derived from baseline:
# phi=0.518 at m_fuel=0.03 kg/s, m_air=1.0 kg/s
_FAR_STOICH = 0.03 / (1.0 * 0.518)  # ~0.0579

M_AIR = 1.0  # kg/s — baseline air mass flow


def load_runner() -> NetworkRunner:
    """Load combustor.json and inject node labels."""
    with (HERE / "combustor.json").open() as fh:
        schema = json.load(fh)

    label_map = {
        _AIR_INLET: "air_inlet",
        _FUEL_INLET: "fuel_inlet",
        _COMBUSTOR: "combustor",
        _OUTLET: "outlet",
    }
    for node in schema["nodes"]:
        if node["id"] in label_map:
            node["data"]["label"] = label_map[node["id"]]

    return NetworkRunner.from_dict(schema)


def main() -> None:
    runner = load_runner()

    phi_range = np.linspace(0.30, 0.90, 13)
    params = pd.DataFrame({"fuel_inlet.m_dot": phi_range * _FAR_STOICH * M_AIR})

    print(f"Sweeping {len(phi_range)} equivalence ratios …")
    sweep_df = runner.sweep(params, metrics=["combustor.T", "combustor.phi"])

    # phi from solver (more precise than our input estimate)
    phi_solved = sweep_df["combustor.phi"].values
    T_exit = sweep_df["combustor.T"].values
    converged = sweep_df["success"].values

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(phi_solved[converged], T_exit[converged], "o-", label="converged")
    if not converged.all():
        ax.plot(phi_solved[~converged], T_exit[~converged], "x", color="red", label="not converged")
        ax.legend()
    ax.set_xlabel("Equivalence ratio $\\phi$")
    ax.set_ylabel("Combustor exit temperature $T$ [K]")
    ax.set_title("Combustion sweep — compressible CH4/air network")
    ax.grid(True, alpha=0.4)
    fig.tight_layout()

    show_or_save(fig, "network_runner_combustion_sweep.png")


if __name__ == "__main__":
    main()
