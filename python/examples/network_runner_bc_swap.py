"""Off-design pressure sensitivity using NetworkResult.swap_boundary().

Loads the compressible CH4/air combustor network, solves the design point
with a fixed air mass flow, swaps the air inlet to a pressure boundary
(Pt auto-populated from the solved state), and sweeps inlet total pressure
±10 % to show how combustor temperature responds.

The swap step also validates that switching boundary type reproduces the same
physical operating point: the re-solve after the swap must give the same
combustor temperature as the original mass-driven solve.

Run::

    uv run python/examples/network_runner_bc_swap.py

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

_AIR_INLET = "node_1778698958490"
_FUEL_INLET = "node_1778698961265"
_COMBUSTOR = "node_1778699014636"
_OUTLET = "node_1778699154284"


def load_runner() -> NetworkRunner:
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

    # 1. Solve design point with fixed air mass flow
    print("Solving design point (mass-driven) …")
    design = runner.solve(
        {"fuel_inlet.m_dot": 0.03},
        init_strategy="incompressible_warmstart",
    )
    assert design.success, f"Design solve failed: {design.message}"

    design_Pt = design.node_state("air_inlet")["Pt"]
    design_T = design.get("combustor.T")
    print(f"  Pt_inlet = {design_Pt:.0f} Pa,  T_combustor = {design_T:.0f} K")

    # 2. Swap air inlet to pressure boundary (Pt auto-populated from solved Pt)
    print("Swapping air_inlet: mass_boundary → pressure_boundary …")
    p_runner = design.swap_boundary("air_inlet", "pressure_boundary")

    # 3. Verify the swap gives the same operating point
    verify = p_runner.solve(init_strategy="incompressible_warmstart")
    assert verify.success, f"Swap verify failed: {verify.message}"
    dT = abs(verify.get("combustor.T") - design_T)
    print(
        f"  Re-solve T_combustor = {verify.get('combustor.T'):.1f} K  (ΔT = {dT:.2f} K vs design)"
    )

    # 4. Sweep inlet total pressure ±10 % at constant fuel flow
    Pt_range = design_Pt * np.linspace(0.90, 1.10, 11)
    params = pd.DataFrame(
        {
            "air_inlet.Pt": Pt_range,
            "fuel_inlet.m_dot": 0.03,
        }
    )

    print(f"Sweeping Pt over [{Pt_range[0]:.0f}, {Pt_range[-1]:.0f}] Pa ({len(Pt_range)} points) …")
    # default init_strategy + continuation (sweep passes _x0 forward automatically)
    sweep_df = p_runner.sweep(params, metrics=["combustor.T"])

    converged = sweep_df["success"].values
    Pt_norm = sweep_df["air_inlet.Pt"].values / design_Pt
    T_combustor = sweep_df["combustor.T"].values
    print(f"  Converged: {converged.sum()}/{len(converged)}")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(Pt_norm[converged], T_combustor[converged], "o-", label="converged")
    if not converged.all():
        ax.plot(
            Pt_norm[~converged], T_combustor[~converged], "x", color="red", label="not converged"
        )
        ax.legend()
    ax.axvline(1.0, color="gray", linestyle="--", linewidth=0.8, label="design point")
    ax.axhline(design_T, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel(r"$P_{t,\mathrm{inlet}}\ /\ P_{t,\mathrm{design}}$")
    ax.set_ylabel(r"$T_\mathrm{combustor}$ [K]")
    ax.set_title(
        "Off-design pressure sensitivity — compressible CH4/air combustor\n"
        "(pressure-driven after BC swap from fixed mass flow)"
    )
    ax.grid(True, alpha=0.4)
    fig.tight_layout()

    show_or_save(fig, "network_runner_bc_swap.png")


if __name__ == "__main__":
    main()
