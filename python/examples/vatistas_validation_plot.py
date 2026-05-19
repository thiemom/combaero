"""Vatistas n-vortex model -- validation plots against digitised paper data.

Recreates Figs 1b, 2a, 2b from:
  Vatistas, Kozel, Mih (1991), Exp. Fluids 11, 73-76.
  DOI: 10.1007/BF00198434
"""

import math
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plot_utils import show_or_save

import combaero as cb

DATA_DIR = pathlib.Path(__file__).parent.parent / "tests" / "data" / "vatistas"


def load(filename: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(DATA_DIR / filename)
    return df["x"].to_numpy(), df[" y"].to_numpy()


def main() -> None:
    r_dense = np.linspace(0.0, 4.0, 400)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    fig.suptitle(
        "Vatistas n-vortex model vs. digitised data\n"
        "Vatistas, Kozel, Mih (1991) Exp. Fluids 11, 73-76",
        fontsize=11,
    )

    # ------------------------------------------------------------------
    # Fig 2a -- tangential velocity (n=2)
    # ------------------------------------------------------------------
    ax = axes[0]
    r_exp, v_exp = load("fig2a_exp_data.csv")
    ax.scatter(r_exp, v_exp, s=14, color="steelblue", zorder=3, label="Exp. data (Fig 2a)")
    v_model = np.array([cb.vatistas_v0_bar(r, 2.0) for r in r_dense])
    ax.plot(r_dense, v_model, "k-", lw=1.5, label="Model n=2")
    ax.set_xlabel(r"$\bar{r} = r / r_c$")
    ax.set_ylabel(r"$\overline{V}_0 = V_\theta \cdot 2\pi r_c / \Gamma$")
    ax.set_title("Tangential velocity (Fig 2a)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, None)

    # ------------------------------------------------------------------
    # Fig 2b -- normalised pressure (n=2)
    # ------------------------------------------------------------------
    ax = axes[1]
    r_exp, dp_exp = load("fig2b_exp_data.csv")
    ax.scatter(r_exp, dp_exp, s=14, color="tomato", zorder=3, label="Exp. data (Fig 2b)")
    dp_model = np.array([(2.0 / math.pi) * math.atan(r**2) for r in r_dense])
    ax.plot(r_dense, dp_model, "k-", lw=1.5, label="Model n=2")
    ax.set_xlabel(r"$\bar{r} = r / r_c$")
    ax.set_ylabel(r"$\overline{\Delta P}$  (normalised, 0$\to$1)")
    ax.set_title("Pressure profile (Fig 2b)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1.05)

    # ------------------------------------------------------------------
    # Fig 1b -- radial velocity magnitude (n=1, 2, 3)
    # ------------------------------------------------------------------
    ax = axes[2]
    configs = [
        ("fig_1b_n=1.csv", 1.0, "o", "steelblue"),
        ("fig1b_n=2.csv", 2.0, "s", "tomato"),
        ("fig1b_n=3.csv", 3.0, "^", "seagreen"),
    ]
    for filename, n, marker, color in configs:
        r_exp, vr_exp = load(filename)
        ax.scatter(
            r_exp,
            vr_exp,
            s=14,
            marker=marker,
            color=color,
            zorder=3,
            label=f"Exp. data n={int(n)}",
        )
        vr_model = np.array([abs(cb.vatistas_vr_bar(r, n)) for r in r_dense])
        ax.plot(r_dense, vr_model, color=color, lw=1.5, label=f"Model n={int(n)}")

    ax.set_xlabel(r"$\bar{r} = r / r_c$")
    ax.set_ylabel(r"$|\overline{V}_r|$  (normalised magnitude)")
    ax.set_title("Radial velocity (Fig 1b)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, None)

    plt.tight_layout()
    show_or_save(fig, "vatistas_validation_plot.png")


if __name__ == "__main__":
    main()
