"""Shared plot utilities for combaero examples."""

import os
import pathlib

import matplotlib.pyplot as plt


def show_or_save(fig: plt.Figure, filename: str) -> None:
    """Show plot interactively or save to file next to the calling script.

    Decision priority:
      1. COMBAERO_SAVE_PLOTS=1  -> always save (headless CI / Cascade)
      2. COMBAERO_SAVE_PLOTS=0  -> always show (force interactive)
      3. backend == 'agg'       -> save (no display available)
      4. otherwise              -> show interactively

    Args:
        fig:      The matplotlib Figure to show or save.
        filename: Output filename (e.g. 'gas_turbine_cycle.png').
                  Saved next to this module in python/examples/.
    """
    env = os.environ.get("COMBAERO_SAVE_PLOTS", "").strip()
    force_save = env == "1"
    force_show = env == "0"
    headless = plt.get_backend().lower() == "agg"

    if force_show or (not force_save and not headless):
        plt.show()
    else:
        out = pathlib.Path(__file__).parent / filename
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Plot saved to '{out}'")
