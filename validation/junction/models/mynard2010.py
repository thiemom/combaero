"""
Mynard & Valen-Sendstad 2010 Unified0D junction-loss model.

The canonical implementation lives in `combaero.network._mynard2010` --
`MPCEv2Element`/`ConstantKTeeElement` depend on it at runtime, so it must
ship inside the `combaero` wheel rather than this dev-only `validation/`
tree (see CHANGELOG: PyPI install of combaero-gui 0.4.0 hit
`ModuleNotFoundError: No module named 'validation'` before this move).
This module re-exports it so existing validation scripts/tests keep
their import path unchanged.
"""

from __future__ import annotations

from combaero.network._mynard2010 import (
    MynardResult,
    junction_loss_coefficient,
)

__all__ = ["MynardResult", "junction_loss_coefficient"]
