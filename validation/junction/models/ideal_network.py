"""The ideal momentum tee (#517) in the junction harness's networks.

Identical to `MPCENetwork` -- same skeletons, axes and extraction -- with the
junction built as `IdealMomentumTeeElement`. That element has no declared flow
direction and no tuning knobs, so the Mynard-specific constructor arguments
are dropped.
"""

from __future__ import annotations

from typing import Any

from combaero.network import IdealMomentumTeeElement
from validation.junction.models.mpce_network import MPCENetwork


class IdealMomentumNetwork(MPCENetwork):
    """`MPCENetwork` with the parameter-free ideal momentum tee."""

    name = "ideal_momentum_network"

    def _junction(self, flow_direction: str, **kw: Any) -> Any:
        return IdealMomentumTeeElement(**kw)
