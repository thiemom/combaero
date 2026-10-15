"""Every momentum chamber a GUI network builds runs its global regime (#588).

Drawn or inserted -- a junction's port chambers are not on the canvas -- an
incompressible network must not keep a compressible chamber, and its sonic
limit, hidden inside.
"""

from __future__ import annotations

import pytest
from test_gui_tee_port_step import _schema

from combaero.network import MomentumChamberNode
from gui.backend.graph_builder import build_network_from_schema
from gui.backend.schemas import NetworkGraphSchema


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_every_chamber_takes_the_global_regime(regime: str) -> None:
    schema = NetworkGraphSchema.model_validate(_schema("mynard", 0.05))
    schema.solver_settings.global_regime = regime
    net = build_network_from_schema(schema)
    chambers = [n for n in net.nodes.values() if isinstance(n, MomentumChamberNode)]
    # The tee's three port chambers at least, none of them drawn.
    assert len(chambers) >= 3
    assert {c.regime for c in chambers} == {regime}
