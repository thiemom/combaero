"""The second law refuses a root (#598), at the incompressible closure's resolution.

An element whose ``sigma_norm_high`` is below ``-SIGMA_TOL`` demotes the solve
to INCONSISTENT. On the incompressible closure ``sigma_norm`` is resolved only
to the model's own order ``M^4``: each port carries its energy at its own
constant density, where the gas has ``R T ln(Pt)``. A split that is lossless
in its own model (an ideal tee's ``loss_factor = 0``) sits on the bound, so
the band is widened by that resolution, and on compressible ports not at all.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

from combaero.network import IdealMomentumTeeElement, MomentumChamberNode, NetworkSolver
from combaero.network.second_law import SIGMA_TOL

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test_ideal_tee_element import _branch_net  # noqa: E402

#: The outlets' head [Pa] on the theta = 60 lossless split below.
OUTLET_HEAD = 36.7


def _split(regime: str, gain: float = 0.0) -> dict:
    """The f = 0 split; ``gain`` raises each outlet's static P by that
    fraction of its head -- total pressure created at the junction."""
    net = _branch_net(theta=60.0, loss_factor=0.0)
    for node in net.nodes.values():
        if isinstance(node, MomentumChamberNode):
            node.regime = regime
    if gain:

        class Gain(IdealMomentumTeeElement):
            def row_scale_kinds(self) -> list[str]:
                return super().row_scale_kinds()

            def residuals(self, states, P_jct, port_mdots):
                res, jac = super().residuals(states, P_jct, port_mdots)
                for i in range(self.N):
                    if self._port_signs[i] > 0:
                        res[i] -= gain * OUTLET_HEAD
                return res, jac

            def verify_solution_consistent(self, sol, energy_rel_tol=1e-9) -> bool:
                # the tee's own volumetric check would refuse it first: the
                # second law must do it alone
                return True

        net.elements["jct"].__class__ = Gain
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NetworkSolver(net).solve(timeout=60.0, auto_retry=False)


def test_a_lossless_incompressible_split_passes_at_the_models_resolution() -> None:
    """Its sigma_norm is -2e-7 (outlet Pt 0.06 Pa above the inlet's): the
    densities' 6e-4 spread, not flow work. The band (1.1e-6) covers it."""
    sol = _split("incompressible")
    v = sol["__second_law__"]["jct"]
    assert sol["__success__"], sol["__message__"]
    assert -1e-6 < v["sigma_norm"] < -SIGMA_TOL  # below the bare bound
    assert v["sigma_norm_resolution"] > -v["sigma_norm"]
    assert v["verdict"] == "pass"


def test_compressible_ports_carry_no_band() -> None:
    sol = _split("compressible")
    assert sol["__success__"], sol["__message__"]
    assert "sigma_norm_resolution" not in sol["__second_law__"]["jct"]


@pytest.mark.parametrize("gain", [0.01, 0.1])
def test_total_pressure_created_above_the_resolution_is_refused(gain: float) -> None:
    """Falsification: 1% of the outlets' head created at the junction is
    refused by the second law alone, the tee's own volumetric check switched
    off (it refuses 0.3% itself). Below ~0.5% the band holds it: the model's
    resolution at port Mach 0.033."""
    sol = _split("incompressible", gain)
    assert sol["__converged__"]
    assert sol["__second_law__"]["jct"]["verdict"] == "fail"
    assert not sol["__success__"]
    assert sol["__outcome__"] == "inconsistent"
    assert "second law" in sol["__message__"]


def test_a_refused_lossless_split_at_mach_names_its_loss_factor() -> None:
    """At port Mach 0.27 a loss_factor = 0 split creates flow work (its loss
    is written at low Mach, ~M^2/4 short), and sigma refuses it; the message
    says so, not just a sigma value. At f = 1 the same flow is accepted."""
    from test_ideal_tee_element import _branch_net

    def solve(f: float) -> dict:
        net = _branch_net(m_in=0.4, theta=60.0, regime="compressible", loss_factor=f)
        for node in net.nodes.values():
            if isinstance(node, MomentumChamberNode):
                node.regime = "compressible"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return NetworkSolver(net).solve(timeout=60.0)

    refused = solve(0.0)
    assert refused["__converged__"] and not refused["__success__"]
    assert refused["__second_law__"]["jct"]["verdict"] == "fail"
    assert "loss_factor 0 is below ~M^2/4" in refused["__message__"]
    assert solve(1.0)["__success__"]
