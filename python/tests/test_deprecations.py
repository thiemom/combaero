"""Every deprecated name warns, and the warning names its removal release.

A deprecation that only lives in a comment or a docstring gives a user no
notice before the name disappears. This is the inventory for 0.7.0: each
entry is called (or accessed) and must raise a DeprecationWarning that names
0.8.0, and the replacement it points to must give the same answer.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import combaero as cb
import combaero.network as net
from combaero.network import pressure_loss

REMOVAL = "removed in 0.8.0"


def _methane_steam() -> np.ndarray:
    X = np.zeros(cb.num_species())
    X[cb.species_index_from_name("CH4")] = 0.10
    X[cb.species_index_from_name("H2O")] = 0.30
    X[cb.species_index_from_name("N2")] = 0.60
    return X


def test_standard_dry_air_composition_warns_and_matches_its_replacement() -> None:
    with pytest.warns(DeprecationWarning, match=REMOVAL):
        old = cb.standard_dry_air_composition()
    assert np.array_equal(np.asarray(old), np.asarray(cb.species.dry_air()))


@pytest.mark.parametrize(
    ("old", "new", "args"),
    [
        ("smr_wgs_equilibrium", "reforming_equilibrium", (1200.0,)),
        ("smr_wgs_equilibrium_adiabatic", "reforming_equilibrium_adiabatic", (900.0,)),
    ],
)
def test_smr_wgs_warns_and_its_replacement_agrees(old: str, new: str, args: tuple[float]) -> None:
    X = _methane_steam()
    with pytest.warns(DeprecationWarning, match=REMOVAL):
        r_old = getattr(cb, old)(*args, X)
    r_new = getattr(cb, new)(*args, X)
    assert np.allclose(r_old.state.X, r_new.state.X, atol=1e-6)


@pytest.mark.parametrize(
    ("alias", "target"),
    [("ConstantLoss", "ConstantFractionLoss"), ("LinearThetaLoss", "LinearThetaFractionLoss")],
)
@pytest.mark.parametrize("module", [net, pressure_loss], ids=["network", "pressure_loss"])
def test_pressure_loss_aliases_warn(module: object, alias: str, target: str) -> None:
    with pytest.warns(DeprecationWarning, match=REMOVAL):
        assert getattr(module, alias) is getattr(pressure_loss, target)


def test_the_aliases_are_not_star_exported() -> None:
    """`from combaero.network import *` must not trigger the warnings."""
    assert "ConstantLoss" not in net.__all__
    assert "LinearThetaLoss" not in net.__all__


def test_an_unknown_network_attribute_still_raises() -> None:
    with pytest.raises(AttributeError):
        net.NoSuchElement  # noqa: B018


def test_the_current_names_do_not_warn() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        cb.species.dry_air()
        cb.reforming_equilibrium(1200.0, _methane_steam())
        assert net.ConstantFractionLoss is pressure_loss.ConstantFractionLoss
