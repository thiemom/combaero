import sys
import warnings
from pathlib import Path

import pytest

import combaero as ca
from combaero.species import SpeciesLocator

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from combustor_liner_cooling_network import solve_operating_point_network_coupled


def test_coupled_combustor_network_converges_without_warning() -> None:
    sp = SpeciesLocator.from_core()

    PR = 20.0
    P2 = PR * 101325.0
    T2 = 667.0
    RH_in = 0.60
    COT_target = 1758.0
    dP_burner_frac = 0.02
    D_ann = 0.30
    u_hot_ann = 25.0

    X_cool = ca.humid_air_composition(288.15, 101325.0, RH_in)
    X_air = ca.species.dry_air()
    X_ch4 = sp.empty()
    X_ch4[sp.indices["CH4"]] = 1.0

    mdot_total = 1.0
    D_ch = 0.008
    L_ch = 0.35
    N_ch = 8

    k_inconel = ca.k_inconel718(900.0)
    wall_layers: list[tuple[float, float]] = [(0.003, k_inconel)]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = solve_operating_point_network_coupled(
            f_cool=0.15,
            mdot_total=mdot_total,
            T2=T2,
            P2=P2,
            X_cool=X_cool,
            X_air=X_air,
            X_ch4=X_ch4,
            COT_target=COT_target,
            dP_burner_frac=dP_burner_frac,
            D_ch=D_ch,
            L_ch=L_ch,
            N_ch=N_ch,
            D_ann=D_ann,
            u_hot_ann=u_hot_ann,
            ch_name="Ribbed",
            wall_layers=wall_layers,
        )

    bad_warnings = [w for w in caught if "did not converge" in str(w.message)]
    assert not bad_warnings, "Coupled network emitted non-convergence warning"

    assert result["h_cool"] > 0.0
    assert result["h_hot"] > 0.0
    assert result["dP_cool"] > 0.0
    assert result["q_wall"] > 0.0

    assert result["T_cool_exit"] > T2
    assert result["FAR"] > 0.0

    assert result["T_hot"] == pytest.approx(COT_target, rel=1e-2)
