"""Chen et al. 2013 Figs. 5a and 7a (validation/ejector/data/chen2013_figs.py).

The digitization's own checks, kept so a re-digitization cannot drift silently:
Fig. 5a's measured points ARE Henry 2007 Fig. 5's data (a third channel),
the P_p reading of each rig follows from the efficiency-free ratio test, and
the entrainment ratio (Huang Eqs. 1-8) reproduces Chen's own model plateaus --
a second author's evaluation of the same equations (fidelity).
"""

from __future__ import annotations

import pytest

from combaero.network._ejector_huang1999 import EjectorGeometry, entrainment_ratio
from validation.ejector.data import chen2013_figs as chen
from validation.ejector.data import henry2007_fig5 as henry

T_K = 293.15
ATM_PA = 101325.0


def _plateau(fig: list, kind: str, pp: float) -> float:
    pts = sorted((r for r in fig if r[1] == kind and r[2] == pp), key=lambda r: r[3])
    return pts[0][4]


def test_fig5a_measured_points_are_henrys_data() -> None:
    measured = [r for r in chen.FIG5A if r[1] == "measured"]
    assert len(measured) == 34
    for rid, _, pp, pc, om, _ in measured:
        hx, hy = min(henry.FIG5_OMEGA_VS_PB[int(pp)], key=lambda q: abs(q[0] - pc))
        assert abs(pc - hx) < 0.013 and abs(om - hy) < 0.011, rid


def _ratio_rms(geom: EjectorGeometry, measured: dict[float, float], offset: float) -> float:
    pps = sorted(measured)
    om = {pp: entrainment_ratio(pp * 1e5 + offset, T_K, ATM_PA, T_K, geom, 1.4).omega for pp in pps}
    d = [om[pp] / om[pps[0]] - measured[pp] / measured[pps[0]] for pp in pps]
    return (sum(v * v for v in d) / len(d)) ** 0.5


@pytest.mark.parametrize(
    ("fig", "reading"), [("7a", "absolute"), ("5a", "gauge")], ids=["rect-absolute", "henry-gauge"]
)
def test_the_primary_pressure_reading_follows_from_plateau_ratios(fig: str, reading: str) -> None:
    """Ratios of measured plateau omega across P_p cancel every efficiency
    coefficient; only the pressure reading moves them. Chen's rectangular rig
    reads absolute, the Hemidi/Henry rig gauge (Henry's "effective pressure"),
    though Chen ran his model on the latter at absolute pressure."""
    a1, a2 = chen.AREA_RATIOS[fig]
    geom = EjectorGeometry(area_ratio_nozzle=a1, area_ratio_mix=a2)
    data = chen.FIG7A if fig == "7a" else chen.FIG5A
    measured = {pp: _plateau(data, "measured", pp) for pp in sorted({r[2] for r in data})}
    absolute, gauge = _ratio_rms(geom, measured, 0.0), _ratio_rms(geom, measured, ATM_PA)
    assert (absolute < gauge / 2) if reading == "absolute" else (gauge < absolute / 2)


@pytest.mark.parametrize("fig", ["5a", "7a"])
def test_the_entrainment_ratio_reproduces_chens_model_plateaus(fig: str) -> None:
    """Fidelity of Eqs. 1-8 on a second author's own outputs: Chen's model
    plateaus at P_p absolute, suction 1 bar (the reading his model used), to
    within 1-2.5% (a uniform +1.5%: his property values are not printed)."""
    a1, a2 = chen.AREA_RATIOS[fig]
    geom = EjectorGeometry(area_ratio_nozzle=a1, area_ratio_mix=a2)
    data = chen.FIG5A if fig == "5a" else chen.FIG7A
    for pp in sorted({r[2] for r in data}):
        om = entrainment_ratio(pp * 1e5, T_K, 1.0e5, T_K, geom, 1.4).omega
        assert 0.0 < om / _plateau(data, "model", pp) - 1.0 < 0.025, pp


def test_the_printed_chain_misses_both_authors_own_critical_pressure() -> None:
    """The documented fidelity gap (README "Why the printed chain misses its
    authors' own numbers"): Huang's printed Eqs. 9-18 land 16-29% below Huang's
    own T_c*-derived P_c* and 9-17% below the plateau end of every Chen model
    curve. If an implementation change closes it, this test says so."""
    import json
    from pathlib import Path

    from validation.ejector import printed_chain_check as pcc

    cases = json.loads(Path(pcc.FIXTURE).read_text())["cases"]
    for c in cases:
        i = c["inputs"]
        geom = EjectorGeometry(i["area_ratio_nozzle"], i["area_ratio_mix"])
        pc = pcc.printed_chain(
            i["p_g_pa"], i["t_g_k"], i["p_e_pa"], i["t_e_k"], geom, i["gamma"], i["r_gas"]
        )
        assert -0.30 < pc / c["paper_theory"]["p_c_pa"] - 1.0 < -0.15, c["id"]
    for fig, pp, geom, lo, hi, _ in pcc.chen_curves():
        pc = pcc.printed_chain(pp * 1e5, T_K, 1.0e5, T_K, geom, 1.4, 287.05) / 1e5
        assert -0.18 < pc / (0.5 * (lo + hi)) - 1.0 < -0.08, (fig, pp)
