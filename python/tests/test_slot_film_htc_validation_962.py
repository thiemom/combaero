"""#962 validation of the film-side h near a slot: three separate questions.

FIDELITY (does the implementation mirror the paper?) on Seban's OWN data:
the digitised Fig. 12 points (m > 1, x/s < 40; signed off 2026-11-03,
validation/cooling/extractions/seban_1960_slot_film_h.md) against the
implemented curve, Fig. 12's ordinate (u_s s/nu_s)^0.3 St_s x 100. They are
the points Seban drew curve B through. 38 points, of which only 5 are
cleanly sampled (26 stacked, 7 ambiguous): an overplotted figure, so the
pooled score is an upper bound on fidelity. Read uncertainty 0.004 dec
(clean) to 0.010 dec (stacked), about 1-2 %: the 11 % RMS is the runs'
scatter about one curve (different m and s), not the reading. Curve C
(x/s > 40) is not scored: no point beyond x/s 40 was digitised.

ACCURACY (another rig), CROSS-SOURCE: the Cranfield SME4 Fig. 1 points
(Hartnett, Birkebak & Eckert 1961; Pai & Whitelaw 1971; m < 1.3, other
labs), against curve B as Nu = 0.135 Pr [Re_s x/s]^0.7 at Pr 0.69. Most of
these runs are at m < 1, where the shipped element hands back to the duct's
h: this scores the curve's constant, not what the element does there.

METHOD COMPARISON (not fidelity for our curve): the book's film-cooled and
TBC worked examples (Lefebvre & Ballal 3rd ed. pp. 338-340, 350-351) as
combaero networks with #960's radiation. The book uses the m > 1.3 wall jet
(Eq. 8.44) at m = 2.54, x/s 18; with its own printed h (1926) the chain
reproduces its answers (test_liner_radiation_793.py). Here: the named
option (Eq. 8.44 on combaero's slot-air properties) and the shipped Seban
curve. No constant is tuned to any of this.
"""

from __future__ import annotations

import csv
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    FlowNetwork,
    GasRadiation,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    SlotFilmElement,
    ThermalWall,
    WallLayer,
)
from validation.cooling.schema import DATA_ROOT

sys.path.insert(0, str(Path(__file__).parent))
import test_liner_radiation_793 as book  # noqa: E402

PR_SEBAN = 0.69  # the Prandtl number Seban's 0.10 x Pr = 0.069 used


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def _seban() -> list[dict[str, str]]:
    return _rows(DATA_ROOT / "seban1960" / "fig12_h_m_gt_1.csv")


def _fig12_ordinate(x_over_s: float, C_near: float = cb.SEBAN1960_C_NEAR) -> float:
    """(u_s s / nu_s)^0.3 St_s x 100: Re_s drops out of Seban's form."""
    C_far = cb.SEBAN1960_C_FAR if C_near == cb.SEBAN1960_C_NEAR else cb.SEBAN1960_C_FAR_PRINTED
    return 100.0 * cb.slot_film_stanton_seban_1960(x_over_s, x_over_s, 1.0, C_near, C_far).St


def _log_errors(rows: list[dict[str, str]], C_near: float = cb.SEBAN1960_C_NEAR) -> np.ndarray:
    return np.array(
        [
            math.log10(float(r["y_group"]) / _fig12_ordinate(float(r["x_over_s"]), C_near))
            for r in rows
        ]
    )


# --- fidelity on Seban's own Fig. 12 ----------------------------------------


def test_the_digitised_set_is_the_signed_off_one() -> None:
    rows = _seban()
    assert len(rows) == 38
    q = [r["quality"] for r in rows]
    assert (q.count("clean"), q.count("stacked"), q.count("ambiguous")) == (5, 26, 7)
    assert max(float(r["x_over_s"]) for r in rows) < 40.0


def test_fidelity_on_seban_fig12() -> None:
    """FIDELITY, Seban's own data, upper bound (overplotted). The drawn
    curve: bias -1.0 % (geometric mean), RMS 0.046 dec (11 %), worst 33 %,
    38 points; -1.2 % without the 7 ambiguous reads; +3.5 % on the 5 clean
    points alone. The same points sit 34 % above the printed Eq. 13 (0.10): the
    misprint, and why the drawn 0.135 ships."""
    e = _log_errors(_seban())
    assert 10 ** e.mean() - 1.0 == pytest.approx(-0.0105, abs=1e-3)
    assert math.sqrt((e**2).mean()) == pytest.approx(0.0456, abs=5e-4)
    assert 10 ** np.abs(e).max() - 1.0 == pytest.approx(0.330, abs=2e-3)
    unamb = _log_errors([r for r in _seban() if r["quality"] != "ambiguous"])
    assert 10 ** unamb.mean() - 1.0 == pytest.approx(-0.0124, abs=1e-3)
    clean = _log_errors([r for r in _seban() if r["quality"] == "clean"])
    assert 10 ** clean.mean() - 1.0 == pytest.approx(0.0353, abs=1e-3)
    printed = _log_errors(_seban(), cb.SEBAN1960_C_NEAR_PRINTED)
    assert 10 ** printed.mean() - 1.0 == pytest.approx(0.336, abs=2e-3)


def test_the_seban_robustness_table_brackets_the_shipped_constant() -> None:
    """Every data-based variant of the curve-B constant (pooling, stacks,
    leave-one-station-out, bootstrap) lies in 0.1296-0.1339, against the
    drawn 0.135; the printed 0.10 and run 128 alone (m 1.10) are the outliers."""
    rows = _rows(DATA_ROOT / "seban1960" / "near_slot_constant_robustness.csv")
    seban = [r for r in rows if r["figure"] == "Seban"]
    data = [
        float(r["C_fit"])
        for r in seban
        if not r["variant"].startswith(("printed", "run 128", "curve B", "clean"))
    ]
    assert len(data) >= 25
    assert min(data) > 0.129 and max(data) < 0.140
    drawn = next(float(r["C_fit"]) for r in seban if r["variant"].startswith("curve B"))
    assert drawn == pytest.approx(cb.SEBAN1960_C_NEAR, abs=1e-3)


# --- cross-source accuracy on SME4 Fig. 1 ------------------------------------


def _sme4_log_errors(C: float, source: str = "") -> np.ndarray:
    pts = _rows(DATA_ROOT / "ballal_lefebvre1973" / "sme4_fig1_near_slot_nu.csv")
    return np.array(
        [
            math.log10(float(p["Nu"]) / (C * float(p["Re_s_x_over_s"]) ** 0.7))
            for p in pts
            if p["source_ref"].startswith(source)
        ]
    )


def test_cross_source_accuracy_on_sme4_fig1() -> None:
    """ACCURACY, CROSS-SOURCE (Hartnett et al. 1961 and Pai & Whitelaw 1971,
    61 points, m < 1.3). Curve B at Pr 0.69: bias -2.5 %, RMS 17.8 %, 32 of 61
    within 10 %; by source +6.9 % (Hartnett, 20) and -6.7 % (Pai &
    Whitelaw, 41). These data lie 32 % above the book's 0.069 (Eq. 8.41) and 1.0 %
    above SME4's own 0.090."""
    C = cb.SEBAN1960_C_NEAR * PR_SEBAN
    e = _sme4_log_errors(C)
    assert len(e) == 61
    assert 10 ** e.mean() - 1.0 == pytest.approx(-0.025, abs=1e-3)
    assert 10 ** math.sqrt((e**2).mean()) - 1.0 == pytest.approx(0.178, abs=2e-3)
    assert int(np.sum(np.abs(10**e - 1.0) <= 0.10)) == 32
    assert 10 ** _sme4_log_errors(C, "ref3").mean() - 1.0 == pytest.approx(0.069, abs=2e-3)
    assert 10 ** _sme4_log_errors(C, "ref7").mean() - 1.0 == pytest.approx(-0.067, abs=2e-3)
    assert 10 ** _sme4_log_errors(0.069).mean() - 1.0 == pytest.approx(0.317, abs=2e-3)
    assert 10 ** _sme4_log_errors(0.090).mean() - 1.0 == pytest.approx(0.010, abs=2e-3)


# --- the book's worked examples as networks: a method comparison -----------

S_SLOT = 0.00145  # p. 338: x/s 18, t/s 0.4, s 1.45 mm
A_SLOT = 5.95e-4  # pi D_L s
X_FILM = 18.0 * S_SLOT  # 0.0261 m
M_SLOT = 0.289  # kg/s, rho_a U_a A_s
A_LINER = 0.0137  # m^2
METAL = [WallLayer(0.0012, 26.0)]
TBC = [WallLayer(0.0005, 2.6), WallLayer(0.0012, 26.0)]


def _book_net(model: str, layers: list[WallLayer]) -> FlowNetwork:
    """The film-cooled liner of pp. 338-340: products at 2280 K, the slot's
    0.289 kg/s of 880 K air, the wall at x/s 18 (the book's point value),
    the annulus and radiation of the uncooled example (pp. 325-328). The
    mainstream is 2.62 - 0.289 kg/s so that the liner's mass flux is the
    book's 191 kg/(m^2 s) and m = 2.54 as printed. Film effectiveness is
    combaero's Marek-Tacina at Tu 14 % (0.801; the book's Eq. 8.38 gives
    0.789, and m 2.54 is outside Marek-Tacina's runs)."""
    Lseg = 0.01
    A = math.pi * 0.132 * Lseg
    D_L = math.sqrt(4.0 * A_LINER / math.pi)
    Yp = book._products(0.05)
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("g_in", m_dot=2.62 - M_SLOT, Tt=book.TG, Y=Yp))
    n.add_node(MassFlowBoundary("c_sup", m_dot=M_SLOT, Tt=book.T3, Y=book.AIR))
    n.add_node(PlenumNode("mix"))
    n.add_node(PressureBoundary("g_out", Pt=book.P3, Tt=book.TG, Y=Yp))
    n.add_node(MassFlowBoundary("a_in", m_dot=7.074, Tt=book.T3, Y=book.AIR))
    n.add_node(PressureBoundary("a_out", Pt=book.P3, Tt=book.T3, Y=book.AIR))
    n.add_element(ChannelElement("g_feed", "g_in", "mix", length=0.01, diameter=D_L))
    n.add_element(
        SlotFilmElement(
            "slot",
            "c_sup",
            "mix",
            slot_height=S_SLOT,
            span=A_SLOT / S_SLOT,
            lip_thickness=0.4 * S_SLOT,
            film_htc=model,
        )
    )
    n.add_element(
        ChannelElement(
            "liner",
            "mix",
            "g_out",
            length=Lseg,
            diameter=D_L,
            surface=ConvectiveSurface(area=A, model=SmoothModel()),
        )
    )
    n.add_element(
        ChannelElement(
            "annulus",
            "a_in",
            "a_out",
            length=Lseg,
            diameter=math.sqrt(4 * 0.01476 / math.pi),
            Dh=0.0576,
            surface=ConvectiveSurface(area=math.pi * 0.1344 * Lseg, model=SmoothModel()),
        )
    )
    n.add_wall(
        ThermalWall(
            "w",
            "liner",
            "annulus",
            layers=layers,
            contact_area=A,
            film="slot",
            film_start=X_FILM,
            film_end=X_FILM,
            radiation_a=GasRadiation(
                beam_length=book.L_B,
                emissivity_wall=book.EPS_W,
                luminosity=book.LUM,
                T_gas=book.TG,
                fuel_air_ratio=book.Q_PZ,
            ),
            radiation_b=CasingRadiation(
                area_ratio=book.AREA_RATIO,
                emissivity_wall=book.EPS_W,
                emissivity_casing=book.EPS_C,
            ),
        )
    )
    return n


def _solve(n: FlowNetwork) -> dict:
    s = NetworkSolver(n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return r


@pytest.mark.parametrize(
    ("model", "layers", "expected", "h"),
    [
        # Printed: 1283/1265 (film), 1304/1236/1220 (TBC).
        ("ballal_lefebvre1973", METAL, (1269.2, 1251.6), 2067.6),
        ("ballal_lefebvre1973", TBC, (1288.3, 1224.0, 1208.6), 2067.6),
        ("seban1960", METAL, (1285.1, 1266.7), 1632.3),
        ("seban1960", TBC, (1308.1, 1240.5, 1224.3), 1632.3),
        # Without a film-side h (the default): the duct's 509 W/m^2K.
        ("none", METAL, (1358.1, 1336.1), 509.4),
    ],
)
def test_book_film_and_tbc_examples_as_networks(
    model: str, layers: list[WallLayer], expected: tuple, h: float
) -> None:
    """METHOD COMPARISON, not fidelity. Against the printed 1283/1265 and
    1304/1236/1220 K:

    - Eq. 8.44 (named option) on combaero's slot air: h 2068 against the
      book's 1926 (k_a 0.060 against the book's 0.0553, Re_s 1.78e4 against
      1.81e4), 14 K cooler.
    - Seban (shipped): h 1632, 15 % below the book, Re_s 1.78e4 above his
      runs (in_box 0). With the book's T_aw it would be 10 K hotter; the
      network's Marek-Tacina eta (T_aw 1158 against 1176) and annulus h (911
      against 921) take 8 K back, so it lands 2 K above print (4 K with the
      TBC).
    - Without a film-side h the wall is 73 K hotter than with Seban's.
    """
    r = _solve(_book_net(model, layers))
    assert r["w.film_blowing_ratio"] == pytest.approx(2.54, abs=0.005)
    assert r["w.h_a"] == pytest.approx(h, abs=0.5)
    for got, want in zip(r["w.T_interface"], expected, strict=True):
        assert got == pytest.approx(want, abs=0.5)
    if model == "seban1960":
        assert r["w.film_h_Re_in_box"] == 0.0 and r["w.film_h_in_box"] == 0.0
        assert r["w.film_h_x_in_box"] == 1.0 and r["w.film_h_M_in_box"] == 1.0


def test_book_printed_inputs_isolate_the_h() -> None:
    """At the book's own T_aw (1176 K) and wall chain, only h differs: Eq.
    8.44 with the book's printed k_a and Re_s gives its 1926 and its answer;
    on combaero's properties 2068 (-4.4 K); Seban's 1632 (+10.2 K)."""
    n = cb.slot_film_nusselt_ballal_lefebvre_1973(18.0, 18.0, 1.81e4)
    assert n.Nu_s * 0.0553 / S_SLOT == pytest.approx(1926.0, abs=1.0)
    X_air = cb.species.dry_air()
    G = M_SLOT / A_SLOT
    out = {}
    for name, kernel in (
        ("bl", cb.slot_film_htc_ballal_lefebvre_1973),
        ("seban", cb.slot_film_htc_seban_1960),
    ):
        h = kernel(X_FILM, X_FILM, S_SLOT, G, 2.54, 0.0, book.T3, book.P3, X_air).h
        w = book._book_wall(METAL)
        w.compute_coupling(h, book.T_AW_FILM, 1.0, book.H2, book.T3, 1.0)
        out[name] = (h, w._last_profile[0])
    assert out["bl"][0] == pytest.approx(2067.5, abs=0.5)
    assert out["seban"][0] == pytest.approx(1632.3, abs=0.5)
    assert out["bl"][1] == pytest.approx(1278.6, abs=0.2)
    assert out["seban"][1] == pytest.approx(1293.2, abs=0.2)
