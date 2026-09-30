"""Tests for EffusionPlateElement (#387).

Geometry ground truth is Andrews et al. (1988), ASME 88-GT-290, Table 1,
read off the page rendered at 400 dpi. The composability claim is checked by
building van de Noort and Ireland's (2022) channel topology out of the
element and confirming the coolant flow falls along the wall.
"""

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

# Andrews et al. (1988) Table 1, effusion plate B. X = 0.6 in = 15.24 mm,
# which is what makes N = 1/X^2 come out at the tabulated 4306 m^-2.
ANDREWS_B = {"hole_diameter": 2.16e-3, "wall_thickness": 6.3e-3, "pitch": 15.24e-3}
ANDREWS_C = {"hole_diameter": 3.27e-3, "wall_thickness": 6.3e-3, "pitch": 15.24e-3}


def _air_boundary(name: str, Pt: float, Tt: float = 300.0) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt = Pt
    b.Tt = Tt
    b.Y = cb.species.dry_air_mass()
    return b


class TestGeometryAgainstAndrewsTable1:
    """The hole count must follow from pitch, as the source's own table does."""

    @pytest.mark.parametrize(
        "geom,n_per_m2_ref,a_over_ah_ref",
        [(ANDREWS_B, 4306.0, 5.28), (ANDREWS_C, 4306.0, 3.4)],
    )
    def test_hole_density_and_area_ratio(self, geom, n_per_m2_ref, a_over_ah_ref):
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **geom)

        # N = 1/X^2 for a square array.
        assert e.n_holes / e.panel_area == pytest.approx(n_per_m2_ref, rel=0.002)

        # A/A_h, Andrews' nomenclature: hole approach area over hole internal
        # area. Agreement is ~1.5%, consistent with D and t being nominal.
        a_approach = e.pitch_x * e.pitch_y - math.pi / 4.0 * e.hole_diameter**2
        a_hole = math.pi * e.hole_diameter * e.hole_length
        assert a_approach / a_hole == pytest.approx(a_over_ah_ref, rel=0.02)

    def test_rounding_is_recorded_not_hidden(self):
        # 0.1 x 0.1 m at this pitch is 43.056 holes. The element flows 43 and
        # says so, rather than silently flowing a fractional hole.
        e = EffusionPlateElement("p", "c", "g", panel_area=0.01, **ANDREWS_B)
        assert e.n_holes == 43
        assert e.hole_count_exact == pytest.approx(43.056, rel=1e-3)
        # The whole-hole count implies a slightly different pitch.
        assert e.pitch_actual == pytest.approx(math.sqrt(0.01 / 43), rel=1e-12)
        assert e.pitch_actual != e.pitch_x

        # A case where rounding and truncation DISAGREE -- 43.056 does not
        # distinguish them, so on its own it would let a truncating
        # implementation pass. 43.7 holes is 44, not 43.
        cell = 15.24e-3 * 15.24e-3
        e2 = EffusionPlateElement("p", "c", "g", panel_area=43.7 * cell, **ANDREWS_B)
        assert e2.hole_count_exact == pytest.approx(43.7, rel=1e-6)
        assert e2.n_holes == 44

    def test_total_area_is_the_holes_not_the_panel(self):
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **ANDREWS_B)
        assert e.area == pytest.approx(e.n_holes * math.pi / 4.0 * e.hole_diameter**2)
        assert e.porosity == pytest.approx(e.area / e.panel_area)
        assert 0.0 < e.porosity < 0.05

    def test_inclined_hole_is_longer_than_the_wall_is_thick(self):
        """L = t/sin(alpha). This is why effusion holes are inclined: more
        internal surface for the same wall thickness."""
        straight = EffusionPlateElement(
            "s",
            "c",
            "g",
            hole_diameter=0.5e-3,
            wall_thickness=1.0e-3,
            pitch=3e-3,
            panel_area=0.0025,
        )
        angled = EffusionPlateElement(
            "a",
            "c",
            "g",
            hole_diameter=0.5e-3,
            wall_thickness=1.0e-3,
            pitch=3e-3,
            panel_area=0.0025,
            angle_deg=30.0,
        )
        assert straight.hole_length == pytest.approx(1.0e-3)
        assert angled.hole_length == pytest.approx(2.0e-3)
        assert angled.hole_length / angled.hole_diameter == pytest.approx(4.0)
        # Same wall, same holes: only the drilled length differs.
        assert angled.n_holes == straight.n_holes
        assert angled.area == pytest.approx(straight.area)

    def test_rejects_geometry_that_cannot_be_built(self):
        with pytest.raises(ValueError, match="rounds to none"):
            EffusionPlateElement(
                "p",
                "c",
                "g",
                hole_diameter=1e-4,
                wall_thickness=1e-3,
                pitch=10e-3,
                panel_area=1e-6,
            )
        with pytest.raises(ValueError, match="angle_deg"):
            EffusionPlateElement("p", "c", "g", angle_deg=0.0, panel_area=1.0, **ANDREWS_B)
        with pytest.raises(ValueError, match="pitch_x and pitch_y"):
            EffusionPlateElement(
                "p", "c", "g", hole_diameter=1e-3, wall_thickness=1e-3, panel_area=1.0
            )


class TestCorrelationSelection:
    def test_normed_metering_correlations_are_refused(self):
        """A panel has no pipe and no beta, so an ISO 5167 correlation is not
        merely inaccurate here -- it is undefined."""
        for bad in ("ReaderHarrisGallagher", "Stolz", "Miller"):
            e = EffusionPlateElement("p", "c", "g", panel_area=1.0, correlation=bad, **ANDREWS_B)
            with pytest.raises(ValueError, match="NORMED metering correlation"):
                e.validate()

    def test_defaults_to_the_plenum_fed_wall_correlation(self):
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **ANDREWS_B)
        assert e.correlation == "IdelchikThick"
        e.validate()

    def test_correlation_sees_one_hole_not_the_equivalent_bore(self):
        """The flow equation uses the summed area; the correlation must still
        be handed a single real hole. Conflating them would feed a 2.16 mm
        hole correlation a 100 mm equivalent bore."""
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **ANDREWS_B)
        e.resolve_topology(FlowNetwork())
        assert e._orifice_geom.d == pytest.approx(e.hole_diameter)
        assert e.equivalent_bore > 50.0 * e.hole_diameter
        assert e._orifice_geom.D == 0.0
        assert e.beta == 0.0


class TestHoleReynoldsNumber:
    """The discharge-hole correlations are functions of the HOLE Reynolds
    number. Feeding them the pipe Re_D was wrong by up to 38% for a
    plenum-fed hole, where D_up = 0 froze Re_D at a 1e5 fallback."""

    def test_hole_reynolds_divides_by_the_hole_count(self):
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **ANDREWS_B)
        e.resolve_topology(FlowNetwork())

        class _S:
            T, P, m_dot = 300.0, 1e5, 1.0
            X = cb.species.dry_air()

        re_panel = e._hole_reynolds(_S())
        mu = cb.transport_state(300.0, 1e5, cb.species.dry_air()).mu
        re_expected = 4.0 * (1.0 / e.n_holes) / (math.pi * e.hole_diameter * mu)
        assert re_panel == pytest.approx(re_expected, rel=1e-9)

    def test_a_plain_orifice_has_one_hole(self):
        o = OrificeElement("o", "a", "b", diameter=2e-3, correlation="IdelchikThick")
        assert o._hole_count() == 1.0

    def test_reynolds_responds_to_flow_rather_than_freezing(self):
        """The bug this pins: with no upstream channel, Re_D fell back to a
        constant, so Cd stopped responding to the flow entirely."""
        e = EffusionPlateElement("p", "c", "g", panel_area=1.0, **ANDREWS_B)
        e.resolve_topology(FlowNetwork())

        def _re(m_dot):
            class _S:
                T, P = 300.0, 1e5
                X = cb.species.dry_air()

            s = _S()
            s.m_dot = m_dot
            return e._hole_reynolds(s)

        assert _re(10.0) == pytest.approx(10.0 * _re(1.0), rel=1e-9)
        assert _re(0.1) < _re(1.0) < _re(10.0)


class TestPlateInANetwork:
    def test_plenum_to_plenum_panel_solves_and_conserves(self):
        g = FlowNetwork()
        g.add_node(_air_boundary("coolant", 1.10e5))
        g.add_node(_air_boundary("gas", 1.00e5))
        panel = EffusionPlateElement("panel", "coolant", "gas", panel_area=0.01, **ANDREWS_B)
        g.add_element(panel)
        sol = NetworkSolver(g).solve(method="lm")
        assert sol["__converged__"]
        m = sol["panel.m_dot"]
        assert m > 0.0

        # Against the ideal incompressible discharge through the same area.
        rho = cb.density(300.0, 1.0e5, cb.species.dry_air())
        m_ideal = panel.area * math.sqrt(2.0 * rho * 1.0e4)
        assert 0.4 < m / m_ideal < 1.0

    def test_more_holes_pass_more_flow(self):
        flows = []
        for pitch in (20e-3, 15.24e-3, 10e-3):
            g = FlowNetwork()
            g.add_node(_air_boundary("coolant", 1.10e5))
            g.add_node(_air_boundary("gas", 1.00e5))
            g.add_element(
                EffusionPlateElement(
                    "panel",
                    "coolant",
                    "gas",
                    hole_diameter=2.16e-3,
                    wall_thickness=6.3e-3,
                    pitch=pitch,
                    panel_area=0.01,
                )
            )
            sol = NetworkSolver(g).solve(method="lm")
            assert sol["__converged__"]
            flows.append(sol["panel.m_dot"])
        assert flows[0] < flows[1] < flows[2]

    def test_diagnostics_report_the_effusion_parameters(self):
        g = FlowNetwork()
        g.add_node(_air_boundary("coolant", 1.10e5))
        g.add_node(_air_boundary("gas", 1.00e5))
        panel = EffusionPlateElement("panel", "coolant", "gas", panel_area=0.01, **ANDREWS_B)
        g.add_element(panel)
        sol = NetworkSolver(g).solve(method="lm")
        assert sol["__converged__"]
        d = {k.split(".", 1)[1]: v for k, v in sol.items() if k.startswith("panel.")}

        assert d["n_holes"] == 43.0
        assert d["porosity"] == pytest.approx(panel.porosity)
        assert d["L_over_d"] == pytest.approx(6.3e-3 / 2.16e-3, rel=1e-9)
        # Andrews correlates on coolant mass flow per unit PLATE area.
        assert d["G_coolant"] == pytest.approx(sol["panel.m_dot"] / 0.01, rel=1e-6)
        assert d["is_ingesting"] == 0.0

    def test_ingestion_is_reported_not_averaged_away(self):
        """Gas pressure above coolant pressure means hot gas coming in. A
        homogenised panel must say so; van de Noort's CMF > 0.5 case."""
        g = FlowNetwork()
        g.add_node(_air_boundary("coolant", 1.00e5))
        g.add_node(_air_boundary("gas", 1.05e5))
        g.add_element(EffusionPlateElement("panel", "coolant", "gas", panel_area=0.01, **ANDREWS_B))
        sol = NetworkSolver(g).solve(method="lm")
        if sol["__converged__"]:
            assert sol["panel.is_ingesting"] == 1.0


class TestChannelByComposition:
    """The effusion CHANNEL, built from panels hung off channel segments --
    van de Noort and Ireland's Level 4 with lateral links. No channel-specific
    element exists, and this test is what says none is needed."""

    @staticmethod
    def _build(n_seg: int, gas_pressures: list[float]) -> tuple[FlowNetwork, list[str]]:
        g = FlowNetwork()
        g.add_node(_air_boundary("feed", 1.20e5))
        panels = []
        for i in range(n_seg):
            g.add_node(PlenumNode(f"n{i}"))
            g.add_node(_air_boundary(f"gas{i}", gas_pressures[i]))
            g.add_element(
                ChannelElement(
                    f"seg{i}",
                    "feed" if i == 0 else f"n{i - 1}",
                    f"n{i}",
                    length=0.02,
                    diameter=6e-3,
                )
            )
            pid = f"panel{i}"
            g.add_element(
                EffusionPlateElement(
                    pid,
                    f"n{i}",
                    f"gas{i}",
                    hole_diameter=0.6e-3,
                    wall_thickness=1.0e-3,
                    pitch=3.0e-3,
                    panel_area=0.02 * 0.02,
                    angle_deg=30.0,
                )
            )
            panels.append(pid)
        g.add_node(_air_boundary("dump", 1.02e5))
        g.add_element(ChannelElement("segN", f"n{n_seg - 1}", "dump", length=0.02, diameter=6e-3))
        return g, panels

    def test_coolant_flow_falls_along_the_wall(self):
        """m_channel_in = m_channel_out + m_effusion at every node, so the
        coolant mass flow is f(x) -- the defining feature of the channel."""
        n = 4
        g, panels = self._build(n, [1.00e5] * n)
        sol = NetworkSolver(g).solve(method="lm")
        assert sol["__converged__"]

        seg_flows = [sol[f"seg{i}.m_dot"] for i in range(n)]
        assert all(a > b for a, b in zip(seg_flows, seg_flows[1:], strict=False)), seg_flows

        # And the drop across each node is exactly what that panel bled off.
        for i in range(n - 1):
            bled = seg_flows[i] - seg_flows[i + 1]
            assert bled == pytest.approx(sol[f"{panels[i]}.m_dot"], rel=1e-6)

    def test_a_uniform_gas_side_still_gives_a_non_uniform_wall(self):
        """The coolant side has its own gradient, and it is not small.

        Channel friction plus the mass being bled off drops the coolant
        pressure along the wall, so upstream panels see more drive even when
        the gas pressure is perfectly uniform. Measured here: 114.2 kPa at the
        first node against 104.0 at the last, and the first panel passes ~1.9x
        what the last one does. This is the coolant-side counterpart of van de
        Noort's migration, and it is invisible to a single homogenised panel.
        """
        n = 4
        g, panels = self._build(n, [1.00e5] * n)
        sol = NetworkSolver(g).solve(method="lm")
        assert sol["__converged__"]

        flows = [sol[f"{p}.m_dot"] for p in panels]
        drives = [sol[f"{p}.dP_drive"] for p in panels]

        # Monotone decreasing in both, and a factor approaching two end to end.
        assert all(a > b for a, b in zip(drives, drives[1:], strict=False)), drives
        assert all(a > b for a, b in zip(flows, flows[1:], strict=False)), flows
        assert flows[0] / flows[-1] > 1.5

    def test_a_falling_external_pressure_evens_the_wall_out(self):
        """A falling gas pressure works AGAINST the coolant-side drop.

        Downstream panels lose coolant pressure but also see less back
        pressure, so the two gradients partly cancel and the distribution
        becomes more uniform -- the opposite of the intuition that an external
        gradient always worsens maldistribution. Worth having a test say so,
        because it is a design lever: matching the external gradient to the
        channel loss is what makes the bleed uniform.
        """
        n = 4
        flat, panels = self._build(n, [1.00e5] * n)
        graded, _ = self._build(n, [1.03e5, 1.01e5, 0.99e5, 0.97e5])
        sol_flat = NetworkSolver(flat).solve(method="lm")
        sol_grad = NetworkSolver(graded).solve(method="lm")
        assert sol_flat["__converged__"] and sol_grad["__converged__"]

        def spread(sol):
            f = [sol[f"{p}.m_dot"] for p in panels]
            shares = [x / sum(f) for x in f]
            return max(shares) - min(shares)

        assert spread(sol_grad) < 0.5 * spread(sol_flat)
        # Still not uniform: the compensation is partial, not exact.
        assert spread(sol_grad) > 0.02


class TestInternalHeatTransfer:
    """The coolant-side coefficient, Andrews 86-GT-225 wired into the element.

    The correlations themselves are pinned in
    `tests/test_effusion_internal.cpp` against the paper's algebra, and
    scored against Andrews Fig. 8 by
    `validation/cooling/effusion_internal_runner.py`. What these tests guard
    is the element's geometry bookkeeping -- which areas, which pitch, which
    length -- because that is where the two can silently disagree.
    """

    class _State:
        """Minimal stand-in for a solved NetworkMixtureState."""

        def __init__(self, T, P, X, m_dot):
            self.T, self.P, self.X, self.m_dot = T, P, X, m_dot

    def _panel(self, panel_area=0.01, **overrides):
        geom = {**ANDREWS_C, **overrides}
        return EffusionPlateElement("panel", "cool", "gas", panel_area=panel_area, **geom)

    def _state(self, elem, G, T=300.0):
        return self._State(T, 101325.0, cb.standard_dry_air_composition(), G * elem.panel_area)

    def test_matches_the_runner_that_was_scored_against_figure_eight(self):
        """Two independent paths to the same coefficient must agree.

        The runner builds it from the dataset metadata and the raw
        correlation; the element from its own constructor geometry. They
        agree to 0.02%, the residual being the element's rounding to a whole
        43 holes -- a real and deliberate difference, not drift.
        """
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from validation.cooling import effusion_internal_runner as er
        from validation.cooling.schema import load_dataset

        series = next(s for s in load_dataset() if s.label.endswith("fig8_h_effusionC"))
        elem = self._panel(panel_area=0.01)
        for G in (0.1, 0.4, 1.0, 1.8):
            got = elem.internal_heat_transfer(self._state(elem, G, er.ASSUMED_T))
            assert got["h_plate_area"] == pytest.approx(er.predict(series, G), rel=1e-3), f"G = {G}"

    def test_the_two_areas_are_both_reported_and_differ_by_the_ratio(self):
        """Which area an h belongs to is the thing that goes wrong silently.

        For plate C the hole internal area and the approach area differ by
        3.46, so using one where the other is meant is a factor-of-three
        error. Both are returned, and the identity between them is the check
        that neither has drifted.
        """
        elem = self._panel()
        out = elem.internal_heat_transfer(self._state(elem, 0.5))

        assert out["area_ratio"] == pytest.approx(3.46, abs=0.03)
        assert out["h_hole_area"] == pytest.approx(
            out["h_plate_area"] * out["area_ratio"], rel=1e-12
        )
        # The same heat, whichever pair is used -- which is what makes the
        # two interchangeable only in matched pairs.
        assert out["h_hole_area"] * out["area_hole_total"] == pytest.approx(
            out["h_plate_area"] * out["area_approach_total"], rel=1e-12
        )

    def test_the_sum_is_the_approach_plus_the_throat(self):
        elem = self._panel()
        out = elem.internal_heat_transfer(self._state(elem, 0.5))
        assert out["Nu_internal"] == pytest.approx(out["Nu_approach"] + out["Nu_throat"])
        assert out["Nu_approach"] > 0 and out["Nu_throat"] > 0

    def test_the_approach_term_leads_at_effusion_reynolds_numbers(self):
        """Andrews' headline, at the flows his own rig ran.

        "the hole approach flow heat transfer is much larger than the
        internal hole heat transfer". It is a claim about exponents (0.476
        against 0.8), so it holds at low flow and reverses at high -- which
        is why a throat-only treatment under-predicts and both terms are
        carried.
        """
        elem = self._panel()
        low = elem.internal_heat_transfer(self._state(elem, 0.1))
        assert low["Nu_approach"] > low["Nu_throat"]

        high = elem.internal_heat_transfer(self._state(elem, 5.0))
        assert high["Nu_approach"] < high["Nu_throat"]

    def test_an_inclined_hole_transfers_more_heat(self):
        """The reason effusion holes are inclined: more internal surface for
        the same wall thickness, and a longer entry length."""
        normal = self._panel(angle_deg=90.0)
        angled = self._panel(angle_deg=30.0)
        assert angled.hole_length == pytest.approx(2.0 * normal.hole_length)

        a = angled.internal_heat_transfer(self._state(angled, 0.5))
        n = normal.internal_heat_transfer(self._state(normal, 0.5))
        assert a["area_hole_total"] > n["area_hole_total"]
        # More surface at the same coefficient scale means more heat, even
        # though the longer hole has a LOWER entry-length enhancement.
        assert (
            a["h_plate_area"] * a["area_approach_total"]
            > n["h_plate_area"] * n["area_approach_total"]
        )
        assert cb.mills_entry_length_factor(
            angled.hole_length / angled.hole_diameter
        ) < cb.mills_entry_length_factor(normal.hole_length / normal.hole_diameter)

    def test_no_flow_means_no_coefficient_rather_than_a_zero(self):
        """A zero would read as "computed, and it is nothing". Nothing was
        computed, so nothing is reported -- the same discipline the scorecard
        uses for an unscored series."""
        elem = self._panel()
        assert elem.internal_heat_transfer(self._state(elem, 0.0)) == {}

    def test_it_rises_with_coolant_flow(self):
        elem = self._panel()
        previous = None
        for G in (0.05, 0.2, 0.6, 1.5):
            v = elem.internal_heat_transfer(self._state(elem, G))["h_plate_area"]
            if previous is not None:
                assert v > previous, f"G = {G}"
            previous = v

    def test_diagnostics_surfaces_it_after_a_solve(self):
        """The coefficient has to reach the user, not just exist."""
        net = FlowNetwork()
        net.add_node(_air_boundary("cool", 120e3))
        net.add_node(_air_boundary("gas", 100e3))
        net.add_element(EffusionPlateElement("panel", "cool", "gas", panel_area=0.01, **ANDREWS_C))
        result = NetworkSolver(net).solve()
        m_dot = result["panel.m_dot"]
        assert m_dot > 0.0, "the panel should blow, not ingest"

        diag = result["__element_diag__"]["panel"]
        for key in (
            "Re_hole",
            "Nu_approach",
            "Nu_throat",
            "h_hole_area",
            "h_plate_area",
            "area_ratio",
        ):
            assert key in diag, f"{key} missing from diagnostics"
        assert diag["h_plate_area"] > 0.0

        # It must describe the flow the solver actually found rather than a
        # default: rebuild Re from the reported mass flow.
        elem = net.elements["panel"]
        mu = cb.transport_state(300.0, 120e3, cb.standard_dry_air_composition()).mu
        assert diag["Re_hole"] == pytest.approx(
            4.0 * (m_dot / elem.n_holes) / (math.pi * elem.hole_diameter * mu),
            rel=0.05,
        )

    def test_the_approach_term_uses_the_drilled_length_not_the_wall_thickness(
        self,
    ):
        """Eq. (18)'s X/(pi L) is the DRILLED length, not the wall thickness.

        Added because a falsification found nothing: swapping `hole_length`
        for `wall_thickness` in the approach term passed every other test in
        this class, since they all use 90-degree holes where the two are
        equal. At 30 degrees they differ by a factor of two, so the
        substitution would double the approach Nusselt number silently.

        Checked against the correlation called directly, which is the only
        way to pin WHICH length reaches it.
        """
        elem = self._panel(angle_deg=30.0)
        assert elem.hole_length == pytest.approx(2.0 * elem.wall_thickness)

        out = elem.internal_heat_transfer(self._state(elem, 0.5))
        expected = cb.effusion_approach_nusselt(
            out["Re_hole"], 0.7272, elem.pitch_actual / elem.hole_length
        )
        assert out["Nu_approach"] == pytest.approx(expected, rel=1e-3)

        # And it is NOT what the wall thickness would give -- the thing the
        # other tests could not distinguish.
        wrong = cb.effusion_approach_nusselt(
            out["Re_hole"], 0.7272, elem.pitch_actual / elem.wall_thickness
        )
        assert wrong == pytest.approx(2.0 * expected, rel=1e-3)
        assert out["Nu_approach"] != pytest.approx(wrong, rel=0.1)

        # The throat term likewise: R_Nu is on the drilled L/D.
        assert out["Nu_throat"] == pytest.approx(
            cb.effusion_throat_nusselt(
                out["Re_hole"], 0.7272, elem.hole_length / elem.hole_diameter
            ),
            rel=1e-3,
        )
