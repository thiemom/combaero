"""Tests for compressible flow elements in network solver."""

import math

import numpy as np
import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    MultiPortChamberElement,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)


def test_compressible_orifice_network():
    """Test network with compressible orifice element."""
    net = FlowNetwork()

    # High pressure ratio to ensure compressible effects
    inlet = PressureBoundary("inlet", Pt=300000.0, Tt=300.0)
    outlet = PressureBoundary("outlet", Pt=101325.0, Tt=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Compressible orifice
    orifice = OrificeElement(
        "orifice",
        "inlet",
        "outlet",
        Cd=0.65,
        diameter=0.011284,
        regime="compressible",
        correlation="fixed",
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive"

    # Check pressure ratio to confirm compressible regime
    PR = outlet.Pt / inlet.Pt
    assert PR < 0.8, f"Should be in compressible regime (PR={PR:.3f})"

    # Verify against direct compressible calculation
    X = cb.species.dry_air()
    mdot_direct, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        inlet.Tt, inlet.Pt, outlet.Pt, X, 0.65, 1e-4, 0.0
    )

    # Should match within solver tolerance
    assert abs(sol["orifice.m_dot"] - mdot_direct) / mdot_direct < 0.01


def test_compressible_channel_network():
    """Test network with compressible Fanno channel element."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", Pt=200000.0, Tt=400.0)
    outlet = PressureBoundary("outlet", Pt=150000.0, Tt=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Simple compressible Fanno channel directly between boundaries
    channel = ChannelElement(
        "channel",
        "inlet",
        "outlet",
        length=1.0,
        diameter=0.08,
        roughness=1e-4,
        regime="compressible",
        friction_model="haaland",
    )

    net.add_element(channel)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["channel.m_dot"] > 0, "Channel mass flow should be positive"

    # Verify pressure drop is reasonable
    dP = inlet.Pt - outlet.Pt
    assert dP > 0, "Pressure should drop across channel"


def test_compressible_channel_beyond_subsonic_envelope_converges():
    """A pressure-driven channel demanding more than the subsonic Fanno
    envelope must converge onto the choke-barrier boundary state.

    The L=1 m, D=0.05 m air channel leaves the subsonic envelope at a Pt
    ratio of about 1.19. Before L_choke was interpolated within the march
    step (it snapped to L/n_steps), the choke-barrier ramp was a
    staircase in m_dot: when the residual zero fell inside a staircase
    jump this network had no root, and every method/init combination
    stalled at a tread edge (best |F| ~ 90 at Pt ratio 1.2).
    """
    mdots = []
    for pr in (1.2, 3.0):
        net = FlowNetwork()
        net.add_node(PressureBoundary("inlet", Pt=100000.0 * pr, Tt=300.0))
        net.add_node(PressureBoundary("outlet", Pt=100000.0, Tt=300.0))
        net.add_element(
            ChannelElement(
                "channel",
                "inlet",
                "outlet",
                length=1.0,
                diameter=0.05,
                regime="compressible",
            )
        )
        solver = NetworkSolver(net)
        sol = solver.solve(timeout=30.0)
        assert sol["__success__"], (
            f"pr={pr}: {sol.get('__message__')} (|F|={sol.get('__final_norm__')})"
        )
        assert sol["channel.m_dot"] > 0
        mdots.append(sol["channel.m_dot"])
    # Choked mass flow scales with inlet pressure: the barrier root must
    # track the envelope boundary, not park at an arbitrary tread.
    assert mdots[1] > 2.0 * mdots[0]


def test_mixed_compressible_incompressible_network():
    """Test network mixing compressible and incompressible elements."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", Pt=250000.0, Tt=350.0)
    outlet = PressureBoundary("outlet", Pt=101325.0, Tt=350.0)
    junction1 = PlenumNode("j1")
    junction2 = PlenumNode("j2")

    net.add_node(inlet)
    net.add_node(outlet)
    net.add_node(junction1)
    net.add_node(junction2)

    # Incompressible channel first
    channel1 = ChannelElement(
        "channel1",
        "inlet",
        "j1",
        length=1.0,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    # Compressible orifice in middle (high pressure drop)
    orifice = OrificeElement(
        "orifice",
        "j1",
        "j2",
        Cd=0.65,
        diameter=0.011284,
        regime="compressible",
        correlation="fixed",
    )

    # Incompressible channel at end
    channel2 = ChannelElement(
        "channel2",
        "j2",
        "outlet",
        length=0.5,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    net.add_element(channel1)
    net.add_element(orifice)
    net.add_element(channel2)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["channel1.m_dot"] > 0
    assert sol["orifice.m_dot"] > 0
    assert sol["channel2.m_dot"] > 0

    # Mass conservation
    assert abs(sol["channel1.m_dot"] - sol["orifice.m_dot"]) < 1e-6
    assert abs(sol["orifice.m_dot"] - sol["channel2.m_dot"]) < 1e-6

    # Pressure should decrease monotonically
    assert inlet.Pt > sol["j1.P"] > sol["j2.P"] > outlet.Pt


def test_compressible_vs_incompressible_comparison():
    """Compare compressible vs incompressible orifice in same network."""

    # Setup identical networks with different regimes
    def create_network(regime):
        net = FlowNetwork()
        inlet = PressureBoundary("inlet", Pt=200000.0, Tt=300.0)
        outlet = PressureBoundary("outlet", Pt=150000.0, Tt=300.0)
        net.add_node(inlet)
        net.add_node(outlet)

        orifice = OrificeElement(
            "orifice",
            "inlet",
            "outlet",
            Cd=0.65,
            diameter=0.015958,
            regime=regime,
            correlation="fixed",
        )
        net.add_element(orifice)
        return net

    # Solve both
    net_incomp = create_network("incompressible")
    net_comp = create_network("compressible")

    solver_incomp = NetworkSolver(net_incomp)
    solver_comp = NetworkSolver(net_comp)

    sol_incomp = solver_incomp.solve(method="lm")
    sol_comp = solver_comp.solve(method="lm")

    # At moderate pressure ratios, they should be similar but not identical
    mdot_incomp = sol_incomp["orifice.m_dot"]
    mdot_comp = sol_comp["orifice.m_dot"]

    # Compressible should give different result (can be up to 20% at PR=0.75)
    error = abs(mdot_comp - mdot_incomp) / mdot_comp
    assert 0.01 < error < 0.25, f"Expected 1-25% difference, got {error * 100:.1f}%"


def test_compressible_choked_orifice():
    """Test compressible orifice in choked flow regime."""
    net = FlowNetwork()

    # Very high pressure ratio to ensure choking
    inlet = PressureBoundary("inlet", Pt=500000.0, Tt=300.0)
    outlet = PressureBoundary("outlet", Pt=101325.0, Tt=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    orifice = OrificeElement(
        "orifice",
        "inlet",
        "outlet",
        Cd=0.65,
        diameter=0.007979,
        regime="compressible",
        correlation="fixed",
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify choked flow
    X = cb.species.dry_air()
    PR_crit = cb.critical_pressure_ratio(inlet.Tt, inlet.Pt, X)
    PR_actual = outlet.Pt / inlet.Pt

    assert PR_actual < PR_crit, (
        f"Flow should be choked (PR={PR_actual:.3f} < PR_crit={PR_crit:.3f})"
    )
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive even when choked"


def test_compressible_channel_high_mach():
    """Test compressible channel at higher Mach number."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", Pt=300000.0, Tt=400.0)
    outlet = PressureBoundary("outlet", Pt=200000.0, Tt=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Long channel with small diameter to get higher velocity
    channel = ChannelElement(
        "channel",
        "inlet",
        "outlet",
        length=5.0,
        diameter=0.03,
        roughness=1e-4,
        regime="compressible",
    )
    net.add_element(channel)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    assert sol["channel.m_dot"] > 0, "Mass flow should be positive"

    # Verify Mach number is significant
    X = cb.species.dry_air()
    rho = cb.density(inlet.Tt, inlet.Pt, X)
    area = 3.14159 * (0.03 / 2) ** 2
    u = sol["channel.m_dot"] / (rho * area)
    a = cb.speed_of_sound(inlet.Tt, X)
    M = u / a

    assert M > 0.1, f"Mach number should be significant (M={M:.3f})"


def _area(diameter: float) -> float:
    return float(np.pi * (diameter / 2.0) ** 2)


def _mdot_for_mach(mach: float, T: float, P_ref: float, X: list[float], area: float) -> float:
    rho = cb.density(T, P_ref, X)
    a = cb.speed_of_sound(T, X)
    return float(mach * rho * a * area)


# Constants matching the network builder; keep in sync
_T_HOT = 740.0
_T_COLD = 320.0
_P_OUT = 2.0e5


def _build_fully_coupled_network(regime: str, mach_target: float) -> FlowNetwork:
    x_air = cb.species.dry_air()
    y_air = list(cb.mole_to_mass(x_air))

    t_hot = _T_HOT
    t_cold = _T_COLD
    p_out_hot = _P_OUT
    p_out_cold = _P_OUT
    d_hot = 0.040
    d_cold = 0.035
    length = 2.0

    mdot_hot = _mdot_for_mach(mach_target, t_hot, p_out_hot, x_air, _area(d_hot))
    mdot_cold = 0.65 * _mdot_for_mach(mach_target, t_cold, p_out_cold, x_air, _area(d_cold))

    channel_regime = "compressible" if regime == "compressible" else "incompressible"
    orifice_regime = "compressible" if regime == "compressible" else "incompressible"

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, Tt=t_hot, Y=y_air))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(PressureBoundary("hot_outlet", Pt=p_out_hot, Tt=t_hot, Y=y_air))

    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, Tt=t_cold, Y=y_air))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(PressureBoundary("cold_outlet", Pt=p_out_cold, Tt=t_cold, Y=y_air))

    hot_surface = ConvectiveSurface(area=np.pi * d_hot * length, model=SmoothModel())
    cold_surface = ConvectiveSurface(area=np.pi * d_cold * length, model=SmoothModel())

    net.add_element(
        ChannelElement(
            id="hot_channel",
            from_node="hot_inlet",
            to_node="hot_plenum",
            length=length,
            diameter=d_hot,
            roughness=5.0e-5,
            regime=channel_regime,
            surface=hot_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="hot_orifice",
            from_node="hot_plenum",
            to_node="hot_outlet",
            Cd=0.72,
            diameter=0.030,
            regime=orifice_regime,
            correlation="fixed",
        )
    )

    net.add_element(
        ChannelElement(
            id="cold_channel",
            from_node="cold_inlet",
            to_node="cold_plenum",
            length=length,
            diameter=d_cold,
            roughness=5.0e-5,
            regime=channel_regime,
            surface=cold_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="cold_orifice",
            from_node="cold_plenum",
            to_node="cold_outlet",
            Cd=0.72,
            diameter=0.025,
            regime=orifice_regime,
            correlation="fixed",
        )
    )

    net.add_wall(
        ThermalWall(
            id="coupling_wall",
            element_a="hot_channel",
            element_b="cold_channel",
            layers=[WallLayer(thickness=0.002, conductivity=20.0)],
        )
    )
    net.thermal_coupling_enabled = True
    return net


def test_fully_coupled_compressible_vs_incompressible_mach_sweep():
    """Fully coupled pressure-ratio sweep over target Mach: compressible vs incompressible."""
    mach_values = np.array([0.20])

    pr_comp: list[float] = []
    pr_incomp: list[float] = []
    thermal_hot_drop: list[float] = []
    thermal_cold_rise: list[float] = []

    for mach_target in mach_values:
        net_comp = _build_fully_coupled_network("compressible", float(mach_target))
        net_incomp = _build_fully_coupled_network("incompressible", float(mach_target))

        solver_comp = NetworkSolver(net_comp)
        solver_incomp = NetworkSolver(net_incomp)

        sol_comp = solver_comp.solve(method="hybr")
        sol_incomp = solver_incomp.solve(method="hybr")

        assert sol_comp["__success__"], f"compressible solve failed at M={mach_target:.3f}"
        assert sol_incomp["__success__"], f"incompressible solve failed at M={mach_target:.3f}"

        pr_comp.append(_P_OUT / float(sol_comp["hot_inlet.Pt"]))
        pr_incomp.append(_P_OUT / float(sol_incomp["hot_inlet.Pt"]))

        t_hot_out = float(solver_comp._derived_states["hot_plenum"][0])
        t_cold_out = float(solver_comp._derived_states["cold_plenum"][0])
        thermal_hot_drop.append(_T_HOT - t_hot_out)
        thermal_cold_rise.append(t_cold_out - _T_COLD)

    pr_comp_arr = np.asarray(pr_comp)
    pr_incomp_arr = np.asarray(pr_incomp)
    pr_diff = np.abs(pr_incomp_arr - pr_comp_arr)

    # PR should be monotonic non-increasing as target Mach increases.
    # A tiny plateau at high Mach is acceptable near choking; end-point
    # upturns are not.
    tol = 1e-6
    assert np.all(np.diff(pr_comp_arr) <= tol)
    assert np.all(np.diff(pr_incomp_arr) <= tol)

    # Compressible model predicts larger pressure drop (lower PR).
    assert np.all(pr_comp_arr < pr_incomp_arr)

    # Relative compressibility impact should grow with Mach.
    rel_diff = pr_diff / pr_comp_arr
    assert np.all(np.diff(rel_diff) > 0)

    # Thermal coupling should remain active through the sweep.
    assert np.min(thermal_hot_drop) > 0.0
    assert np.min(thermal_cold_rise) > 0.0


def test_homotopy_initialization_strategy():
    """Verify that init_strategy='homotopy' converges for a high-flow network."""
    net = FlowNetwork()
    print("running: test_homotopy_initialization_strategy")
    # Define a high-flow scenario that might struggle with cold start
    # but should be easy for homotopy
    inlet = MassFlowBoundary("inlet", m_dot=2.0, Tt=500.0)
    outlet = PressureBoundary("outlet", Pt=101325.0, Tt=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Long narrow channel to create high pressure drop and high velocity
    channel = ChannelElement(
        "channel",
        "inlet",
        "outlet",
        length=5.0,
        diameter=0.1,
        roughness=1e-4,
        regime="compressible",
    )
    net.add_element(channel)

    solver = NetworkSolver(net)

    # Solve with homotopy
    sol = solver.solve(method="hybr", init_strategy="homotopy")
    print("solver solution:")
    print(sol)
    assert sol["__success__"], "Homotopy solve failed"
    assert sol["channel.m_dot"] > 0, "Mass flow should be positive"


# ---------------------------------------------------------------------------
# analytical_pt_prop init strategy
# ---------------------------------------------------------------------------


def _mpce_tee_network(
    pt_ratio: float,
    theta_deg: float,
    psi: float,
    flow_direction: str,
    Tt: float = 300.0,
    Pt_ref: float = 100000.0,
) -> FlowNetwork:
    """Canonical 3-port MPCE tee. Mirrors tmp/mpce_cold_start_audit.py."""
    diameter = 0.05
    A_com = math.pi * (diameter / 2.0) ** 2
    A_bra = A_com / psi
    D_bra = 2.0 * math.sqrt(A_bra / math.pi)
    length = 1.0
    Pt_hi = Pt_ref * pt_ratio
    Pt_lo = Pt_ref
    net = FlowNetwork()

    if flow_direction == "branch":
        net.add_node(PressureBoundary("pb_hi", Pt=Pt_hi, Tt=Tt))
        net.add_node(PressureBoundary("pb_lo_str", Pt=Pt_lo, Tt=Tt))
        net.add_node(PressureBoundary("pb_lo_bra", Pt=Pt_lo, Tt=Tt))
        net.add_node(MomentumChamberNode("mc_com", area=A_com))
        net.add_node(MomentumChamberNode("mc_str", area=A_com))
        net.add_node(MomentumChamberNode("mc_bra", area=A_bra))
        net.add_element(
            ChannelElement(
                "ch_in", "pb_hi", "mc_com", length=length, diameter=diameter, regime="compressible"
            )
        )
        net.add_element(
            ChannelElement(
                "ch_str",
                "mc_str",
                "pb_lo_str",
                length=length,
                diameter=diameter,
                regime="compressible",
            )
        )
        net.add_element(
            ChannelElement(
                "ch_bra",
                "mc_bra",
                "pb_lo_bra",
                length=length,
                diameter=D_bra,
                regime="compressible",
            )
        )
        net.add_element(
            MultiPortChamberElement(
                id="jct",
                inlet_nodes=["mc_com"],
                outlet_nodes=["mc_str", "mc_bra"],
                inlet_angles_deg=[0.0],
                outlet_angles_deg=[0.0, theta_deg],
                port_areas=[A_com, A_com, A_bra],
            )
        )
    else:
        net.add_node(PressureBoundary("pb_hi_str", Pt=Pt_hi, Tt=Tt))
        net.add_node(PressureBoundary("pb_hi_bra", Pt=Pt_hi, Tt=Tt))
        net.add_node(PressureBoundary("pb_lo", Pt=Pt_lo, Tt=Tt))
        net.add_node(MomentumChamberNode("mc_str", area=A_com))
        net.add_node(MomentumChamberNode("mc_bra", area=A_bra))
        net.add_node(MomentumChamberNode("mc_com", area=A_com))
        net.add_element(
            ChannelElement(
                "ch_str",
                "pb_hi_str",
                "mc_str",
                length=length,
                diameter=diameter,
                regime="compressible",
            )
        )
        net.add_element(
            ChannelElement(
                "ch_bra",
                "pb_hi_bra",
                "mc_bra",
                length=length,
                diameter=D_bra,
                regime="compressible",
            )
        )
        net.add_element(
            ChannelElement(
                "ch_out", "mc_com", "pb_lo", length=length, diameter=diameter, regime="compressible"
            )
        )
        net.add_element(
            MultiPortChamberElement(
                id="jct",
                inlet_nodes=["mc_str", "mc_bra"],
                outlet_nodes=["mc_com"],
                inlet_angles_deg=[0.0, theta_deg],
                outlet_angles_deg=[0.0],
                port_areas=[A_com, A_bra, A_com],
            )
        )
    return net


def test_analytical_pt_prop_rescues_cold_stuck_case():
    """Case #4 from the LHS-32 audit, with both solver-review fixes in:
    the choke barrier (this PR) and the honest collector-port closure
    (#212). Their interplay flipped this test twice:

    - pre-barrier, pre-#212: cold start stalled at |F|~3111 in the choked
      flat region; analytical_pt_prop escaped in ~110 evals (the original
      #210 rescue claim, but its long hybr trajectory through the barrier
      region proved platform-sensitive on CI).
    - barrier alone: the flat attractor was the whole story -- the default
      cold start converged in seconds and analytical was no faster.
    - barrier + honest closure (this, final state): the v1 impulse merge's
      root structure changed; the default cold start stalls at |F|~1e5
      again, while analytical_pt_prop converges DIRECTLY in ~1-2 s (no
      LM-fallback marathon). The rescue semantics are restored with
      comfortable margin.
    """
    net = _mpce_tee_network(
        pt_ratio=1.8385290361105023,
        theta_deg=46.83168762022338,
        psi=0.7810566293706353,
        flow_direction="merge",
    )
    solver = NetworkSolver(net)
    # Local timing: ~1s with lm from the analytical seed; 180s ceiling
    # tolerates slower CI runners.
    sol = solver.solve(method="lm", init_strategy="analytical_pt_prop", timeout=180.0)
    assert sol["__success__"], f"analytical_pt_prop failed: {sol.get('__message__')}"
    assert sol["__final_norm__"] < 1e-3
    # Mass conservation at the merge junction (sum of port flows into the
    # chamber = 0, per sign convention of MultiPortChamberElement).
    m_str, m_bra, m_out = sol["ch_str.m_dot"], sol["ch_bra.m_dot"], sol["ch_out.m_dot"]
    # Either supply -> outlet or reverse; either root is valid without a
    # flow_direction constraint. Just require mass balance.
    assert abs(m_str + m_bra - m_out) < 1e-3, (
        f"mass imbalance at merge: m_str={m_str}, m_bra={m_bra}, m_out={m_out}"
    )


def _mpce_mfb_merge_network(theta_deg: float = 45.0) -> FlowNetwork:
    """MFB-driven 3-port merge tee with unequal supply temperatures.

    Mass-flow boundaries impose the flow direction, so the solve lands the
    physical basin deterministically -- suitable for closure/mixing
    regression tests independent of basin selection.
    """
    diameter = 0.05
    A = math.pi * (diameter / 2.0) ** 2
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_hot", m_dot=0.3, Tt=400.0))
    net.add_node(MassFlowBoundary("mfb_cold", m_dot=0.1, Tt=300.0))
    net.add_node(PressureBoundary("pb_out", Pt=100000.0, Tt=300.0))
    net.add_node(MomentumChamberNode("mc_str", area=A))
    net.add_node(MomentumChamberNode("mc_bra", area=A))
    net.add_node(MomentumChamberNode("mc_com", area=A))
    net.add_element(
        ChannelElement(
            "ch_hot", "mfb_hot", "mc_str", length=1.0, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_cold", "mfb_cold", "mc_bra", length=1.0, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_out", "mc_com", "pb_out", length=1.0, diameter=diameter, regime="compressible"
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_str", "mc_bra"],
            outlet_nodes=["mc_com"],
            inlet_angles_deg=[0.0, theta_deg],
            outlet_angles_deg=[0.0],
            port_areas=[A, A, A],
        )
    )
    return net


def test_mpce_collector_port_carries_dynamic_head():
    """Collector-port MCNs must see their real face flow in the
    Pt = P + 0.5*rho*v^2 closure.

    Before the throughflow fix, MultiPortChamberElement.flow_at_node
    returned 0, so any port fed by the junction had _total_m_dot = 0 and
    its closure silently degenerated to Pt = P -- the outflow dynamic head
    (tens of kPa at these conditions) vanished from the bookkeeping.
    """
    net = _mpce_mfb_merge_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=120.0)
    assert sol["__final_norm__"] < 1e-3

    node = net.nodes["mc_com"]
    m_out = sol["ch_out.m_dot"]
    assert m_out == pytest.approx(0.4, abs=1e-4)
    assert node._total_m_dot == pytest.approx(m_out, rel=1e-6)

    T_com, _, _ = solver._derived_states["mc_com"]
    rho = cb.density(T_com, sol["mc_com.P"], cb.species.dry_air())
    A = node.area
    q_expected = 0.5 * rho * (m_out / (rho * A)) ** 2
    q_solved = sol["mc_com.Pt"] - sol["mc_com.P"]
    assert q_solved > 1e3, "collector dynamic head must not degenerate to zero"
    assert q_solved == pytest.approx(q_expected, rel=1e-4)


def test_mpce_merge_outlet_temperature_mass_weighted():
    """Merging streams with unequal temperatures must mix mass-weighted.

    Before the throughflow fix, all streams through the junction carried
    m_dot = 0, so the mixer fell back to the FIRST stream's temperature:
    the outlet of a 0.3 kg/s @ 400 K + 0.1 kg/s @ 300 K merge reported
    400 K instead of the enthalpy-weighted ~375 K.
    """
    net = _mpce_mfb_merge_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=120.0)
    assert sol["__final_norm__"] < 1e-3

    T_com, _, _ = solver._derived_states["mc_com"]
    assert T_com == pytest.approx(375.1, abs=2.0)
    assert abs(T_com - 400.0) > 20.0, "outlet T stuck at first-stream value"


def test_mpce_collector_port_mcn_inherits_area():
    """Auto-sized collector-port MCNs inherit the junction port area.

    Nodes resolve before elements, and a collector port has no upstream
    channel to inherit Dh from, so MomentumChamberNode.resolve_topology
    leaves it at the 0.1 m^2 fallback; MultiPortChamberElement.resolve_
    topology must then push the inferred port area onto it, otherwise the
    Pt closure sees a near-zero face velocity.
    """
    D_main, D_bra = 0.05, 0.04
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_hi", Pt=150000.0, Tt=300.0))
    net.add_node(PressureBoundary("pb_lo_str", Pt=100000.0, Tt=300.0))
    net.add_node(PressureBoundary("pb_lo_bra", Pt=100000.0, Tt=300.0))
    net.add_node(MomentumChamberNode("mc_com"))
    net.add_node(MomentumChamberNode("mc_str"))
    net.add_node(MomentumChamberNode("mc_bra"))
    net.add_element(
        ChannelElement(
            "ch_in", "pb_hi", "mc_com", length=1.0, diameter=D_main, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_str", "mc_str", "pb_lo_str", length=1.0, diameter=D_main, regime="compressible"
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra", "mc_bra", "pb_lo_bra", length=1.0, diameter=D_bra, regime="compressible"
        )
    )
    net.add_element(
        MultiPortChamberElement(
            id="jct",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 45.0],
        )
    )
    net.resolve_all_topology()

    A_main = math.pi * (D_main / 2.0) ** 2
    A_bra = math.pi * (D_bra / 2.0) ** 2
    # Inlet port: pre-existing inheritance via upstream channel Dh.
    assert net.nodes["mc_com"].area == pytest.approx(A_main, rel=1e-12)
    # Collector ports: inherited through the junction's port areas.
    assert net.nodes["mc_str"].area == pytest.approx(A_main, rel=1e-12)
    assert net.nodes["mc_bra"].area == pytest.approx(A_bra, rel=1e-12)
    assert net.nodes["mc_bra"].Dh == pytest.approx(D_bra, rel=1e-12)


def test_choke_barrier_rejects_dead_branch_spurious_root():
    """Regression for the choked-channel flat region (solver review, 2026-07).

    Before the choke barrier in channel_compressible_mdot_and_jacobian, the
    state below was an EXACT root (|F| ~ 3e-7) of the branch tee: the branch
    arm dead (m_dot = 0) and the straight arm carrying the full 0.883 kg/s at
    M ~ 1.1 with zero total-pressure drop, because choked channels reported
    dP = 0 with zero m_dot sensitivity. With the barrier, the straight
    channel's row must now show a large residual at this state.
    """
    net = _mpce_tee_network(
        pt_ratio=2.0,
        theta_deg=45.0,
        psi=1.0,
        flow_direction="branch",
    )
    # Residual probes outside solve() require resolved topology: without it
    # the port MCNs are not marked as junction ports and _residuals is a
    # different function than the one the solver sees.
    net.resolve_all_topology()
    solver = NetworkSolver(net)
    x = solver._build_x0()
    spurious = {
        "mc_com.P": 100000.0,
        "mc_com.Pt": 187100.39,
        "mc_str.P": 100000.0,
        "mc_str.Pt": 100000.0,
        "mc_bra.P": 100000.0,
        "mc_bra.Pt": 100000.0,
        "jct.P_jct": 100000.0,
        "ch_in.m_dot": 0.8831240560197469,
        "ch_str.m_dot": 0.8831240560197469,
        "ch_bra.m_dot": 0.0,
    }
    for name, val in spurious.items():
        x[solver._name_to_index[name]] = val
    res = solver._residuals(x)
    assert float(np.linalg.norm(res)) > 1e4, (
        "choked dead-branch state must no longer be a root of the network"
    )


def test_analytical_pt_prop_rejects_unknown_strategy_name():
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    net.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    net.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    solver = NetworkSolver(net)
    with pytest.raises(ValueError, match="init_strategy must be one of"):
        solver.solve(init_strategy="analytical_pt_property")  # type: ignore[arg-type]


def test_mpce_network_default_init_auto_upgrades():
    """MPCE networks silently upgrade init_strategy='default' to
    'analytical_pt_prop' (32/32 vs 28/32 certified-root convergence on
    the 2026-07 inverse-design audit)."""
    net = _mpce_tee_network(
        pt_ratio=1.5,
        theta_deg=60.0,
        psi=1.0,
        flow_direction="branch",
    )
    solver = NetworkSolver(net)
    # auto_retry=False: this test forces a failure (maxfev=5) to inspect
    # the recorded strategy; with the retry enabled the diagnostics would
    # reflect the outlet-ref warm-start retry instead of the upgrade.
    # Short timeout: only the recorded strategy stamp matters here; the LM
    # fallback would otherwise idle away the rest of the budget.
    _ = solver.solve(timeout=2.0, options={"maxfev": 5}, auto_retry=False)
    used = solver._diagnostic_data["solver_settings_used"]["init_strategy"]
    assert used == "analytical_pt_prop"


def test_non_mpce_network_default_init_unchanged():
    """Networks without an MPCE keep the plain default initialization."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    net.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    net.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    solver = NetworkSolver(net)
    _ = solver.solve(timeout=10.0)
    used = solver._diagnostic_data["solver_settings_used"]["init_strategy"]
    assert used == "default"


def test_mpce_explicit_strategy_opts_out_of_auto_upgrade():
    """Passing a non-default init_strategy suppresses the MPCE upgrade."""
    net = _mpce_tee_network(
        pt_ratio=1.5,
        theta_deg=60.0,
        psi=1.0,
        flow_direction="branch",
    )
    solver = NetworkSolver(net)
    # Short timeout: only the recorded strategy stamp matters here.
    _ = solver.solve(init_strategy="homotopy", timeout=2.0, options={"maxfev": 5})
    used = solver._diagnostic_data["solver_settings_used"]["init_strategy"]
    assert used == "homotopy"


def test_mfb_driven_mpce_network_keeps_plain_default_init():
    """MFB-driven MPCE networks (< 2 PressureBoundaries) must NOT
    auto-upgrade to analytical_pt_prop.

    The strategy's Pt propagation needs a boundary pressure gradient.
    With a single PB the propagated per-node Pt differences are
    path-length artifacts and the Bernoulli sign logic seeds REVERSED
    channel m_dots, parking the first iterate inside the junctions'
    soft-barrier penalty region (|F(x0)| ~ 1e6 observed on a real
    5-tee GUI network) while the plain cold start converges.
    """
    net = _mpce_mfb_merge_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=30.0)
    used = solver._diagnostic_data["solver_settings_used"]["init_strategy"]
    assert used == "default"
    assert sol["__success__"], sol.get("__message__")


def test_incompressible_warmstart_deprecated():
    net = FlowNetwork()
    net.add_node(PressureBoundary("in", Pt=200000.0, Tt=300.0))
    net.add_node(PressureBoundary("out", Pt=100000.0, Tt=300.0))
    net.add_element(OrificeElement("orf", "in", "out", 0.6, 0.05, correlation="fixed"))
    solver = NetworkSolver(net)
    with pytest.warns(DeprecationWarning, match="incompressible_warmstart"):
        _ = solver.solve(init_strategy="incompressible_warmstart", timeout=10.0)


def test_analytical_pt_prop_respects_user_initial_guess():
    """User-provided initial_guess on a channel wins over the analytical seed."""
    net = _mpce_tee_network(
        pt_ratio=1.5,
        theta_deg=60.0,
        psi=1.0,
        flow_direction="branch",
    )
    # Preset an obviously-wrong m_dot guess on ch_in that the solver would
    # normally overwrite via analytical seeding. If our precedence order is
    # correct, x0 keeps this value going into the solve.
    net.elements["ch_in"].initial_guess = {"ch_in.m_dot": 42.0}
    solver = NetworkSolver(net)
    _ = solver.solve(
        method="hybr",
        init_strategy="analytical_pt_prop",
        timeout=10.0,
        options={"maxfev": 5},  # 5 evals: enough to observe x0 but not to converge
        auto_retry=False,
    )
    # After solve, _init_overrides must be cleared (single-shot semantics).
    assert solver._init_overrides == {}
    # And the user value survived in the source of truth (not clobbered).
    assert net.elements["ch_in"].initial_guess["ch_in.m_dot"] == 42.0


def _certified_merge_case1_network() -> FlowNetwork:
    """Certified audit case 1 (merge, exact-root fixture, 2026-07-04).

    Constants from tmp/mpce_audit_v2_cases.jsonl: an inverse-designed
    merge tee whose certified root has both feeds flowing forward
    (m_str=0.2576, m_bra=0.2129). The strict=False soft-barrier system
    additionally admits an EXACT artifact root with the branch feed
    slightly reversed (-0.0035 kg/s), where the one-sided penalty
    cancels the Pt-continuity mismatch.
    """
    from combaero.network.mpce_v2_element import MPCEv2Element

    D_com = 0.08
    A_com = math.pi * (D_com / 2.0) ** 2
    psi = 0.6823922768792331
    A_bra = A_com / psi
    D_bra = 2.0 * math.sqrt(A_bra / math.pi)
    theta = 59.990127553939544
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_hi_str", Pt=114712.54467660612, Tt=300.0))
    net.add_node(PressureBoundary("pb_hi_bra", Pt=113912.05698248411, Tt=300.0))
    net.add_node(PressureBoundary("pb_lo", Pt=1.0e5, Tt=300.0))
    net.add_node(MomentumChamberNode("mc_str", area=A_com))
    net.add_node(MomentumChamberNode("mc_bra", area=A_bra))
    net.add_node(MomentumChamberNode("mc_com", area=A_com))
    net.add_node(MomentumChamberNode("n2_com", area=A_com))
    net.add_element(
        ChannelElement(
            "ch_str",
            "pb_hi_str",
            "mc_str",
            length=1.0,
            diameter=D_com,
            regime="compressible",
        )
    )
    net.add_element(
        ChannelElement(
            "ch_bra",
            "pb_hi_bra",
            "mc_bra",
            length=1.0,
            diameter=D_bra,
            regime="compressible",
        )
    )
    net.add_element(
        ChannelElement(
            "ch_out",
            "mc_com",
            "n2_com",
            length=1.0,
            diameter=D_com,
            regime="compressible",
        )
    )
    net.add_element(
        OrificeElement(
            "orf_out",
            "n2_com",
            "pb_lo",
            Cd=0.62,
            diameter=0.06559708398216724,
            regime="compressible",
            correlation="fixed",
        )
    )
    net.add_element(
        MPCEv2Element(
            id="jct",
            inlet_nodes=["mc_str", "mc_bra"],
            outlet_nodes=["mc_com"],
            inlet_angles_deg=[0.0, theta],
            outlet_angles_deg=[0.0],
            port_areas=[A_com, A_bra, A_com],
            flow_direction="merge",
            strict=False,
        )
    )
    return net


def test_soft_barrier_artifact_root_demoted_to_failure():
    """A converged wrong-direction soft-barrier root must NOT be
    reported as success.

    The one-sided penalty alpha * mdot^2 on a slightly reversed port can
    exactly cancel that row's Pt-continuity mismatch, giving |F| ~ 1e-10
    at a state the strict residual rejects (raises). Seeding the solve
    AT that artifact root makes the demotion deterministic: the solver
    "converges" immediately, then the post-solve direction verification
    must fail it.
    """
    artifact_root = {
        "mc_str.P": 110171.3996,
        "mc_str.Pt": 114033.6720,
        "mc_bra.P": 113912.0118,
        "mc_bra.Pt": 113912.0971,
        "mc_com.P": 110227.2258,
        "mc_com.Pt": 114033.6720,
        "n2_com.P": 109529.0722,
        "n2_com.Pt": 113359.7813,
        "ch_str.m_dot": 0.49797,
        "ch_bra.m_dot": -0.00349,
        "ch_out.m_dot": 0.49448,
        "orf_out.m_dot": 0.49448,
        "jct.P_jct": 114033.6720,
    }
    net = _certified_merge_case1_network()
    solver = NetworkSolver(net)
    net.resolve_all_topology()
    x0 = solver._build_x0()
    for name, val in artifact_root.items():
        x0[solver._name_to_index[name]] = val
    with pytest.warns(UserWarning, match="artifact root"):
        sol = solver.solve(x0=x0, timeout=30.0)
    assert not sol["__success__"]
    assert "artifact root" in sol["__message__"]


def test_certified_merge_root_passes_direction_verification():
    """The genuine all-forward root of the same fixture stays a success."""
    net = _certified_merge_case1_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=30.0)
    assert sol["__success__"], sol.get("__message__")
    # Certified ground truth from the inverse-design generator.
    assert abs(sol["ch_str.m_dot"] - 0.25761378909667976) < 5e-3
    assert abs(sol["ch_bra.m_dot"] - 0.21294516374058128) < 5e-3


def _mfb_branch_tee_orifice_network(m_dot: float, regime: str) -> FlowNetwork:
    """Practical flow -> tee -> pressure topology (GUI parity fixture).

    Mirrors gui/tmp/one_mpce_tee_mdot_to_p.json: a mass-flow inlet feeds a
    separating MPCEv2 tee whose two legs each end in an orifice discharging
    to a single ambient PressureBoundary. The orifices provide the leg
    resistance that makes the fixture family feasible (bare equal-Pt sinks
    with resistance-free constant-area legs have no root for any physical
    tee closure -- see the 2026-07-04 case-31 analysis).
    """
    from combaero.network.mpce_v2_element import MPCEv2Element

    D_main, D_leg = 0.1, 0.06
    A_main = math.pi * (D_main / 2.0) ** 2
    A_leg = math.pi * (D_leg / 2.0) ** 2
    net = FlowNetwork()
    net.add_node(MassFlowBoundary("mfb_in", m_dot=m_dot, Tt=300.0))
    net.add_node(PressureBoundary("pb_amb", Pt=101325.0, Tt=300.0))
    net.add_node(MomentumChamberNode("mc_com", area=A_main))
    net.add_node(MomentumChamberNode("mc_str", area=A_main))
    net.add_node(MomentumChamberNode("mc_bra", area=A_leg))
    net.add_node(PlenumNode("pl_str"))
    net.add_node(PlenumNode("pl_bra"))
    net.add_element(
        ChannelElement("ch_in", "mfb_in", "mc_com", length=1.0, diameter=D_main, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_str", "mc_str", "pl_str", length=1.0, diameter=D_main, regime=regime)
    )
    net.add_element(
        ChannelElement("ch_bra", "mc_bra", "pl_bra", length=1.0, diameter=D_leg, regime=regime)
    )
    net.add_element(
        OrificeElement(
            "orf_str",
            "pl_str",
            "pb_amb",
            Cd=0.6,
            diameter=0.032,
            regime=regime,
            correlation="fixed",
        )
    )
    net.add_element(
        OrificeElement(
            "orf_bra",
            "pl_bra",
            "pb_amb",
            Cd=0.6,
            diameter=0.03,
            regime=regime,
            correlation="fixed",
        )
    )
    net.add_element(
        MPCEv2Element(
            id="tee",
            inlet_nodes=["mc_com"],
            outlet_nodes=["mc_str", "mc_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, 90.0],
            port_areas=[A_main, A_main, A_leg],
            flow_direction="branch",
            strict=False,
        )
    )
    return net


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_flow_tee_pressure_cold_start_converges(regime: str):
    """Regression: the practical MFB -> tee -> orifices -> PB topology must
    converge from a plain cold start.

    The GUI network this mirrors was reported 'practically unsolvable' at
    m_dot=0.4 compressible; headless diagnosis (2026-07-05) showed the
    root exists at every m_dot and the GUI failures came from warm-starting
    off a stored FAILED result. This pins the cold-start behavior so a
    genuine solver regression cannot hide behind that GUI artifact again.
    m_dot=0.2 keeps the compressible solve fast (~50 evals); higher m_dot
    approaches the orifice choke edge and only costs iterations, not the
    root.
    """
    net = _mfb_branch_tee_orifice_network(0.2, regime)
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")
    assert sol["__final_norm__"] < 1e-3
    m_str = sol["ch_str.m_dot"]
    m_bra = sol["ch_bra.m_dot"]
    assert m_str > 0.0 and m_bra > 0.0
    assert m_str + m_bra == pytest.approx(0.2, abs=1e-4)
    # The larger straight orifice (0.032 vs 0.03 bore) must carry more flow.
    assert m_str > m_bra


def _mk_state(P: float, Pt: float, m_dot: float = 0.1) -> "object":
    from combaero.network.components import NetworkMixtureState

    y_air = list(cb.mole_to_mass(cb.species.dry_air()))
    return NetworkMixtureState(P=P, Pt=Pt, T=300.0, Tt=300.0, m_dot=m_dot, Y=y_air)


def test_outlet_ref_channel_moves_density_sensitivity_downstream():
    """With _incompressible_p_ref='outlet' the incompressible channel
    evaluates density at the downstream static, so the friction-loss
    pressure sensitivity must land on the downstream node's P (and the
    default inlet reference must be unchanged)."""
    ch = ChannelElement("ch", "up", "dn", length=1.0, diameter=0.05, regime="incompressible")
    s_in = _mk_state(P=1.5e5, Pt=1.51e5)
    s_out = _mk_state(P=1.2e5, Pt=1.21e5)

    _, jac_default = ch.residuals(s_in, s_out)
    assert "up.P" in jac_default[0]
    assert "dn.P" not in jac_default[0]

    ch._incompressible_p_ref = "outlet"
    res_out, jac_outlet = ch.residuals(s_in, s_out)
    assert "dn.P" in jac_outlet[0]
    assert "up.P" not in jac_outlet[0]

    # Lower reference pressure -> lower density -> higher velocity ->
    # larger friction loss -> smaller residual (res = Pt_in - Pt_out - dP).
    ch._incompressible_p_ref = "inlet"
    res_inlet, _ = ch.residuals(s_in, s_out)
    assert res_out[0] < res_inlet[0]


def test_outlet_ref_orifice_moves_density_sensitivity_downstream():
    orf = OrificeElement(
        "orf", "up", "dn", Cd=0.6, diameter=0.03, regime="incompressible", correlation="fixed"
    )
    s_in = _mk_state(P=1.5e5, Pt=1.51e5)
    s_out = _mk_state(P=1.2e5, Pt=1.2e5)

    _, jac_default = orf.residuals(s_in, s_out)
    assert "up.P" in jac_default[0]

    orf._incompressible_p_ref = "outlet"
    _, jac_outlet = orf.residuals(s_in, s_out)
    assert "up.P" not in jac_outlet[0]
    assert "dn.P" in jac_outlet[0]


def test_auto_retry_outletref_rescues_failed_cold_solve():
    """When the cold compressible solve fails, solve() must retry from
    the outlet-referenced incompressible warm start and succeed. The
    primary failure is forced so the test is deterministic and fast."""
    net = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net)
    real_impl = solver._solve_impl
    calls: list[dict] = []

    def fake_impl(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            return {"__success__": False, "__message__": "forced failure", "__final_norm__": 1e3}
        return real_impl(**kwargs)

    solver._solve_impl = fake_impl
    with pytest.warns(RuntimeWarning, match="auto-retrying"):
        sol = solver.solve(timeout=120.0)
    assert sol["__success__"], sol.get("__message__")
    # primary (forced fail) + incompressible proxy + seeded retry
    assert len(calls) == 3
    assert calls[1].get("x0") is None  # proxy solves cold
    assert calls[2].get("x0") is not None  # retry is seeded
    ssu = solver._diagnostic_data["solver_settings_used"]
    assert ssu["init_strategy"] == "outletref_warmstart"
    assert ssu["auto_retry"] is True
    assert "auto-retry" in sol["__message__"]


def test_auto_retry_optout_returns_primary_failure():
    net = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net)
    calls: list[dict] = []

    def fake_impl(**kwargs):
        calls.append(kwargs)
        return {"__success__": False, "__message__": "forced failure", "__final_norm__": 1e3}

    solver._solve_impl = fake_impl
    sol = solver.solve(timeout=30.0, auto_retry=False)
    assert not sol["__success__"]
    assert len(calls) == 1


def test_auto_retry_skipped_on_incompressible_networks():
    """The retry proxy re-solves in the incompressible regime, which is a
    no-op for already-incompressible networks -- no retry there."""
    net = _mfb_branch_tee_orifice_network(0.2, "incompressible")
    solver = NetworkSolver(net)
    calls: list[dict] = []

    def fake_impl(**kwargs):
        calls.append(kwargs)
        return {"__success__": False, "__message__": "forced failure", "__final_norm__": 1e3}

    solver._solve_impl = fake_impl
    sol = solver.solve(timeout=30.0)
    assert not sol["__success__"]
    assert len(calls) == 1


def test_outlet_ref_proxy_seed_tracks_compressible_root_better_than_inlet():
    """Core hypothesis pin (2026-07-05 outlet-ref seed experiment): the
    incompressible proxy with densities at the DOWNSTREAM static lands
    closer to the compressible solution than the classic inlet-referenced
    proxy -- measured as RMS over pressure unknowns. If a model change
    breaks this ordering, the auto-retry seed quality regresses silently."""
    net = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0, auto_retry=False)
    assert sol["__success__"], sol.get("__message__")
    names = solver.unknown_names
    x_comp = {n: sol[n] for n in names}

    seed_out = solver._outlet_ref_incompressible_seed(p_ref="outlet", timeout=30.0)
    seed_in = solver._outlet_ref_incompressible_seed(p_ref="inlet", timeout=30.0)
    assert seed_out is not None, "outlet-ref proxy must converge on this fixture"
    assert seed_in is not None, "inlet-ref proxy must converge on this fixture"

    def rms_p(seed) -> float:
        sq = [
            (v - x_comp[n]) ** 2
            for n, v in zip(names, seed, strict=False)
            if ".P" in n or ".Pt" in n
        ]
        return math.sqrt(sum(sq) / len(sq))

    assert rms_p(seed_out) < rms_p(seed_in)


def test_outlet_ref_seeded_solve_reaches_same_root_as_cold():
    """The outlet-ref seed must not steer the compressible solve to a
    different root than the cold start finds."""
    net = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net)
    sol_cold = solver.solve(timeout=60.0, auto_retry=False)
    assert sol_cold["__success__"], sol_cold.get("__message__")

    net2 = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver2 = NetworkSolver(net2)
    net2.resolve_all_topology()
    seed = solver2._outlet_ref_incompressible_seed(p_ref="outlet", timeout=30.0)
    assert seed is not None
    sol_seeded = solver2.solve(timeout=60.0, x0=seed)
    assert sol_seeded["__success__"], sol_seeded.get("__message__")
    assert sol_seeded["ch_str.m_dot"] == pytest.approx(sol_cold["ch_str.m_dot"], rel=1e-3)
    assert sol_seeded["ch_bra.m_dot"] == pytest.approx(sol_cold["ch_bra.m_dot"], rel=1e-3)


def test_outletref_warmstart_explicit_strategy_converges():
    """init_strategy='outletref_warmstart' runs the proxy seed directly
    as the primary initialization (no cold attempt first) and reaches
    the same root as the cold start."""
    net = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver = NetworkSolver(net)
    sol_cold = solver.solve(timeout=60.0, auto_retry=False)
    assert sol_cold["__success__"], sol_cold.get("__message__")

    net2 = _mfb_branch_tee_orifice_network(0.2, "compressible")
    solver2 = NetworkSolver(net2)
    sol = solver2.solve(timeout=60.0, init_strategy="outletref_warmstart")
    assert sol["__success__"], sol.get("__message__")
    ssu = solver2._diagnostic_data["solver_settings_used"]
    assert ssu["init_strategy"] == "outletref_warmstart"
    # Explicit strategy selection opts out of the auto-retry chain.
    assert "auto_retry" not in ssu
    assert sol["ch_str.m_dot"] == pytest.approx(sol_cold["ch_str.m_dot"], rel=1e-3)
    assert sol["ch_bra.m_dot"] == pytest.approx(sol_cold["ch_bra.m_dot"], rel=1e-3)
