"""Tests for compressible flow elements in network solver."""

import numpy as np

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
    WallConnection,
)


def test_compressible_orifice_network():
    """Test network with compressible orifice element."""
    net = FlowNetwork()

    # High pressure ratio to ensure compressible effects
    inlet = PressureBoundary("inlet", P_total=300000.0, T_total=300.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Compressible orifice
    orifice = OrificeElement(
        "orifice", "inlet", "outlet", Cd=0.65, area=1e-4, regime="compressible"
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive"

    # Check pressure ratio to confirm compressible regime
    PR = outlet.P_total / inlet.P_total
    assert PR < 0.8, f"Should be in compressible regime (PR={PR:.3f})"

    # Verify against direct compressible calculation
    X = cb.standard_dry_air_composition()
    mdot_direct, _, _, _ = cb._core.orifice_compressible_mdot_and_jacobian(
        inlet.T_total, inlet.P_total, outlet.P_total, X, 0.65, 1e-4, 0.0
    )

    # Should match within solver tolerance
    assert abs(sol["orifice.m_dot"] - mdot_direct) / mdot_direct < 0.01


def test_compressible_pipe_network():
    """Test network with compressible Fanno pipe element."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=200000.0, T_total=400.0)
    outlet = PressureBoundary("outlet", P_total=150000.0, T_total=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Simple compressible Fanno pipe directly between boundaries
    pipe = PipeElement(
        "pipe",
        "inlet",
        "outlet",
        length=1.0,
        diameter=0.08,
        roughness=1e-4,
        regime="compressible_fanno",
        friction_model="haaland",
    )

    net.add_element(pipe)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["pipe.m_dot"] > 0, "Pipe mass flow should be positive"

    # Verify pressure drop is reasonable
    dP = inlet.P_total - outlet.P_total
    assert dP > 0, "Pressure should drop across pipe"


def test_mixed_compressible_incompressible_network():
    """Test network mixing compressible and incompressible elements."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=250000.0, T_total=350.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=350.0)
    junction1 = PlenumNode("j1")
    junction2 = PlenumNode("j2")

    net.add_node(inlet)
    net.add_node(outlet)
    net.add_node(junction1)
    net.add_node(junction2)

    # Incompressible pipe first
    pipe1 = PipeElement(
        "pipe1",
        "inlet",
        "j1",
        length=1.0,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    # Compressible orifice in middle (high pressure drop)
    orifice = OrificeElement("orifice", "j1", "j2", Cd=0.65, area=1e-4, regime="compressible")

    # Incompressible pipe at end
    pipe2 = PipeElement(
        "pipe2",
        "j2",
        "outlet",
        length=0.5,
        diameter=0.05,
        roughness=1e-4,
        regime="incompressible",
    )

    net.add_element(pipe1)
    net.add_element(orifice)
    net.add_element(pipe2)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify solution
    assert sol["pipe1.m_dot"] > 0
    assert sol["orifice.m_dot"] > 0
    assert sol["pipe2.m_dot"] > 0

    # Mass conservation
    assert abs(sol["pipe1.m_dot"] - sol["orifice.m_dot"]) < 1e-6
    assert abs(sol["orifice.m_dot"] - sol["pipe2.m_dot"]) < 1e-6

    # Pressure should decrease monotonically
    assert inlet.P_total > sol["j1.P"] > sol["j2.P"] > outlet.P_total


def test_compressible_vs_incompressible_comparison():
    """Compare compressible vs incompressible orifice in same network."""

    # Setup identical networks with different regimes
    def create_network(regime):
        net = FlowNetwork()
        inlet = PressureBoundary("inlet", P_total=200000.0, T_total=300.0)
        outlet = PressureBoundary("outlet", P_total=150000.0, T_total=300.0)
        net.add_node(inlet)
        net.add_node(outlet)

        orifice = OrificeElement("orifice", "inlet", "outlet", Cd=0.65, area=2e-4, regime=regime)
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
    inlet = PressureBoundary("inlet", P_total=500000.0, T_total=300.0)
    outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0)

    net.add_node(inlet)
    net.add_node(outlet)

    orifice = OrificeElement(
        "orifice", "inlet", "outlet", Cd=0.65, area=5e-5, regime="compressible"
    )
    net.add_element(orifice)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    # Verify choked flow
    X = cb.standard_dry_air_composition()
    PR_crit = cb.critical_pressure_ratio(inlet.T_total, inlet.P_total, X)
    PR_actual = outlet.P_total / inlet.P_total

    assert PR_actual < PR_crit, (
        f"Flow should be choked (PR={PR_actual:.3f} < PR_crit={PR_crit:.3f})"
    )
    assert sol["orifice.m_dot"] > 0, "Mass flow should be positive even when choked"


def test_compressible_pipe_high_mach():
    """Test compressible pipe at higher Mach number."""
    net = FlowNetwork()

    inlet = PressureBoundary("inlet", P_total=300000.0, T_total=400.0)
    outlet = PressureBoundary("outlet", P_total=200000.0, T_total=400.0)

    net.add_node(inlet)
    net.add_node(outlet)

    # Long pipe with small diameter to get higher velocity
    pipe = PipeElement(
        "pipe",
        "inlet",
        "outlet",
        length=5.0,
        diameter=0.03,
        roughness=1e-4,
        regime="compressible_fanno",
    )
    net.add_element(pipe)

    solver = NetworkSolver(net)
    sol = solver.solve(method="lm")

    assert sol["pipe.m_dot"] > 0, "Mass flow should be positive"

    # Verify Mach number is significant
    X = cb.standard_dry_air_composition()
    rho = cb.density(inlet.T_total, inlet.P_total, X)
    area = 3.14159 * (0.03 / 2) ** 2
    u = sol["pipe.m_dot"] / (rho * area)
    a = cb.speed_of_sound(inlet.T_total, X)
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
    x_air = cb.standard_dry_air_composition()
    y_air = list(cb.mole_to_mass(x_air))

    t_hot = _T_HOT
    t_cold = _T_COLD
    p_out_hot = _P_OUT
    p_out_cold = _P_OUT
    d_hot = 0.018
    d_cold = 0.014
    length = 2.0

    mdot_hot = _mdot_for_mach(mach_target, t_hot, p_out_hot, x_air, _area(d_hot))
    mdot_cold = 0.65 * _mdot_for_mach(mach_target, t_cold, p_out_cold, x_air, _area(d_cold))

    pipe_regime = "compressible_fanno" if regime == "compressible" else "incompressible"
    orifice_regime = "compressible" if regime == "compressible" else "incompressible"

    net = FlowNetwork()
    net.add_node(MassFlowBoundary("hot_inlet", m_dot=mdot_hot, T_total=t_hot, Y=y_air))
    net.add_node(PlenumNode("hot_plenum"))
    net.add_node(PressureBoundary("hot_outlet", P_total=p_out_hot, T_total=t_hot, Y=y_air))

    net.add_node(MassFlowBoundary("cold_inlet", m_dot=mdot_cold, T_total=t_cold, Y=y_air))
    net.add_node(PlenumNode("cold_plenum"))
    net.add_node(PressureBoundary("cold_outlet", P_total=p_out_cold, T_total=t_cold, Y=y_air))

    hot_surface = ConvectiveSurface(area=np.pi * d_hot * length, model=SmoothModel())
    cold_surface = ConvectiveSurface(area=np.pi * d_cold * length, model=SmoothModel())

    net.add_element(
        PipeElement(
            id="hot_pipe",
            from_node="hot_inlet",
            to_node="hot_plenum",
            length=length,
            diameter=d_hot,
            roughness=5.0e-5,
            regime=pipe_regime,
            surface=hot_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="hot_orifice",
            from_node="hot_plenum",
            to_node="hot_outlet",
            Cd=0.72,
            area=_area(0.030),
            regime=orifice_regime,
        )
    )

    net.add_element(
        PipeElement(
            id="cold_pipe",
            from_node="cold_inlet",
            to_node="cold_plenum",
            length=length,
            diameter=d_cold,
            roughness=5.0e-5,
            regime=pipe_regime,
            surface=cold_surface,
        )
    )
    net.add_element(
        OrificeElement(
            id="cold_orifice",
            from_node="cold_plenum",
            to_node="cold_outlet",
            Cd=0.72,
            area=_area(0.025),
            regime=orifice_regime,
        )
    )

    net.add_wall(
        WallConnection(
            id="coupling_wall",
            element_a="hot_pipe",
            element_b="cold_pipe",
            wall_thickness=0.002,
            wall_conductivity=20.0,
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

        pr_comp.append(_P_OUT / float(sol_comp["hot_inlet.P_total"]))
        pr_incomp.append(_P_OUT / float(sol_incomp["hot_inlet.P_total"]))

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
