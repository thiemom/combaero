"""Tests for EnergyBoundary integration with network nodes."""

import pytest

import combaero as cb
from combaero.network import (
    CombustorNode,
    EnergyBoundary,
    FlowNetwork,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _air_Y() -> list[float]:
    ns = cb.num_species()
    Y = [0.0] * ns
    for k in range(ns):
        if cb.species_name(k) == "N2":
            Y[k] = 0.767
        elif cb.species_name(k) == "O2":
            Y[k] = 0.233
    return Y


def _ch4_Y() -> list[float]:
    ns = cb.num_species()
    Y = [0.0] * ns
    for k in range(ns):
        if cb.species_name(k) == "CH4":
            Y[k] = 1.0
    return Y


# =============================================================================
# 1. EnergyBoundary class basics
# =============================================================================


class TestEnergyBoundaryClass:
    def test_default_Q_zero(self) -> None:
        eb = EnergyBoundary("eb1")
        assert eb.Q == 0.0
        assert eb.id == "eb1"

    def test_positive_Q(self) -> None:
        eb = EnergyBoundary("heater", Q=5000.0)
        assert eb.Q == 5000.0

    def test_negative_Q(self) -> None:
        eb = EnergyBoundary("cooler", Q=-3000.0)
        assert eb.Q == -3000.0

    def test_mutable_Q(self) -> None:
        eb = EnergyBoundary("eb", Q=1000.0)
        eb.Q = -2000.0
        assert eb.Q == -2000.0


# =============================================================================
# 2. PlenumNode with EnergyBoundary
# =============================================================================


class TestPlenumEnergyBoundary:
    def test_plenum_no_energy_boundary_unchanged(self) -> None:
        """PlenumNode with no energy boundaries behaves as before."""
        node = PlenumNode("p1")
        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=300.0,
            T_total=300.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        T, Y, _ = node.compute_derived_state([up])
        assert pytest.approx(300.0, abs=0.1) == T

    def test_plenum_heating_raises_temperature(self) -> None:
        """Positive Q raises the mixed temperature."""
        node = PlenumNode("p1")
        eb = EnergyBoundary("heater", Q=50000.0)  # 50 kW
        node.add_energy_boundary(eb)

        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=300.0,
            T_total=300.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        T_heated, _, _ = node.compute_derived_state([up])

        # Without energy boundary
        node2 = PlenumNode("p2")
        T_base, _, _ = node2.compute_derived_state([up])

        assert T_heated > T_base
        # delta_h = 50000 W / 1 kg/s = 50000 J/kg, ΔT ≈ 50 K for air
        assert T_heated - T_base == pytest.approx(50000.0 / 1005.0, rel=0.05)

    def test_plenum_cooling_lowers_temperature(self) -> None:
        """Negative Q lowers the mixed temperature."""
        node = PlenumNode("p1")
        eb = EnergyBoundary("cooler", Q=-80000.0)  # -80 kW
        node.add_energy_boundary(eb)

        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=600.0,
            T_total=600.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        T_cooled, _, _ = node.compute_derived_state([up])

        node2 = PlenumNode("p2")
        T_base, _, _ = node2.compute_derived_state([up])

        assert T_cooled < T_base

    def test_plenum_multiple_energy_boundaries(self) -> None:
        """Multiple boundaries sum their Q values."""
        node = PlenumNode("p1")
        node.add_energy_boundary(EnergyBoundary("h1", Q=30000.0))
        node.add_energy_boundary(EnergyBoundary("h2", Q=20000.0))

        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=300.0,
            T_total=300.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        T_two, _, _ = node.compute_derived_state([up])

        # Compare with single 50 kW
        node2 = PlenumNode("p2")
        node2.add_energy_boundary(EnergyBoundary("h_total", Q=50000.0))
        T_one, _, _ = node2.compute_derived_state([up])

        assert T_two == pytest.approx(T_one, abs=0.01)

    def test_plenum_enthalpy_conservation(self) -> None:
        """h_out = h_in + Q/mdot."""
        Q = 40000.0
        mdot = 1.5

        node = PlenumNode("p1")
        node.add_energy_boundary(EnergyBoundary("eb", Q=Q))

        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=400.0,
            T_total=400.0,
            m_dot=mdot,
            Y=_air_Y(),
        )
        T_out, Y_out, _ = node.compute_derived_state([up])

        X_out = cb.mass_to_mole(cb.normalize_fractions(Y_out))
        X_in = cb.mass_to_mole(cb.normalize_fractions(up.Y))

        h_in = cb.h_mass(up.T_total, X_in)
        h_out = cb.h_mass(T_out, X_out)
        delta_h = Q / mdot

        assert h_out == pytest.approx(h_in + delta_h, rel=1e-6)

    def test_plenum_zero_Q_noop(self) -> None:
        """Zero-Q boundary doesn't change temperature."""
        node = PlenumNode("p1")
        node.add_energy_boundary(EnergyBoundary("eb", Q=0.0))

        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=500.0,
            T_total=500.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        T_with, _, _ = node.compute_derived_state([up])

        node2 = PlenumNode("p2")
        T_without, _, _ = node2.compute_derived_state([up])

        assert T_with == pytest.approx(T_without, abs=1e-10)


# =============================================================================
# 3. CombustorNode with EnergyBoundary
# =============================================================================


class TestCombustorEnergyBoundary:
    def _make_combustor_upstream(self):
        """Returns (fuel_state, air_state) as MixtureState objects."""
        from combaero.network import MixtureState

        fuel = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=300.0,
            T_total=300.0,
            m_dot=0.05,
            Y=_ch4_Y(),
        )
        air = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=600.0,
            T_total=600.0,
            m_dot=1.0,
            Y=_air_Y(),
        )
        return fuel, air

    def test_combustor_cooling_lowers_T(self) -> None:
        fuel, air = self._make_combustor_upstream()

        node_base = CombustorNode("comb_base", method="complete")
        T_base, _, _ = node_base.compute_derived_state([fuel, air])

        node_cooled = CombustorNode("comb_cooled", method="complete")
        node_cooled.add_energy_boundary(EnergyBoundary("cool", Q=-200000.0))
        T_cooled, _, _ = node_cooled.compute_derived_state([fuel, air])

        assert T_cooled < T_base

    def test_combustor_zero_Q_noop(self) -> None:
        fuel, air = self._make_combustor_upstream()

        node_zero = CombustorNode("comb_zero", method="complete")
        node_zero.add_energy_boundary(EnergyBoundary("eb", Q=0.0))
        T_zero, _, _ = node_zero.compute_derived_state([fuel, air])

        node_none = CombustorNode("comb_none", method="complete")
        T_none, _, _ = node_none.compute_derived_state([fuel, air])

        assert T_zero == pytest.approx(T_none, abs=1e-6)

    def test_combustor_equilibrium_cooling(self) -> None:
        fuel, air = self._make_combustor_upstream()

        node_base = CombustorNode("comb_base", method="equilibrium")
        T_base, _, _ = node_base.compute_derived_state([fuel, air])

        node_cooled = CombustorNode("comb_cooled", method="equilibrium")
        node_cooled.add_energy_boundary(EnergyBoundary("cool", Q=-100000.0))
        T_cooled, _, _ = node_cooled.compute_derived_state([fuel, air])

        assert T_cooled < T_base


# =============================================================================
# 4. MomentumChamberNode with EnergyBoundary
# =============================================================================


class TestMomentumChamberEnergyBoundary:
    def test_momentum_chamber_heating(self) -> None:
        node = MomentumChamberNode("mc1")
        node.add_energy_boundary(EnergyBoundary("heater", Q=60000.0))

        from combaero.network import MixtureState

        up = MixtureState(
            P=200000.0,
            P_total=200000.0,
            T=400.0,
            T_total=400.0,
            m_dot=2.0,
            Y=_air_Y(),
        )
        T_heated, _, _ = node.compute_derived_state([up])

        node2 = MomentumChamberNode("mc2")
        T_base, _, _ = node2.compute_derived_state([up])

        assert T_heated > T_base
        # delta_h = 60000/2 = 30000 J/kg, ΔT ≈ 30 K
        assert T_heated - T_base == pytest.approx(30000.0 / 1014.0, rel=0.05)


# =============================================================================
# 5. Sign convention verification
# =============================================================================


class TestSignConvention:
    def test_positive_Q_heats(self) -> None:
        """Positive Q always increases temperature."""
        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=300.0,
            T_total=300.0,
            m_dot=1.0,
            Y=_air_Y(),
        )

        node_base = PlenumNode("base")
        T_base, _, _ = node_base.compute_derived_state([up])

        for Q in [1000.0, 10000.0, 100000.0]:
            node = PlenumNode(f"h_{Q}")
            node.add_energy_boundary(EnergyBoundary("eb", Q=Q))
            T, _, _ = node.compute_derived_state([up])
            assert T_base < T, f"Q={Q} should heat, but T={T} <= T_base={T_base}"

    def test_negative_Q_cools(self) -> None:
        """Negative Q always decreases temperature."""
        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=800.0,
            T_total=800.0,
            m_dot=1.0,
            Y=_air_Y(),
        )

        node_base = PlenumNode("base")
        T_base, _, _ = node_base.compute_derived_state([up])

        for Q in [-1000.0, -10000.0, -100000.0]:
            node = PlenumNode(f"c_{Q}")
            node.add_energy_boundary(EnergyBoundary("eb", Q=Q))
            T, _, _ = node.compute_derived_state([up])
            assert T_base > T, f"Q={Q} should cool, but T={T} >= T_base={T_base}"

    def test_opposite_Q_cancel(self) -> None:
        """Equal and opposite Q boundaries cancel out."""
        from combaero.network import MixtureState

        up = MixtureState(
            P=101325.0,
            P_total=101325.0,
            T=500.0,
            T_total=500.0,
            m_dot=1.0,
            Y=_air_Y(),
        )

        node = PlenumNode("p")
        node.add_energy_boundary(EnergyBoundary("h", Q=50000.0))
        node.add_energy_boundary(EnergyBoundary("c", Q=-50000.0))
        T, _, _ = node.compute_derived_state([up])

        node2 = PlenumNode("p2")
        T_base, _, _ = node2.compute_derived_state([up])

        assert pytest.approx(T_base, abs=0.01) == T


# =============================================================================
# 6. Full network solve with energy boundary
# =============================================================================


class TestNetworkWithEnergyBoundary:
    def test_simple_network_with_heated_plenum(self) -> None:
        """PressureBoundary -> Orifice -> heated Plenum -> Orifice -> PressureBoundary.
        The heated plenum should have a higher temperature than without heating."""
        graph = FlowNetwork()

        inlet = PressureBoundary("inlet", P_total=200000.0, T_total=300.0, Y=_air_Y())
        outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=_air_Y())

        plenum = PlenumNode("plenum")
        plenum.add_energy_boundary(EnergyBoundary("heater", Q=50000.0))

        orf1 = OrificeElement("orf1", "inlet", "plenum", Cd=0.8, area=0.01)
        orf2 = OrificeElement("orf2", "plenum", "outlet", Cd=0.8, area=0.01)

        graph.add_node(inlet)
        graph.add_node(plenum)
        graph.add_node(outlet)
        graph.add_element(orf1)
        graph.add_element(orf2)

        solver = NetworkSolver(graph)
        solution = solver.solve()

        # Should converge
        assert solution is not None
        mdot = solution["orf1.m_dot"]
        assert mdot > 0

        # Plenum T should be above inlet T
        T_plenum = solver._derived_states["plenum"][0]
        assert T_plenum > 300.0

        # Verify: delta_h = Q / mdot
        expected_delta_h = 50000.0 / mdot
        X_air = cb.mass_to_mole(cb.normalize_fractions(_air_Y()))
        h_in = cb.h_mass(300.0, X_air)
        h_out = cb.h_mass(T_plenum, X_air)
        assert h_out == pytest.approx(h_in + expected_delta_h, rel=1e-4)

    def test_cooled_plenum_network(self) -> None:
        """Network with a cooled plenum: temperature should drop."""
        graph = FlowNetwork()

        inlet = PressureBoundary("inlet", P_total=300000.0, T_total=800.0, Y=_air_Y())
        outlet = PressureBoundary("outlet", P_total=101325.0, T_total=300.0, Y=_air_Y())

        plenum = PlenumNode("plenum")
        # Use moderate Q relative to expected mdot (~10 kg/s with large orifice)
        plenum.add_energy_boundary(EnergyBoundary("cooler", Q=-50000.0))

        orf1 = OrificeElement("orf1", "inlet", "plenum", Cd=0.8, area=0.02)
        orf2 = OrificeElement("orf2", "plenum", "outlet", Cd=0.8, area=0.02)

        graph.add_node(inlet)
        graph.add_node(plenum)
        graph.add_node(outlet)
        graph.add_element(orf1)
        graph.add_element(orf2)

        solver = NetworkSolver(graph)
        solution = solver.solve()

        assert solution is not None
        T_plenum = solver._derived_states["plenum"][0]
        assert T_plenum < 800.0

        # Verify: delta_h = Q / mdot
        mdot = solution["orf1.m_dot"]
        expected_delta_h = -50000.0 / mdot
        X_air = cb.mass_to_mole(cb.normalize_fractions(_air_Y()))
        h_in = cb.h_mass(800.0, X_air)
        h_out = cb.h_mass(T_plenum, X_air)
        assert h_out == pytest.approx(h_in + expected_delta_h, rel=1e-4)


def test_heat_exchange_node():
    """Verify Q parameter changes outlet temperature correctly."""
    graph = FlowNetwork()

    # Simple network: MassFlowBoundary → Plenum(Q=1000W) → PressureBoundary
    Y_air = _air_Y()
    inlet = MassFlowBoundary("inlet", m_dot=0.1, T_total=600.0, Y=Y_air)
    plenum = PlenumNode("plenum")
    outlet = PressureBoundary("outlet", P_total=101325.0)

    # Add heat exchange: 1000W heating
    heater = EnergyBoundary("heater", Q=1000.0)
    plenum.add_energy_boundary(heater)

    orf1 = OrificeElement("orf1", "inlet", "plenum", Cd=0.6, area=0.001)
    orf2 = OrificeElement("orf2", "plenum", "outlet", Cd=0.6, area=0.001)

    graph.add_node(inlet)
    graph.add_node(plenum)
    graph.add_node(outlet)
    graph.add_element(orf1)
    graph.add_element(orf2)

    solver = NetworkSolver(graph)
    sol = solver.solve()

    assert sol["__success__"]

    # Debug: Check if EnergyBoundary is recognized
    print(f"Number of energy boundaries on plenum: {len(plenum.energy_boundaries)}")
    print(f"Heater Q: {plenum.energy_boundaries[0].Q}")

    # Get inlet and outlet states
    T_in = 600.0
    mdot = 0.1
    Q = 1000.0

    # Get actual outlet temperature from derived state
    T_out = solver._derived_states["plenum"][0]

    # Verify exact energy balance: h_out = h_in + Q/mdot
    X_air = cb.mass_to_mole(cb.normalize_fractions(Y_air))
    h_in = cb.h_mass(T_in, X_air)
    h_out_actual = cb.h_mass(T_out, X_air)

    # Expected enthalpy after heat addition
    expected_h_out = h_in + Q / mdot

    # Verify enthalpy balance (this is the fundamental conservation law)
    assert h_out_actual == pytest.approx(expected_h_out, rel=1e-3), (
        f"Energy balance failed: h_out={h_out_actual:.1f} J/kg, expected={expected_h_out:.1f} J/kg"
    )

    # Also verify temperature increased (sanity check)
    delta_T = T_out - T_in
    assert delta_T > 0, f"Temperature should increase: ΔT = {delta_T:.1f}K"

    print(f"Heat exchange test passed: T_in={T_in:.1f}K → T_out={T_out:.1f}K (ΔT={delta_T:.1f}K)")
    print(
        f"Energy balance verified: h_in={h_in:.1f} J/kg + Q/mdot={Q / mdot:.1f} J/kg = h_out={h_out_actual:.1f} J/kg"
    )

    # Test that Q parameter affects the solution by comparing with Q=0
    # Create reference network with no heat exchange
    ref_graph = FlowNetwork()
    ref_inlet = MassFlowBoundary("inlet", m_dot=0.1, T_total=600.0, Y=Y_air)
    ref_plenum = PlenumNode("plenum")
    ref_outlet = PressureBoundary("outlet", P_total=101325.0)

    ref_orf1 = OrificeElement("orf1", "inlet", "plenum", Cd=0.6, area=0.001)
    ref_orf2 = OrificeElement("orf2", "plenum", "outlet", Cd=0.6, area=0.001)

    ref_graph.add_node(ref_inlet)
    ref_graph.add_node(ref_plenum)
    ref_graph.add_node(ref_outlet)
    ref_graph.add_element(ref_orf1)
    ref_graph.add_element(ref_orf2)

    ref_solver = NetworkSolver(ref_graph)
    ref_solver.solve()
    T_ref = ref_solver._derived_states["plenum"][0]

    # Heated plenum should be hotter than reference
    assert T_out > T_ref, (
        f"Heated plenum ({T_out:.1f}K) should be hotter than reference ({T_ref:.1f}K)"
    )

    print(f"Reference test passed: Reference T={T_ref:.1f}K, Heated T={T_out:.1f}K")
