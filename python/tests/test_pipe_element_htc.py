"""Test PipeElement htc_and_T() method."""

import combaero as cb
from combaero.network import ConvectiveSurface, MixtureState, PipeElement, SmoothModel


def test_pipe_element_default_surface():
    """Test PipeElement with default surface (area=0)."""
    pipe = PipeElement(
        id="pipe1",
        from_node="node1",
        to_node="node2",
        length=1.0,
        diameter=0.05,
        roughness=1e-5,
    )

    # Default surface has area=0
    assert pipe.surface.area == 0.0

    # Create test state
    state = MixtureState(
        T=300.0,
        P=101325.0,
        P_total=101325.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.mole_to_mass(cb.standard_dry_air_composition()),
    )

    # htc_and_T should return None for area=0
    result = pipe.htc_and_T(state)
    assert result is None


def test_pipe_element_with_surface():
    """Test PipeElement with custom ConvectiveSurface."""
    surface = ConvectiveSurface(
        area=0.1,
        model=SmoothModel(correlation="gnielinski"),
    )

    pipe = PipeElement(
        id="pipe2",
        from_node="node1",
        to_node="node2",
        length=1.0,
        diameter=0.05,
        roughness=1e-5,
        surface=surface,
    )

    # Surface should be set
    assert pipe.surface.area == 0.1
    assert pipe.surface.model.correlation == "gnielinski"

    # Create test state
    state = MixtureState(
        T=300.0,
        P=101325.0,
        P_total=101325.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.mole_to_mass(cb.standard_dry_air_composition()),
    )

    # htc_and_T should compute values
    result = pipe.htc_and_T(state)
    assert result is not None

    h, T_aw, A = result
    assert h > 0  # HTC should be positive
    assert T_aw > 0  # Adiabatic wall temperature should be positive
    assert A == 0.1  # Area should match surface area

    # T_aw should be close to T for low velocity
    assert abs(T_aw - state.T) < 10.0


def test_pipe_element_with_t_wall():
    """Test PipeElement with specified wall temperature."""
    surface = ConvectiveSurface(
        area=0.1,
        model=SmoothModel(correlation="gnielinski"),
    )

    pipe = PipeElement(
        id="pipe3",
        from_node="node1",
        to_node="node2",
        length=1.0,
        diameter=0.05,
        roughness=1e-5,
        surface=surface,
        t_wall=400.0,  # Wall temperature specified
    )

    # Create test state
    state = MixtureState(
        T=300.0,
        P=101325.0,
        P_total=101325.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.mole_to_mass(cb.standard_dry_air_composition()),
    )

    # htc_and_T should compute values with wall temperature
    result = pipe.htc_and_T(state)
    assert result is not None

    h, T_aw, A = result
    assert h > 0
    assert T_aw > 0
    assert A == 0.1


def test_network_element_default():
    """Test NetworkElement base class default implementation."""
    from combaero.network.components import NetworkElement

    class TestElement(NetworkElement):
        def unknowns(self):
            return []

        def residuals(self, state_in, state_out):
            return [], {}

        def n_equations(self):
            return 0

        def resolve_topology(self, graph):
            pass

    element = TestElement("test", "node1", "node2")

    # Create test state
    state = MixtureState(
        T=300.0,
        P=101325.0,
        P_total=101325.0,
        T_total=300.0,
        m_dot=0.1,
        Y=cb.mole_to_mass(cb.standard_dry_air_composition()),
    )

    # Base class should return None
    result = element.htc_and_T(state)
    assert result is None
