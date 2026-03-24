"""Test PipeElement htc_and_T() method."""

import math

import pytest

import combaero as cb
from combaero.network import (
    ConvectiveSurface,
    DimpledModel,
    MixtureState,
    PipeElement,
    RibbedModel,
    SmoothModel,
)

# Jacobian fields that every non-None htc_and_T() result must expose.
# These live on ChannelResult and are required for wall-coupling Jacobian relay.
CHANNEL_RESULT_JACOBIANS = {
    "dh_dmdot",
    "dh_dT",
    "ddP_dmdot",
    "ddP_dT",
    "dq_dmdot",
    "dq_dT",
    "dq_dT_wall",
}


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
        Y=cb.mole_to_mass(cb.species.dry_air()),
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
        Y=cb.mole_to_mass(cb.species.dry_air()),
    )

    # htc_and_T should compute values
    result = pipe.htc_and_T(state)
    assert result is not None

    assert result.h > 0  # HTC should be positive
    assert result.T_aw > 0  # Adiabatic wall temperature should be positive
    assert pipe.surface.area == 0.1  # Area should match surface area

    # T_aw should be close to T for low velocity
    assert abs(result.T_aw - state.T) < 10.0


def _make_state(m_dot: float = 0.5) -> MixtureState:
    """Helper: representative flowing state with non-trivial Re."""
    return MixtureState(
        T=600.0,
        P=200000.0,
        P_total=210000.0,
        T_total=610.0,
        m_dot=m_dot,
        Y=cb.mole_to_mass(cb.species.dry_air()),
    )


def _assert_jacobians_not_discarded(result: object, label: str) -> None:
    """Assert *result* exposes all ChannelResult Jacobian fields.

    Uses ``dir()`` introspection so the check is resilient to future
    return-type changes — any object that carries the fields passes.
    """
    if result is None:
        return  # None means "no surface"; nothing to discard
    attrs = set(dir(result))
    missing = CHANNEL_RESULT_JACOBIANS - attrs
    assert not missing, (
        f"{label} returns {type(result).__name__} missing Jacobians: {missing}. "
        f"Likely discarding ChannelResult and returning a plain tuple."
    )


# ---------------------------------------------------------------------------
# Parametrized: ConvectiveSurface models
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "model_name, model",
    [
        ("Smooth", SmoothModel()),
        ("Ribbed", RibbedModel(e_D=0.1, pitch_to_height=10.0, alpha_deg=90.0)),
        ("Dimpled", DimpledModel(d_Dh=0.2, h_d=0.2, S_d=2.0)),
    ],
    ids=["Smooth", "Ribbed", "Dimpled"],
)
def test_convective_surface_exposes_jacobians(model_name: str, model: object) -> None:
    """ConvectiveSurface.htc_and_T() must not discard ChannelResult Jacobians."""
    surface = ConvectiveSurface(area=0.15, model=model)
    state = _make_state()
    rho = state.density()
    pipe_area = math.pi * (0.04 / 2) ** 2
    velocity = state.m_dot / (rho * pipe_area)

    result = surface.htc_and_T(
        T=state.T,
        P=state.P,
        X=state.X,
        velocity=velocity,
        diameter=0.04,
        length=1.2,
        T_wall=math.nan,
    )
    _assert_jacobians_not_discarded(result, f"ConvectiveSurface({model_name})")


# ---------------------------------------------------------------------------
# Parametrized: PipeElement with different surface models
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "model_name, model",
    [
        ("Smooth", SmoothModel()),
        ("Ribbed", RibbedModel(e_D=0.1, pitch_to_height=10.0, alpha_deg=90.0)),
        ("Dimpled", DimpledModel(d_Dh=0.2, h_d=0.2, S_d=2.0)),
    ],
    ids=["Smooth", "Ribbed", "Dimpled"],
)
def test_pipe_element_exposes_jacobians(model_name: str, model: object) -> None:
    """PipeElement.htc_and_T() must not discard ChannelResult Jacobians."""
    surface = ConvectiveSurface(area=0.15, model=model)
    pipe = PipeElement(
        id=f"pipe_{model_name.lower()}",
        from_node="n_in",
        to_node="n_out",
        length=1.2,
        diameter=0.04,
        roughness=5e-5,
        surface=surface,
    )
    state = _make_state()
    result = pipe.htc_and_T(state)
    _assert_jacobians_not_discarded(result, f"PipeElement({model_name})")


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
        Y=cb.mole_to_mass(cb.species.dry_air()),
    )

    # htc_and_T should compute values with wall temperature
    result = pipe.htc_and_T(state)
    assert result is not None

    assert result.h > 0
    assert result.T_aw > 0
    assert pipe.surface.area == 0.1


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
        Y=cb.mole_to_mass(cb.species.dry_air()),
    )

    # Base class should return None
    result = element.htc_and_T(state)
    assert result is None
