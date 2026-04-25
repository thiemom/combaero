from combaero.network.components import (
    ChannelElement,
    LosslessConnectionElement,
    PlenumNode,
    PressureBoundary,
)
from combaero.network.graph import FlowNetwork


def test_flownetwork_serialization():
    # Build a sample network
    network1 = FlowNetwork()

    # Nodes
    n1 = PressureBoundary(" inlet ", Pt=150000.0, Tt=300.0, Y=[0.0] * 14)
    n2 = PlenumNode(" plenum1")
    n3 = PressureBoundary(" outlet ", Pt=100000.0, Tt=300.0)

    network1.add_node(n1)
    network1.add_node(n2)
    network1.add_node(n3)

    # Elements
    e1 = LosslessConnectionElement(" inlet_conn ", " inlet ", " plenum1")
    e2 = ChannelElement(
        " channel1",
        " plenum1",
        " outlet ",
        length=2.5,
        diameter=0.05,
        roughness=1e-5,
        regime="incompressible",
        friction_model="haaland",
    )

    network1.add_element(e1)
    network1.add_element(e2)

    # Serialize
    data = network1.to_dict()

    # Verify basic structure
    assert "nodes" in data
    assert "elements" in data
    assert " inlet " in data["nodes"]
    assert " channel1" in data["elements"]

    # Verify registry_id is captured
    channel_kwargs = data["elements"][" channel1"]["kwargs"]
    assert channel_kwargs["friction_model"] == "haaland"
    assert channel_kwargs["length"] == 2.5

    # Deserialize
    network2 = FlowNetwork.from_dict(data)

    # Validate reconstructed network
    assert len(network2.nodes) == 3
    assert len(network2.elements) == 2

    channel_reconstructed = network2.elements[" channel1"]
    assert isinstance(channel_reconstructed, ChannelElement)
    assert channel_reconstructed.friction_model == "haaland"
    assert channel_reconstructed.length == 2.5
