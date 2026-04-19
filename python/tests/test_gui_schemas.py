from gui.backend.schemas import ChannelData, MomentumChamberData, PlenumData


def test_channel_data_serialization():
    # Test data reflecting what the frontend sends
    test_data = {
        "L": 5.0,
        "D": 0.2,
        "roughness": 2e-5,
        "Nu_multiplier": 1.25,
        "f_multiplier": 0.9,
        "regime": "compressible",
        "initial_guess": {"m_dot": 1.5},
        "label": "Main Supply Channel",
    }

    channel = ChannelData(**test_data)
    assert channel.L == 5.0
    assert channel.D == 0.2
    assert channel.Nu_multiplier == 1.25
    assert channel.f_multiplier == 0.9
    assert channel.label == "Main Supply Channel"
    print("ChannelData validation successful")


def test_channel_data_multiplier_defaults():
    channel = ChannelData()
    assert channel.Nu_multiplier == 1.0
    assert channel.f_multiplier == 1.0


def test_momentum_chamber_multiplier_defaults():
    data = MomentumChamberData()
    assert data.Nu_multiplier == 1.0
    assert data.f_multiplier == 1.0


def test_momentum_chamber_multiplier_serialization():
    data = MomentumChamberData(area=0.2, Nu_multiplier=1.4, f_multiplier=0.85)
    assert data.area == 0.2
    assert data.Nu_multiplier == 1.4
    assert data.f_multiplier == 0.85


def test_plenum_data_label():
    data = PlenumData(label="Header")
    assert data.label == "Header"


if __name__ == "__main__":
    test_channel_data_serialization()
    test_plenum_data_label()
    print("All GUI schema tests passed!")
