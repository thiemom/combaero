from gui.backend.schemas import PipeData, PlenumData


def test_pipe_data_serialization():
    # Test data reflecting what the frontend sends
    test_data = {
        "L": 5.0,
        "D": 0.2,
        "roughness": 2e-5,
        "regime": "compressible",
        "initial_guess": {"m_dot": 1.5},
        "label": "Main Supply Pipe",
    }

    pipe = PipeData(**test_data)
    assert pipe.L == 5.0
    assert pipe.D == 0.2
    assert pipe.label == "Main Supply Pipe"
    print("PipeData validation successful")


def test_plenum_data_label():
    data = PlenumData(label="Header")
    assert data.label == "Header"


if __name__ == "__main__":
    test_pipe_data_serialization()
    test_plenum_data_label()
    print("All GUI schema tests passed!")
