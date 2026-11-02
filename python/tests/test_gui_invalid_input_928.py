"""#928: invalid input reaches the GUI user as a 400 with its message.

The backend maps ValueError to 400 and every other exception to 500. Before
#928 the humid-air functions refused a relative humidity outside [0, 1] with a
RuntimeError (C++ std::runtime_error), so a typo in a boundary's humidity came
back as a server error. Invalid input now raises ValueError everywhere.
"""

from __future__ import annotations

import pytest
import test_gui_http_routes as http


def _humid_pipe(rh: float) -> dict:
    g = http._pipe()
    g["nodes"][0]["data"]["composition"] = {
        "source": "humid_air",
        "mode": "mole",
        "relative_humidity": rh,
    }
    return g


@pytest.mark.parametrize("path", ["/solve", "/resolve_preview"])
def test_a_relative_humidity_above_one_is_a_400_naming_it(path: str) -> None:
    r = http._post(path, _humid_pipe(1.5))
    http._assert_readable_error(r, "Relative humidity must be between 0 and 1")


def test_the_global_ambient_humidity_is_refused_the_same_way() -> None:
    g = _humid_pipe(0.5)
    del g["nodes"][0]["data"]["composition"]["relative_humidity"]
    g["solver_settings"] = {"ambient_RH": -0.1}
    r = http._post("/solve", g)
    http._assert_readable_error(r, "Relative humidity must be between 0 and 1")


def test_a_valid_humidity_still_solves() -> None:
    body = http._post("/solve", _humid_pipe(0.6)).json()
    assert body["success"], body["message"]
