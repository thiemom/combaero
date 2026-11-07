"""The frontend refuses an explicit null where the backend does (#1038).

``gui/frontend/src/utils/validation.ts`` lists, per node type, the fields
the backend schema types as a plain float with a default (``NULL_REFUSED``):
an absent key takes the default, an explicit null is refused ("Input should
be a valid number"). The client checks read an absent field as the default,
so they must refuse the null first, or a file the GUI calls valid fails at
Solve. This keeps the list true to the schema.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from pydantic import BaseModel, ValidationError

from gui.backend.schemas import DiscreteLossData, EjectorData, MPCETeeData, TurbomachineData

_VALIDATION = (
    Path(__file__).resolve().parents[2] / "gui" / "frontend" / "src" / "utils" / "validation.ts"
)

_SCHEMA_BY_TYPE: dict[str, type[BaseModel]] = {
    "compressor": TurbomachineData,
    "turbine": TurbomachineData,
    "mpce_tee": MPCETeeData,
    "ejector": EjectorData,
    "discrete_loss": DiscreteLossData,
}


def _null_refused() -> dict[str, list[str]]:
    text = _VALIDATION.read_text(encoding="utf-8")
    m = re.search(r"const NULL_REFUSED[^=]*=\s*\{(.*?)\n\};", text, re.S)
    assert m is not None, "NULL_REFUSED not found in validation.ts"
    out: dict[str, list[str]] = {}
    for node_type, inner in re.findall(r"(\w+):\s*\{([^}]*)\}", m.group(1)):
        out[node_type] = re.findall(r"(\w+):\s*\"", inner)
    return out


_CASES = [(t, f) for t, fields in _null_refused().items() for f in fields]


def test_the_list_is_parsed() -> None:
    assert ("compressor", "pressure_ratio") in _CASES
    assert {t for t, _ in _CASES} == set(_SCHEMA_BY_TYPE)


@pytest.mark.parametrize(("node_type", "field"), _CASES)
def test_backend_refuses_null_and_defaults_absent(node_type: str, field: str) -> None:
    schema = _SCHEMA_BY_TYPE[node_type]
    assert field in schema.model_fields, f"{schema.__name__} has no {field}"
    with pytest.raises(ValidationError):
        schema.model_validate({field: None})
    assert getattr(schema.model_validate({}), field) is not None
