from pathlib import Path

import pytest

from crankshaft_balancing import load_shaft

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "example_shaft.json"


@pytest.fixture(scope="session")
def shaft():
    return load_shaft(EXAMPLE)
