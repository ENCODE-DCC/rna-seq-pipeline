from pathlib import Path

import pytest


@pytest.fixture
def test_data_dir():
    return Path("test_data")
