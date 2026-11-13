"""Markers the suite uses.

`installed` marks the tests that check THIS machine rather than the code: they
need TNP, NanoBodyBuilder2, mkdssp 3.x and the Aggrescan3D pair, and they skip
where those are absent. Deselect them with `-m "not installed"`.
"""


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "installed: needs the structural toolchain on this machine")
