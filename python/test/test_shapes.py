import pytest
from ndgrid.shapes import regular_sphere


@pytest.mark.parametrize("level", range(0, 4))
def test_regular_sphere(level):
    grid = regular_sphere(level)

    assert grid.topology_dim == 2
    assert grid.geometry_dim == 3
