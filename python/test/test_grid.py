import numpy as np
from ndgrid.grid import from_raw_data
from ndgrid.shapes import regular_sphere
from ndelement.reference_cell import ReferenceCellType


def test_from_raw_data():
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )
    cells = np.array([[0, 1, 3], [1, 3, 2]])

    grid = from_raw_data(coordinates, cells, ReferenceCellType.Triangle, 1)

    assert grid.topology_dim == 2
    assert grid.geometry_dim == 3


def test_entities():
    grid = regular_sphere(4)
    for i in range(grid.entity_count(ReferenceCellType.Triangle)):
        entity = grid.entity(2, i)
        assert entity.local_index == i
