import numpy as np
from ndgrid.grid import from_raw_data
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
