import pytest
import numpy as np
from ndgrid.shapes import regular_sphere
from ndelement.reference_cell import ReferenceCellType


@pytest.mark.parametrize("level", range(4))
def test_points(level):
    grid = regular_sphere(level)
    pts = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ]
    )
    gmap = grid.geometry_map(ReferenceCellType.Triangle, pts)

    for i in range(grid.entity_count(ReferenceCellType.Triangle)):
        mapped_pts = np.array([[10.0 for j in range(3)] for i in range(3)])
        gmap.points(i, mapped_pts)
        for p in mapped_pts:
            assert np.isclose(p.dot(p), 1.0)


@pytest.mark.parametrize("level", range(4))
def test_jacobians_dets_normals_cell(level):
    grid = regular_sphere(level)
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1 / 3, 1 / 3]])
    gmap = grid.geometry_map(ReferenceCellType.Triangle, pts)

    assert gmap.geometry_dimension == 3
    assert gmap.point_count == 4
    assert gmap.entity_topology_dimension == 2

    mapped_pts = np.empty((4, 3))
    jacobians = np.empty((4, 2, 3))
    jacobians2 = np.empty((4, 2, 3))
    jdets = np.empty(4)
    normals = np.empty((4, 3))
    for i in range(grid.entity_count(ReferenceCellType.Triangle)):
        gmap.points(i, mapped_pts)
        gmap.jacobians_dets_normals(i, jacobians, jdets, normals)
        gmap.jacobians(i, jacobians2)

        assert np.allclose(jacobians, jacobians2)

        for p, n in zip(mapped_pts, normals):
            assert p.dot(n) > 0.0

        for j, jdet, n in zip(jacobians, jdets, normals):
            assert np.allclose(n, np.cross(j[0], j[1]) / jdet)

        for jdet, j in zip(jdets, jacobians):
            assert np.isclose(jdet, np.linalg.norm(np.cross(j[0], j[1])))
