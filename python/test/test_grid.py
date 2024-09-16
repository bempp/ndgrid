import pytest
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


@pytest.mark.parametrize("level", range(4))
def test_entities(level):
    grid = regular_sphere(level)
    for i in range(grid.entity_count(ReferenceCellType.Triangle)):
        entity = grid.entity(2, i)
        assert entity.local_index == i
        if entity.id is not None:
            assert i == grid.entity_from_id(2, entity.id).local_index


@pytest.mark.parametrize("level", range(4))
def test_entity_types(level):
    etypes = [
        ReferenceCellType.Point,
        ReferenceCellType.Interval,
        ReferenceCellType.Triangle,
    ]
    grid = regular_sphere(level)
    for d, e in enumerate(etypes):
        assert grid.entity_types(d) == [e]

        for i in range(grid.entity_count(e)):
            entity = grid.entity(d, i)
            assert entity.entity_type == e


@pytest.mark.parametrize("level", range(4))
def test_ownership(level):
    etypes = [
        ReferenceCellType.Point,
        ReferenceCellType.Interval,
        ReferenceCellType.Triangle,
    ]
    grid = regular_sphere(level)
    for d, e in enumerate(etypes):
        for i in range(grid.entity_count(e)):
            entity = grid.entity(d, i)
            assert entity.ownership.is_owned


@pytest.mark.parametrize("level", range(4))
def test_subentities(level):
    etypes = [
        ReferenceCellType.Point,
        ReferenceCellType.Interval,
        ReferenceCellType.Triangle,
    ]
    grid = regular_sphere(level)
    connected_entities = [
        [[[] for _ in etypes] for i in range(grid.entity_count(e))] for e in etypes
    ]
    for d, etype in enumerate(
        [
            ReferenceCellType.Point,
            ReferenceCellType.Interval,
            ReferenceCellType.Triangle,
        ]
    ):
        for i in range(grid.entity_count(etype)):
            entity = grid.entity(d, i)
            for d2 in range(d + 1):
                for j, e in enumerate(entity.topology.sub_entities(d2)):
                    assert e == entity.topology.sub_entity(d2, j)
                    connected_entities[d2][e][d].append(i)

    for d, ce_d in enumerate(connected_entities):
        for i, ce_di in enumerate(ce_d):
            entity = grid.entity(d, i)
            for d2, ce_did2 in enumerate(ce_di):
                if d2 > d:
                    assert set(ce_did2) == set(entity.topology.connected_entities(d2))


@pytest.mark.parametrize("level", range(4))
def test_geometry(level):
    grid = regular_sphere(level)

    for cell in range(grid.entity_count(ReferenceCellType.Triangle)):
        entity = grid.entity(2, cell)
        assert entity.geometry.degree == 1
        points = entity.geometry.points()
        assert points.shape[0] == entity.geometry.point_count
        for pt in points:
            assert np.isclose(pt.dot(pt), 1.0)
