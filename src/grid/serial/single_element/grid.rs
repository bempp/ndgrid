//! Single element grid
use crate::{
    geometry::{Point, SingleElementEntityGeometry, SingleElementGeometry},
    topology::serial::{SingleTypeEntityTopology, SingleTypeTopology},
    traits::{Entity, Grid},
    types::{Ownership, RealScalar},
};
use ndelement::{traits::FiniteElement, types::ReferenceCellType};
use rlst::RawAccess;

/// Single element grid entity
pub struct SingleElementGridEntity<'a, T: RealScalar, E: FiniteElement> {
    grid: &'a SingleElementGrid<T, E>,
    cell_index: usize,
    entity_dim: usize,
    entity_index: usize,
}

impl<'e, T: RealScalar, E: FiniteElement> SingleElementGridEntity<'e, T, E> {
    /// Create new
    pub fn new(
        grid: &'e SingleElementGrid<T, E>,
        cell_index: usize,
        entity_dim: usize,
        entity_index: usize,
    ) -> Self {
        Self {
            grid,
            cell_index,
            entity_dim,
            entity_index,
        }
    }
}
impl<'e, T: RealScalar, E: FiniteElement> Entity for SingleElementGridEntity<'e, T, E> {
    type EntityDescriptor = ReferenceCellType;
    type Topology<'a> = SingleTypeEntityTopology<'a> where Self: 'a;
    type Geometry<'a> = SingleElementEntityGeometry<'a, T, E> where Self: 'a;
    fn entity_type(&self) -> ReferenceCellType {
        self.grid.topology.entity_types()[self.entity_dim]
    }
    fn local_index(&self) -> usize {
        self.grid
            .topology
            .cell_entity_index(self.cell_index, self.entity_dim, self.entity_index)
    }
    fn global_index(&self) -> usize {
        self.local_index()
    }
    fn geometry(&self) -> Self::Geometry<'_> {
        SingleElementEntityGeometry::new(
            &self.grid.geometry,
            self.cell_index,
            self.entity_dim,
            self.entity_index,
        )
    }
    fn topology(&self) -> Self::Topology<'_> {
        SingleTypeEntityTopology::new(&self.grid.topology, self.entity_type(), self.local_index())
    }
    fn ownership(&self) -> Ownership {
        Ownership::Owned
    }
}

/// Single element grid entity iterator
pub struct SingleElementGridEntityIter<'a, T: RealScalar, E: FiniteElement> {
    grid: &'a SingleElementGrid<T, E>,
    dim: usize,
    index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement> SingleElementGridEntityIter<'a, T, E> {
    /// Create new
    pub fn new(grid: &'a SingleElementGrid<T, E>, dim: usize) -> Self {
        Self {
            grid,
            dim,
            index: 0,
        }
    }
}
impl<'a, T: RealScalar, E: FiniteElement> Iterator for SingleElementGridEntityIter<'a, T, E> {
    type Item = SingleElementGridEntity<'a, T, E>;

    fn next(&mut self) -> Option<SingleElementGridEntity<'a, T, E>> {
        self.index += 1;
        self.grid.entity(self.dim, self.index - 1)
    }
}

/// Serial single element grid
pub struct SingleElementGrid<T: RealScalar, E: FiniteElement> {
    topology: SingleTypeTopology,
    geometry: SingleElementGeometry<T, E>,
}

impl<T: RealScalar, E: FiniteElement> SingleElementGrid<T, E> {
    /// Create new
    pub fn new(topology: SingleTypeTopology, geometry: SingleElementGeometry<T, E>) -> Self {
        Self { topology, geometry }
    }
}
impl<T: RealScalar, E: FiniteElement> Grid for SingleElementGrid<T, E> {
    type T = T;
    type Point<'a> = Point<'a, T> where Self: 'a;
    type Entity<'a> = SingleElementGridEntity<'a, T, E> where Self: 'a;
    type Topology<'a> = SingleTypeEntityTopology<'a> where Self: 'a;
    type Geometry<'a> = SingleElementEntityGeometry<'a, T, E> where Self: 'a;
    type EntityDescriptor = ReferenceCellType;
    type EntityIter<'a> = SingleElementGridEntityIter<'a, T, E>
    where
        Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.geometry.dim()
    }
    fn topology_dim(&self) -> usize {
        self.topology.dim()
    }

    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>> {
        if local_index
            < self
                .topology
                .entity_count(self.topology.entity_types()[dim])
        {
            if dim == self.topology_dim() {
                Some(SingleElementGridEntity::new(self, local_index, dim, 0))
            } else {
                let cell = self.topology.upward_connectivity[dim][self.topology_dim() - dim - 1]
                    [local_index][0];
                let index = self.topology.downward_connectivity[self.topology_dim()][dim]
                    .view()
                    .slice(1, cell)
                    .data()
                    .iter()
                    .position(|&i| i == local_index)
                    .unwrap();
                Some(SingleElementGridEntity::new(self, cell, dim, index))
            }
        } else {
            None
        }
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        SingleElementGridEntityIter::new(self, dim)
    }

    fn entity_from_id(&self, _dim: usize, _id: usize) -> Option<Self::Entity<'_>> {
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::traits::Topology;
    use itertools::izip;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        reference_cell,
        types::Continuity,
    };
    use rlst::{rlst_dynamic_array2, RandomAccessMut};

    fn example_grid_triangle() -> SingleElementGrid<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [3, 4]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        *points.get_mut([0, 3]).unwrap() = 2.0;
        *points.get_mut([1, 3]).unwrap() = 1.0;
        *points.get_mut([2, 3]).unwrap() = 0.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGrid::new(
            SingleTypeTopology::new(&[0, 1, 2, 2, 1, 3], ReferenceCellType::Triangle),
            SingleElementGeometry::<f64, CiarletElement<f64>>::new(
                ReferenceCellType::Triangle,
                points,
                &[0, 1, 2, 2, 1, 3],
                &family,
            ),
        )
    }

    #[test]
    fn test_edges_triangle() {
        let grid = example_grid_triangle();
        let conn = reference_cell::connectivity(
            grid.entity(grid.topology_dim(), 0).unwrap().entity_type(),
        );
        for edge in grid.entity_iter(1) {
            let cell = grid.entity(grid.topology_dim(), edge.cell_index).unwrap();
            for (i, v) in izip!(
                &conn[1][edge.entity_index][0],
                edge.topology().sub_entity_iter(0)
            ) {
                assert_eq!(v, cell.topology().sub_entity(0, *i));
            }
        }
    }
}
