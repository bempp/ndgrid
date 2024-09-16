//! Single element grid
#[cfg(feature = "serde")]
use crate::{
    geometry::single_element::SerializableGeometry,
    topology::serial::single_type::SerializableTopology, traits::ConvertToSerializable,
};
use crate::{
    geometry::{GeometryMap, SingleElementEntityGeometry, SingleElementGeometry},
    topology::serial::{SingleTypeEntityTopology, SingleTypeTopology},
    traits::{Entity, Grid},
    types::{Ownership, RealScalar},
};
use ndelement::{
    ciarlet::{CiarletElement, LagrangeElementFamily},
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::{Continuity, ReferenceCellType},
};
use rlst::{rlst_array_from_slice2, rlst_dynamic_array2, RawAccess, RawAccessMut};

/// Single element grid entity
#[derive(Debug)]
pub struct SingleElementGridEntity<
    'a,
    T: RealScalar,
    E: FiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a SingleElementGrid<T, E>,
    cell_index: usize,
    entity_dim: usize,
    entity_index: usize,
}

impl<'e, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>>
    SingleElementGridEntity<'e, T, E>
{
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
impl<'e, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Entity
    for SingleElementGridEntity<'e, T, E>
{
    type T = T;
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
    fn id(&self) -> Option<usize> {
        self.grid
            .topology
            .entity_id(self.entity_dim, self.local_index())
    }
}

/// Single element grid entity iterator
#[derive(Debug)]
pub struct SingleElementGridEntityIter<
    'a,
    T: RealScalar,
    E: FiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a SingleElementGrid<T, E>,
    dim: usize,
    index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>>
    SingleElementGridEntityIter<'a, T, E>
{
    /// Create new
    pub fn new(grid: &'a SingleElementGrid<T, E>, dim: usize) -> Self {
        Self {
            grid,
            dim,
            index: 0,
        }
    }
}
impl<'a, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Iterator
    for SingleElementGridEntityIter<'a, T, E>
{
    type Item = SingleElementGridEntity<'a, T, E>;

    fn next(&mut self) -> Option<SingleElementGridEntity<'a, T, E>> {
        self.index += 1;
        self.grid.entity(self.dim, self.index - 1)
    }
}

/// Serial single element grid
#[derive(Debug)]
pub struct SingleElementGrid<T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> {
    topology: SingleTypeTopology,
    geometry: SingleElementGeometry<T, E>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
#[serde(bound = "for<'de2> T: serde::Deserialize<'de2>")]
pub struct SerializableGrid<T: RealScalar + serde::Serialize>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    topology: SerializableTopology,
    geometry: SerializableGeometry<T>,
}

#[cfg(feature = "serde")]
impl<T: RealScalar + serde::Serialize> ConvertToSerializable
    for SingleElementGrid<T, CiarletElement<T>>
{
    type SerializableType = SerializableGrid<T>;
    fn to_serializable(&self) -> SerializableGrid<T> {
        SerializableGrid {
            topology: self.topology.to_serializable(),
            geometry: self.geometry.to_serializable(),
        }
    }
    fn from_serializable(s: SerializableGrid<T>) -> Self {
        Self {
            topology: SingleTypeTopology::from_serializable(s.topology),
            geometry: SingleElementGeometry::from_serializable(s.geometry),
        }
    }
}

impl<T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> SingleElementGrid<T, E> {
    /// Create new
    pub fn new(topology: SingleTypeTopology, geometry: SingleElementGeometry<T, E>) -> Self {
        Self { topology, geometry }
    }
}

impl<T: RealScalar> SingleElementGrid<T, CiarletElement<T>> {
    /// Create new from raw data
    pub fn new_from_raw_data(
        coordinates: &[T],
        gdim: usize,
        cells: &[usize],
        cell_type: ReferenceCellType,
        geometry_degree: usize,
    ) -> Self {
        let npts = coordinates.len() / gdim;
        let mut points = rlst_dynamic_array2!(T, [gdim, npts]);
        points.data_mut().copy_from_slice(coordinates);

        let family = LagrangeElementFamily::<T>::new(geometry_degree, Continuity::Standard);

        let geometry =
            SingleElementGeometry::<T, CiarletElement<T>>::new(cell_type, points, cells, &family);

        let points_per_cell = family.element(cell_type).dim();
        let tpoints_per_cell = reference_cell::entity_counts(cell_type)[0];
        let ncells = cells.len() / points_per_cell;

        let mut tcells = vec![0; tpoints_per_cell * ncells];
        for c in 0..ncells {
            tcells[c * tpoints_per_cell..(c + 1) * tpoints_per_cell].copy_from_slice(
                &cells[c * points_per_cell..c * points_per_cell + tpoints_per_cell],
            );
        }

        let topology = SingleTypeTopology::new(&tcells, cell_type, None, None);

        Self { topology, geometry }
    }
}

impl<T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Grid
    for SingleElementGrid<T, E>
{
    type T = T;
    type Entity<'a> = SingleElementGridEntity<'a, T, E> where Self: 'a;
    type GeometryMap<'a> = GeometryMap<'a, T> where Self: 'a;
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

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.topology.entity_types()[dim..dim + 1]
    }

    fn entity_count(&self, entity_type: ReferenceCellType) -> usize {
        self.topology.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        SingleElementGridEntityIter::new(self, dim)
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        // TODO: store a HashMap rather than doing a find every time?
        self.topology.ids[dim].as_ref().map(|ids| {
            ids.iter()
                .find(|&i| *i == id)
                .map(|i| self.entity(dim, *i))?
        })?
    }

    fn geometry_map(&self, entity_type: ReferenceCellType, points: &[T]) -> GeometryMap<'_, T> {
        let entity_dim = reference_cell::dim(entity_type);
        let npoints = points.len() / entity_dim;
        let rlst_points = rlst_array_from_slice2!(points, [entity_dim, npoints]);
        if entity_type == self.topology.entity_types()[self.topology_dim()] {
            GeometryMap::new(
                self.geometry.element(),
                &rlst_points,
                self.geometry.points(),
                self.geometry.cells(),
            )
        } else {
            unimplemented!();
        }
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
            SingleTypeTopology::new(&[0, 1, 2, 2, 1, 3], ReferenceCellType::Triangle, None, None),
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
