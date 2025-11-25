//! Mixed grid
#[cfg(feature = "mpi")]
use crate::ParallelGridImpl;
#[cfg(feature = "mpi")]
use crate::{
    MixedGridBuilder,
    traits::{Builder, DistributableGrid, ParallelBuilder},
    types::GraphPartitioner,
};
#[cfg(feature = "serde")]
use crate::{
    geometry::mixed::SerializableGeometry, topology::mixed::SerializableTopology,
    traits::ConvertToSerializable,
};
use crate::{
    geometry::{GeometryMap, MixedEntityGeometry, MixedGeometry},
    topology::mixed::{MixedEntityTopology, MixedTopology},
    traits::{Entity, Grid},
    types::{Ownership, Scalar},
};
use itertools::izip;
#[cfg(feature = "mpi")]
use mpi::traits::{Communicator, Equivalence};
use ndelement::{
    ciarlet::{CiarletElement, LagrangeElementFamily},
    map::IdentityMap,
    reference_cell,
    traits::{ElementFamily, FiniteElement, MappedFiniteElement},
    types::{Continuity, ReferenceCellType},
};
use rlst::dense::{base_array::BaseArray, data_container::VectorContainer};
use rlst::{SliceArray, rlst_dynamic_array};
use std::collections::HashMap;

/// Mixed grid entity
#[derive(Debug)]
pub struct MixedGridEntity<
    'a,
    T: Scalar,
    E: MappedFiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a MixedGrid<T, E>,
    cell_type: ReferenceCellType,
    cell_index: usize,
    entity_type: ReferenceCellType,
    entity_index: usize,
    geometry_element_index: usize,
    geometry_cell_index: usize,
}

impl<'e, T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>>
    MixedGridEntity<'e, T, E>
{
    /// Create new
    pub fn new(
        grid: &'e MixedGrid<T, E>,
        cell_type: ReferenceCellType,
        cell_index: usize,
        entity_type: ReferenceCellType,
        entity_index: usize,
    ) -> Self {
        Self {
            grid,
            cell_type,
            cell_index,
            entity_type,
            entity_index,
            geometry_element_index: grid.geometry.insertion_indices_to_element_indices
                [grid.topology.insertion_indices[&cell_type][cell_index]],
            geometry_cell_index: grid.geometry.insertion_indices_to_cell_indices
                [grid.topology.insertion_indices[&cell_type][cell_index]],
        }
    }
}
impl<T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>> Entity
    for MixedGridEntity<'_, T, E>
{
    type T = T;
    type EntityDescriptor = ReferenceCellType;
    type Topology<'a>
        = MixedEntityTopology<'a>
    where
        Self: 'a;
    type Geometry<'a>
        = MixedEntityGeometry<'a, T, E>
    where
        Self: 'a;
    fn entity_type(&self) -> ReferenceCellType {
        self.entity_type
    }
    fn local_index(&self) -> usize {
        self.grid.topology.cell_entity_index(
            self.cell_type,
            self.cell_index,
            self.entity_type,
            self.entity_index,
        )
    }
    fn global_index(&self) -> usize {
        self.local_index()
    }
    fn geometry(&self) -> Self::Geometry<'_> {
        MixedEntityGeometry::new(
            &self.grid.geometry,
            self.geometry_element_index,
            self.geometry_cell_index,
            reference_cell::dim(self.entity_type),
            self.entity_index,
        )
    }
    fn topology(&self) -> Self::Topology<'_> {
        MixedEntityTopology::new(&self.grid.topology, self.entity_type(), self.local_index())
    }
    fn ownership(&self) -> Ownership {
        Ownership::Owned
    }
    fn id(&self) -> Option<usize> {
        self.grid
            .topology
            .entity_id(self.entity_type, self.local_index())
    }
}

/// Mixed grid entity iterator
#[derive(Debug)]
pub struct MixedGridEntityIter<
    'a,
    T: Scalar,
    E: MappedFiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a MixedGrid<T, E>,
    entity_type: ReferenceCellType,
    index: usize,
}

impl<'a, T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>>
    MixedGridEntityIter<'a, T, E>
{
    /// Create new
    pub fn new(grid: &'a MixedGrid<T, E>, entity_type: ReferenceCellType) -> Self {
        Self {
            grid,
            entity_type,
            index: 0,
        }
    }
}
impl<'a, T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>> Iterator
    for MixedGridEntityIter<'a, T, E>
{
    type Item = MixedGridEntity<'a, T, E>;

    fn next(&mut self) -> Option<MixedGridEntity<'a, T, E>> {
        self.index += 1;
        self.grid.entity(self.entity_type, self.index - 1)
    }
}

/// Serial mixed element grid
#[derive(Debug)]
pub struct MixedGrid<T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>> {
    topology: MixedTopology,
    geometry: MixedGeometry<T, E>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
#[serde(bound = "for<'de2> T: serde::Deserialize<'de2>")]
pub struct SerializableGrid<T: Scalar + serde::Serialize>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    topology: SerializableTopology,
    geometry: SerializableGeometry<T>,
}

#[cfg(feature = "serde")]
impl<T: Scalar + serde::Serialize> ConvertToSerializable
    for MixedGrid<T, CiarletElement<T, IdentityMap, T>>
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
            topology: MixedTopology::from_serializable(s.topology),
            geometry: MixedGeometry::from_serializable(s.geometry),
        }
    }
}

impl<T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>> MixedGrid<T, E> {
    /// Create new
    pub fn new(topology: MixedTopology, geometry: MixedGeometry<T, E>) -> Self {
        Self { topology, geometry }
    }
}

impl<T: Scalar> MixedGrid<T, CiarletElement<T, IdentityMap, T>> {
    /// Create new from raw data
    pub fn new_from_raw_data(
        coordinates: &[T],
        gdim: usize,
        cells: &[usize],
        cell_types: &[ReferenceCellType],
        cell_degrees: &[usize],
    ) -> Self {
        let npts = coordinates.len() / gdim;
        let mut points = rlst_dynamic_array!(T, [gdim, npts]);
        points.data_mut().unwrap().copy_from_slice(coordinates);

        let mut element_families = vec![];
        let mut element_family_indices = HashMap::new();

        let cell_families = cell_degrees
            .iter()
            .map(|d| {
                *element_family_indices.entry(*d).or_insert_with(|| {
                    let index = element_families.len();
                    element_families
                        .push(LagrangeElementFamily::<T, T>::new(*d, Continuity::Standard));
                    index
                })
            })
            .collect::<Vec<_>>();

        let geometry = MixedGeometry::<T, CiarletElement<T, IdentityMap, T>>::new(
            cell_types,
            points,
            cells,
            &element_families,
            &cell_families,
        );

        let mut start = 0;
        let mut tcells = vec![];
        for (t, f) in izip!(cell_types, &cell_families) {
            let tpoints_per_cell = reference_cell::entity_counts(*t)[0];
            for i in 0..tpoints_per_cell {
                tcells.push(cells[start + i]);
            }
            start += element_families[*f].element(*t).dim();
        }

        let topology = MixedTopology::new(&tcells, cell_types, None, None);

        Self { topology, geometry }
    }
}

impl<T: Scalar, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>> Grid
    for MixedGrid<T, E>
{
    type T = T;
    type Entity<'a>
        = MixedGridEntity<'a, T, E>
    where
        Self: 'a;
    type GeometryMap<'a>
        = GeometryMap<'a, T, BaseArray<VectorContainer<T>, 2>, BaseArray<VectorContainer<usize>, 2>>
    where
        Self: 'a;
    type EntityDescriptor = ReferenceCellType;
    type EntityIter<'a>
        = MixedGridEntityIter<'a, T, E>
    where
        Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.geometry.dim()
    }
    fn topology_dim(&self) -> usize {
        self.topology.dim()
    }

    fn entity(
        &self,
        entity_type: ReferenceCellType,
        local_index: usize,
    ) -> Option<Self::Entity<'_>> {
        let dim = reference_cell::dim(entity_type);
        if local_index < self.topology.entity_count(entity_type) {
            if dim == self.topology_dim() {
                Some(MixedGridEntity::new(
                    self,
                    entity_type,
                    local_index,
                    entity_type,
                    0,
                ))
            } else {
                for t in &self.topology.entity_types()[self.topology_dim()] {
                    if let Some(cell) =
                        self.topology.upward_connectivity[&entity_type][t][local_index].first()
                    {
                        let dc = &self.topology.downward_connectivity[t][&entity_type];
                        if let Some(index) =
                            (0..dc.shape()[0]).position(|i| dc[[i, *cell]] == local_index)
                        {
                            return Some(MixedGridEntity::new(self, *t, *cell, entity_type, index));
                        }
                    }
                }
                None
            }
        } else {
            None
        }
    }

    fn entity_types(&self, dim: usize) -> &[ReferenceCellType] {
        &self.topology.entity_types()[dim]
    }

    fn entity_count(&self, entity_type: ReferenceCellType) -> usize {
        self.topology.entity_count(entity_type)
    }

    fn entity_iter(&self, entity_type: ReferenceCellType) -> Self::EntityIter<'_> {
        MixedGridEntityIter::new(self, entity_type)
    }

    fn entity_from_id(
        &self,
        entity_type: ReferenceCellType,
        id: usize,
    ) -> Option<Self::Entity<'_>> {
        self.topology.ids_to_indices[&entity_type]
            .get(&id)
            .map(|i| self.entity(entity_type, *i))?
    }

    fn geometry_map(
        &self,
        entity_type: ReferenceCellType,
        geometry_degree: usize,
        points: &[T],
    ) -> GeometryMap<'_, T, BaseArray<VectorContainer<T>, 2>, BaseArray<VectorContainer<usize>, 2>>
    {
        let entity_dim = reference_cell::dim(entity_type);
        let npoints = points.len() / entity_dim;
        let rlst_points = SliceArray::<T, 2>::from_shape(points, [entity_dim, npoints]);

        for i in 0..self.geometry.element_count() {
            let e = self.geometry.element(i);
            if e.cell_type() == entity_type && e.lagrange_superdegree() == geometry_degree {
                return GeometryMap::new(
                    e,
                    &rlst_points,
                    self.geometry.points(),
                    self.geometry.cells(i),
                );
            }
        }
        unimplemented!();
    }
}

#[cfg(feature = "mpi")]
impl<T: Scalar + Equivalence, E: MappedFiniteElement<CellType = ReferenceCellType, T = T>>
    DistributableGrid for MixedGrid<T, E>
{
    type ParallelGrid<'a, C: Communicator + 'a> =
        ParallelGridImpl<'a, C, MixedGrid<T, CiarletElement<T, IdentityMap, T>>>;

    fn distribute<'a, C: Communicator>(
        &self,
        comm: &'a C,
        partitioner: GraphPartitioner,
    ) -> Self::ParallelGrid<'a, C> {
        let mut b = MixedGridBuilder::<T>::new(self.geometry.dim());
        let pts = self.geometry.points();
        for p in 0..pts.shape()[1] {
            b.add_point(
                p,
                &pts.data().unwrap()[p * pts.shape()[0]..(p + 1) * pts.shape()[0]],
            );
        }

        for (c, (element_i, cell_i)) in izip!(
            &self.geometry.insertion_indices_to_element_indices,
            &self.geometry.insertion_indices_to_cell_indices
        )
        .enumerate()
        {
            let e = self.geometry.element(*element_i);
            let cells = self.geometry.cells(*element_i);
            b.add_cell(
                c,
                (
                    e.cell_type(),
                    e.lagrange_superdegree(),
                    &cells.data().unwrap()
                        [cell_i * cells.shape()[0]..(cell_i + 1) * cells.shape()[0]],
                ),
            );
        }
        b.create_parallel_grid_root(comm, partitioner)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::traits::{GeometryMap, Topology};
    use approx::*;
    use itertools::izip;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        reference_cell,
        types::Continuity,
    };
    use rlst::rlst_dynamic_array;

    fn example_grid_triangle() -> MixedGrid<f64, CiarletElement<f64, IdentityMap>> {
        let mut points = rlst_dynamic_array!(f64, [3, 4]);
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
        MixedGrid::new(
            MixedTopology::new(
                &[0, 1, 2, 2, 1, 3],
                &[ReferenceCellType::Triangle, ReferenceCellType::Triangle],
                None,
                None,
            ),
            MixedGeometry::<f64, CiarletElement<f64, IdentityMap>>::new(
                &[ReferenceCellType::Triangle, ReferenceCellType::Triangle],
                points,
                &[0, 1, 2, 2, 1, 3],
                &[family],
                &[0, 0],
            ),
        )
    }

    fn example_grid_mixed() -> MixedGrid<f64, CiarletElement<f64, IdentityMap>> {
        let mut points = rlst_dynamic_array!(f64, [2, 8]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 2.0;
        *points.get_mut([1, 2]).unwrap() = 0.0;
        *points.get_mut([0, 3]).unwrap() = 4.0;
        *points.get_mut([1, 3]).unwrap() = 0.0;
        *points.get_mut([0, 4]).unwrap() = 0.0;
        *points.get_mut([1, 4]).unwrap() = 1.0;
        *points.get_mut([0, 5]).unwrap() = 1.0;
        *points.get_mut([1, 5]).unwrap() = 1.0;
        *points.get_mut([0, 6]).unwrap() = 2.0;
        *points.get_mut([1, 6]).unwrap() = 1.0;
        *points.get_mut([0, 7]).unwrap() = 4.0;
        *points.get_mut([1, 7]).unwrap() = 1.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        MixedGrid::new(
            MixedTopology::new(
                &[0, 5, 4, 0, 1, 5, 1, 2, 5, 6, 2, 3, 6, 7],
                &[
                    ReferenceCellType::Triangle,
                    ReferenceCellType::Triangle,
                    ReferenceCellType::Quadrilateral,
                    ReferenceCellType::Quadrilateral,
                ],
                None,
                None,
            ),
            MixedGeometry::<f64, CiarletElement<f64, IdentityMap>>::new(
                &[
                    ReferenceCellType::Triangle,
                    ReferenceCellType::Triangle,
                    ReferenceCellType::Quadrilateral,
                    ReferenceCellType::Quadrilateral,
                ],
                points,
                &[0, 5, 4, 0, 1, 5, 1, 2, 5, 6, 2, 3, 6, 7],
                &[family],
                &[0, 0, 0, 0],
            ),
        )
    }

    #[test]
    fn test_edges_triangle() {
        let grid = example_grid_triangle();
        let conn = reference_cell::connectivity(ReferenceCellType::Triangle);
        for edge in grid.entity_iter(ReferenceCellType::Interval) {
            let cell = grid
                .entity(ReferenceCellType::Triangle, edge.cell_index)
                .unwrap();
            for (i, v) in izip!(
                &conn[1][edge.entity_index][0],
                edge.topology().sub_entity_iter(ReferenceCellType::Point)
            ) {
                assert_eq!(v, cell.topology().sub_entity(ReferenceCellType::Point, *i));
            }
        }
    }

    #[test]
    fn test_geometry_map() {
        let grid = example_grid_mixed();
        let mut mapped_pts = vec![0.0; 4];

        let pts = vec![0.0, 0.0, 0.5, 0.5];
        let map = grid.geometry_map(ReferenceCellType::Triangle, 1, &pts);

        map.physical_points(0, &mut mapped_pts);
        assert_relative_eq!(mapped_pts[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[2], 0.5, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[3], 1.0, epsilon = 1e-10);

        map.physical_points(1, &mut mapped_pts);
        assert_relative_eq!(mapped_pts[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[2], 1.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[3], 0.5, epsilon = 1e-10);

        let pts = vec![0.5, 0.0, 1.0, 1.0];
        let map = grid.geometry_map(ReferenceCellType::Quadrilateral, 1, &pts);

        map.physical_points(0, &mut mapped_pts);
        assert_relative_eq!(mapped_pts[0], 1.5, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[2], 2.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[3], 1.0, epsilon = 1e-10);

        map.physical_points(1, &mut mapped_pts);
        assert_relative_eq!(mapped_pts[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[2], 4.0, epsilon = 1e-10);
        assert_relative_eq!(mapped_pts[3], 1.0, epsilon = 1e-10);
    }
}
