//! Single element grid
use crate::{
    geometry::{GeometryMap, SingleElementEntityGeometryBorrowed, SingleElementGeometryBorrowed},
    topology::single_type::{SingleTypeEntityTopologyBorrowed, SingleTypeTopologyBorrowed},
    traits::{Entity, Grid},
    types::{Array2DBorrowed, Ownership, RealScalar},
};
use ndelement::{reference_cell, traits::FiniteElement, types::ReferenceCellType};
use rlst::{rlst_array_from_slice2, RawAccess};

/// Single element grid entity
#[derive(Debug)]
pub struct SingleElementGridEntityBorrowed<
    'a,
    T: RealScalar,
    E: FiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a SingleElementGridBorrowed<'a, T, E>,
    cell_index: usize,
    entity_dim: usize,
    entity_index: usize,
}

impl<'e, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>>
    SingleElementGridEntityBorrowed<'e, T, E>
{
    /// Create new
    pub fn new(
        grid: &'e SingleElementGridBorrowed<'e, T, E>,
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
impl<T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Entity
    for SingleElementGridEntityBorrowed<'_, T, E>
{
    type T = T;
    type EntityDescriptor = ReferenceCellType;
    type Topology<'a>
        = SingleTypeEntityTopologyBorrowed<'a>
    where
        Self: 'a;
    type Geometry<'a>
        = SingleElementEntityGeometryBorrowed<'a, T, E>
    where
        Self: 'a;
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
        SingleElementEntityGeometryBorrowed::new(
            &self.grid.geometry,
            self.cell_index,
            self.entity_dim,
            self.entity_index,
        )
    }
    fn topology(&self) -> Self::Topology<'_> {
        SingleTypeEntityTopologyBorrowed::new(
            &self.grid.topology,
            self.entity_type(),
            self.local_index(),
        )
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
pub struct SingleElementGridBorrowedEntityIter<
    'a,
    T: RealScalar,
    E: FiniteElement<CellType = ReferenceCellType, T = T>,
> {
    grid: &'a SingleElementGridBorrowed<'a, T, E>,
    dim: usize,
    index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>>
    SingleElementGridBorrowedEntityIter<'a, T, E>
{
    /// Create new
    pub fn new(grid: &'a SingleElementGridBorrowed<'a, T, E>, dim: usize) -> Self {
        Self {
            grid,
            dim,
            index: 0,
        }
    }
}
impl<'a, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Iterator
    for SingleElementGridBorrowedEntityIter<'a, T, E>
{
    type Item = SingleElementGridEntityBorrowed<'a, T, E>;

    fn next(&mut self) -> Option<SingleElementGridEntityBorrowed<'a, T, E>> {
        self.index += 1;
        self.grid.entity(self.dim, self.index - 1)
    }
}

/// Serial single element grid
#[derive(Debug)]
pub struct SingleElementGridBorrowed<
    'a,
    T: RealScalar,
    E: FiniteElement<CellType = ReferenceCellType, T = T>,
> {
    topology: SingleTypeTopologyBorrowed<'a>,
    geometry: SingleElementGeometryBorrowed<'a, T, E>,
}

impl<'a, T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>>
    SingleElementGridBorrowed<'a, T, E>
{
    #[allow(clippy::too_many_arguments)]
    /// Create new from raw data
    pub fn new(
        tdim: usize,
        ids: Vec<Option<&'a [usize]>>,
        entity_types: &'a [ReferenceCellType],
        entity_counts: &'a [usize],
        downward_connectivity: Vec<Vec<Array2DBorrowed<'a, usize>>>,
        upward_connectivity: Vec<Vec<Vec<&'a [usize]>>>,
        points: Array2DBorrowed<'a, T>,
        cells: Array2DBorrowed<'a, usize>,
        elements: Vec<E>,
    ) -> Self {
        let topology = SingleTypeTopologyBorrowed::new(
            tdim,
            ids,
            entity_types,
            entity_counts,
            downward_connectivity,
            upward_connectivity,
        );
        let geometry = SingleElementGeometryBorrowed::new(points, cells, elements);

        Self { topology, geometry }
    }
}

impl<T: RealScalar, E: FiniteElement<CellType = ReferenceCellType, T = T>> Grid
    for SingleElementGridBorrowed<'_, T, E>
{
    type T = T;
    type Entity<'a>
        = SingleElementGridEntityBorrowed<'a, T, E>
    where
        Self: 'a;
    type GeometryMap<'a>
        = GeometryMap<'a, T, Array2DBorrowed<'a, T>, Array2DBorrowed<'a, usize>>
    where
        Self: 'a;
    type EntityDescriptor = ReferenceCellType;
    type EntityIter<'a>
        = SingleElementGridBorrowedEntityIter<'a, T, E>
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
                Some(SingleElementGridEntityBorrowed::new(
                    self,
                    local_index,
                    dim,
                    0,
                ))
            } else {
                let cell = self.topology.upward_connectivity[dim][self.topology_dim() - dim - 1]
                    [local_index][0];
                let index = self.topology.downward_connectivity[self.topology_dim()][dim]
                    .r()
                    .slice(1, cell)
                    .data()
                    .iter()
                    .position(|&i| i == local_index)
                    .unwrap();
                Some(SingleElementGridEntityBorrowed::new(self, cell, dim, index))
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
        SingleElementGridBorrowedEntityIter::new(self, dim)
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        // TODO: store a HashMap rather than doing a find every time?
        self.topology.ids[dim].as_ref().map(|ids| {
            ids.iter()
                .find(|&i| *i == id)
                .map(|i| self.entity(dim, *i))?
        })?
    }

    fn geometry_map(&self, entity_type: ReferenceCellType, points: &[T]) -> Self::GeometryMap<'_> {
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
