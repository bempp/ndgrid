//! Grids
mod builder;
pub mod local_grid;

pub use local_grid::{SingleElementGrid, SingleElementGridBorrowed, SingleElementGridBuilder};

use local_grid::LocalGrid;

#[cfg(feature = "serde")]
use crate::traits::{ConvertToSerializable, RONImportParallel};
use crate::{
    traits::{Grid, ParallelGrid as ParallelGridTrait},
    types::Ownership,
};
use mpi::traits::Communicator;

/// Parallel grid
#[derive(Debug)]
pub struct ParallelGrid<'a, C: Communicator, G: Grid + Sync> {
    comm: &'a C,
    local_grid: LocalGrid<G>,
}

unsafe impl<C: Communicator, GridImpl: Grid + Sync> Sync for ParallelGrid<'_, C, GridImpl> {}

impl<'a, C: Communicator, G: Grid + Sync> ParallelGrid<'a, C, G> {
    /// Create new
    pub fn new(
        comm: &'a C,
        serial_grid: G,
        ownership: Vec<Vec<Ownership>>,
        global_indices: Vec<Vec<usize>>,
    ) -> Self {
        Self {
            comm,
            local_grid: LocalGrid::new(serial_grid, ownership, global_indices),
        }
    }
}

#[cfg(feature = "serde")]
impl<'a, C: Communicator + 'a, G: Grid + Sync + ConvertToSerializable> RONImportParallel<'a, C>
    for ParallelGrid<'a, C, G>
where
    for<'de2> <G as ConvertToSerializable>::SerializableType: serde::Deserialize<'de2>,
    Self: 'a,
{
    fn create_from_ron_info(comm: &'a C, local_grid: LocalGrid<G>) -> Self {
        Self { comm, local_grid }
    }
}

impl<C: Communicator, G: Grid + Sync> ParallelGridTrait<C> for ParallelGrid<'_, C, G> {
    type LocalGrid<'a>
        = LocalGrid<G>
    where
        Self: 'a;
    fn comm(&self) -> &C {
        self.comm
    }
    fn local_grid(&self) -> &LocalGrid<G> {
        &self.local_grid
    }
}
impl<C: Communicator, G: Grid + Sync> Grid for ParallelGrid<'_, C, G> {
    type T = G::T;
    type Entity<'a>
        = <LocalGrid<G> as Grid>::Entity<'a>
    where
        Self: 'a;
    type GeometryMap<'a>
        = <LocalGrid<G> as Grid>::GeometryMap<'a>
    where
        Self: 'a;
    type EntityDescriptor = <LocalGrid<G> as Grid>::EntityDescriptor;
    type EntityIter<'a>
        = <LocalGrid<G> as Grid>::EntityIter<'a>
    where
        Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.local_grid.geometry_dim()
    }
    fn topology_dim(&self) -> usize {
        self.local_grid.topology_dim()
    }

    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>> {
        self.local_grid.entity(dim, local_index)
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.local_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.local_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        self.local_grid.entity_iter(dim)
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        self.local_grid.entity_from_id(dim, id)
    }

    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_> {
        self.local_grid.geometry_map(entity_type, points)
    }
}
