//! Grids
mod builder;
pub mod local_grid;

pub use local_grid::{SingleElementGrid, SingleElementGridBorrowed, SingleElementGridBuilder};

use local_grid::LocalGrid;

#[cfg(feature = "serde")]
use crate::traits::{ConvertToSerializable, RONImportParallel};
use crate::{
    traits::{Grid, ParallelGrid},
    types::Ownership,
};
use mpi::traits::Communicator;

/// Parallel grid
#[derive(Debug)]
pub struct ParallelGridImpl<'a, C: Communicator, G: Grid + Sync> {
    comm: &'a C,
    local_grid: LocalGrid<G>,
    cell_index_layout: std::rc::Rc<bempp_distributed_tools::IndexLayout<'a, C>>,
}

impl<'a, C: Communicator, G: Grid + Sync> ParallelGridImpl<'a, C, G> {
    /// Create new
    pub fn new(
        comm: &'a C,
        serial_grid: G,
        ownership: Vec<Vec<Ownership>>,
        global_indices: Vec<Vec<usize>>,
    ) -> Self {
        let local_grid = LocalGrid::new(serial_grid, ownership, global_indices);
        let owned_cell_count = local_grid.owned_cell_count();
        Self {
            comm,
            local_grid,
            cell_index_layout: std::rc::Rc::new(
                bempp_distributed_tools::IndexLayout::from_local_counts(owned_cell_count, comm),
            ),
        }
    }
}

#[cfg(feature = "serde")]
impl<'a, C: Communicator + 'a, G: Grid + Sync + ConvertToSerializable> RONImportParallel<'a, C>
    for ParallelGridImpl<'a, C, G>
where
    for<'de2> <G as ConvertToSerializable>::SerializableType: serde::Deserialize<'de2>,
    Self: 'a,
{
    fn create_from_ron_info(comm: &'a C, local_grid: LocalGrid<G>) -> Self {
        let owned_cell_count = local_grid.owned_cell_count();
        Self {
            comm,
            local_grid,
            cell_index_layout: std::rc::Rc::new(
                bempp_distributed_tools::IndexLayout::from_local_counts(owned_cell_count, comm),
            ),
        }
    }
}

impl<C: Communicator, G: Grid + Sync> ParallelGrid for ParallelGridImpl<'_, C, G> {
    type LocalGrid = LocalGrid<G>;
    type C = C;

    fn comm(&self) -> &C {
        self.comm
    }
    fn local_grid(&self) -> &LocalGrid<G> {
        &self.local_grid
    }

    fn cell_index_layout(&self) -> std::rc::Rc<bempp_distributed_tools::IndexLayout<'_, C>> {
        self.cell_index_layout.clone()
    }
}
