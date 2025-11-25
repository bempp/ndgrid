//! Parallel grid
#[cfg(feature = "serde")]
use crate::traits::{ConvertToSerializable, RONImportParallel};
use crate::{
    grid::local_grid::LocalGrid,
    traits::{Grid, ParallelGrid},
    types::{Ownership, Scalar},
};
use mpi::traits::Communicator;
use rlst::distributed_tools::IndexLayout;
#[cfg(feature = "serde")]
use std::hash::Hash;
use std::{collections::HashMap, fmt::Debug};

/// Parallel grid
#[derive(Debug)]
pub struct ParallelGridImpl<'a, C: Communicator, G: Grid + Sync> {
    comm: &'a C,
    local_grid: LocalGrid<G>,
    cell_layout: std::rc::Rc<IndexLayout<'a, C>>,
}

impl<'a, C: Communicator, G: Grid + Sync> ParallelGridImpl<'a, C, G> {
    /// Create new
    pub fn new(
        comm: &'a C,
        serial_grid: G,
        ownership: HashMap<G::EntityDescriptor, Vec<Ownership>>,
        global_indices: HashMap<G::EntityDescriptor, Vec<usize>>,
    ) -> Self {
        let local_grid = LocalGrid::new(serial_grid, ownership, global_indices);
        let owned_cell_count = local_grid.owned_cell_count();
        Self {
            comm,
            local_grid,
            cell_layout: std::rc::Rc::new(IndexLayout::from_local_counts(owned_cell_count, comm)),
        }
    }
}

#[cfg(feature = "serde")]
impl<
    'a,
    EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash + serde::Serialize,
    C: Communicator + 'a,
    G: Grid<EntityDescriptor = EntityDescriptor> + Sync + ConvertToSerializable,
> RONImportParallel<'a, C> for ParallelGridImpl<'a, C, G>
where
    for<'de2> <G as ConvertToSerializable>::SerializableType: serde::Deserialize<'de2>,
    for<'de2> EntityDescriptor: serde::Deserialize<'de2>,
    Self: 'a,
{
    fn create_from_ron_info(comm: &'a C, local_grid: LocalGrid<G>) -> Self {
        let owned_cell_count = local_grid.owned_cell_count();
        Self {
            comm,
            local_grid,
            cell_layout: std::rc::Rc::new(IndexLayout::from_local_counts(owned_cell_count, comm)),
        }
    }
}

impl<T: Scalar, C: Communicator, G: Grid<T = T> + Sync> ParallelGrid
    for ParallelGridImpl<'_, C, G>
{
    type LocalGrid = LocalGrid<G>;

    type C = C;

    type T = T;
    fn comm(&self) -> &C {
        self.comm
    }
    fn local_grid(&self) -> &Self::LocalGrid {
        &self.local_grid
    }

    fn cell_layout(&self) -> std::rc::Rc<IndexLayout<'_, C>> {
        self.cell_layout.clone()
    }
}
