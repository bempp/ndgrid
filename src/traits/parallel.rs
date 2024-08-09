//! MPI parallelised grid

use crate::traits::{Builder, Grid};
use mpi::traits::Communicator;
use std::collections::HashMap;

pub trait ParallelBuilder: Builder {
    //! Trait to build a MPI parellelised mesh from a grid builder.

    /// The type of the parallel grid that the builder creates
    type ParallelGrid<'a, C: Communicator + 'a>: ParallelGrid;

    /// Create the grid on the main process
    fn create_parallel_grid<'a, C: Communicator>(
        self,
        comm: &'a C,
        cell_owners: &HashMap<usize, usize>,
    ) -> Self::ParallelGrid<'a, C>;

    /// Create the grid on a subprocess
    fn receive_parallel_grid<C: Communicator>(
        self,
        comm: &C,
        root_process: usize,
    ) -> Self::ParallelGrid<'_, C>;
}

pub trait ParallelGrid {
    //! An MPI parallelised grid
    /// The MPI communicator type
    type Comm: Communicator;
    /// The type of the subgrid on each process
    type LocalGrid: Grid;

    /// The MPI communicator
    fn comm(&self) -> &Self::Comm;

    /// The subgrid on the process
    fn local_grid(&self) -> &Self::LocalGrid;
}
