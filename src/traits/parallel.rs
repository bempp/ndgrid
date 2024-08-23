//! Traits for MPI parallel grids
use super::{Builder, Grid};
use mpi::traits::Communicator;

pub trait ParallelBuilder: Builder {
    //! MPI parallel grid builder

    /// Parallel grid type
    type ParallelGrid<'a, C: Communicator + 'a>: Grid + ParallelGrid<C>
    where
        Self: 'a;

    /// Create a parallel grid
    fn create_parallel_grid<'a, C: Communicator>(&self, comm: &'a C) -> Self::ParallelGrid<'a, C>;

    /// Receive a parallel grid
    fn receive_parallel_grid<'a, C: Communicator>(
        &self,
        comm: &'a C,
        root_rank: i32,
    ) -> Self::ParallelGrid<'a, C>;
}

pub trait ParallelGrid<C: Communicator>: Grid {
    //! MPI parallel grid

    /// Local grid type
    type LocalGrid<'a>: Sync
        + Grid<T = <Self as Grid>::T, EntityDescriptor = <Self as Grid>::EntityDescriptor>
    where
        Self: 'a;

    /// MPI communicator
    fn comm(&self) -> &C;
    /// Local grid on the current process
    fn local_grid(&self) -> &Self::LocalGrid<'_>;
}
