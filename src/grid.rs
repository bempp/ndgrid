//! Grids
#[cfg(feature = "mpi")]
pub mod parallel;
pub mod serial;

#[cfg(feature = "mpi")]
pub use parallel::ParallelGrid;
pub use serial::{SingleElementGrid, SingleElementGridBuilder};
