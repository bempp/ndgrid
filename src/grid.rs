//! Grids
pub mod local_grid;
#[cfg(feature = "mpi")]
mod parallel_builder;
#[cfg(feature = "mpi")]
mod parallel_grid;

pub use local_grid::{MixedGrid, MixedGridBuilder, SingleElementGrid, SingleElementGridBuilder};
#[cfg(feature = "mpi")]
pub use parallel_grid::ParallelGridImpl;
