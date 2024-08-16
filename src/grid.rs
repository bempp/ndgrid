//! Grids
#[cfg(feature = "mpi")]
pub mod parallel;
pub mod serial;
pub use serial::{SingleElementGrid, SingleElementGridBuilder};
