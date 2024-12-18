//! Grids
pub mod parallel;
pub mod serial;

pub use parallel::ParallelGrid;
pub use serial::{SingleElementGrid, SingleElementGridBorrowed, SingleElementGridBuilder};
