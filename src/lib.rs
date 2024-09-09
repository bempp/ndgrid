//! n-dimensional grid
#![cfg_attr(feature = "strict", deny(warnings), deny(unused_crate_dependencies))]
#![warn(missing_docs)]

pub mod bindings;
pub mod geometry;
pub mod grid;
mod io;
pub mod shapes;
pub mod topology;
pub mod traits;
pub mod types;

#[cfg(feature = "mpi")]
pub use grid::ParallelGrid;
pub use grid::{SingleElementGrid, SingleElementGridBuilder};
