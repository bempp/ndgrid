//! A serial single element grid
//!
//! This grid uses a single element to represent the geometry of its cells

use rlst::prelude::*;
mod topology;

use crate::types::{Array2D, RealScalar};
pub use topology::{SingleElementTopology, SingleElementCellTopology};

/// Single element grid
pub struct SingleElementGrid<T: RealScalar> {
    points: DynamicArray<T, 2>,
    cells: Array2D<usize>,
    topology: SingleElementTopology,
}
