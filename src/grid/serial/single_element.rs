//! A serial single element grid
//!
//! This grid uses a single element to represent the geometry of its cells

use rlst::prelude::*;
mod topology;

use crate::types::{Array2D, RealScalar};
pub use topology::{SingleElementCellTopology, SingleElementTopology};

/// Single element grid
pub struct SingleElementGrid<T: RealScalar> {
    _points: DynamicArray<T, 2>,
    _cells: Array2D<usize>,
    _topology: SingleElementTopology,
}
