//! A serial single element grid
//!
//! This grid uses a single element to represent the geometry of its cells

use rlst::prelude::*;
mod topology;

use crate::types::{IntegerArray2, RealScalar};
use topology::SingleElementTopology;

pub struct SingleElementGrid<T: RealScalar> {
    points: DynamicArray<T, 2>,
    cells: IntegerArray2,
    topology: SingleElementTopology,
}
