//! A serial triangle grid

use rlst::prelude::*;

use crate::types::{IntegerArray2, RealScalar};

pub struct TriangleGrid<T: RealScalar> {
    points: DynamicArray<T, 2>,
    cells: IntegerArray2,
    topology: TriangleTopology,
}
