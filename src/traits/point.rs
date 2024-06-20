//! Point

use rlst::RlstScalar;

/// A point
pub trait Point {
    /// The floating point type used for coordinates
    type T: RlstScalar;

    /// Get the coordinates of the point
    fn coords(&self, data: &mut [Self::T]);
}
