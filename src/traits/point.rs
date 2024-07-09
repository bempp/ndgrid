//! Point
use crate::types::RealScalar;

/// A point
pub trait Point {
    /// Scalar type
    type T: RealScalar;

    /// Return the dimension of the point.
    fn dim(&self) -> usize;

    /// Get the coordinates of the point.
    fn coords(&self, data: &mut [Self::T]);
}
