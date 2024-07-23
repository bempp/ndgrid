//! Entity geometry
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

/// The geometry of an entity
pub trait Geometry {
    /// Point Type
    type Point<'a>: Point
    where
        Self: 'a;

    /// Point iterator
    type PointIter<'a>: Iterator<Item = Self::Point<'a>>
    where
        Self: 'a;

    /// Points
    fn points(&self) -> Self::PointIter<'_>;

    /// Point count
    fn point_count(&self) -> usize;
}
