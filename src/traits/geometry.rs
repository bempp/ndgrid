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
    /// Scalar type
    type T: RealScalar;

    /// Point Type
    type Point<'a>: Point<T = Self::T>
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

    /// Embedded superdegree
    fn degree(&self) -> usize;
}
