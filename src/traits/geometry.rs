//! Entity geometry
use crate::types::Scalar;

/// A point
pub trait Point {
    /// Scalar type
    type T: Scalar;

    /// The point's index
    fn index(&self) -> usize;

    /// Return the dimension of the point.
    fn dim(&self) -> usize;

    /// Get the coordinates of the point.
    fn coords(&self, data: &mut [Self::T]);
}

/// The geometry of an entity
///
/// The geometry contains information about all points that make up the entity.
pub trait Geometry {
    /// Scalar type
    type T: Scalar;

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
