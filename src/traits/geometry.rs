//! Entity geometry
use crate::traits::Point;

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

    /// volume
    fn volume(&self) -> usize;
}
