//! Entity geometry

use super::Grid;

/// The geometry of an entity
pub trait Geometry {
    /// Grid type
    type Grid: Grid;

    /// Point iterator
    type PointIter<'a>: Iterator<Item = <Self::Grid as Grid>::Point<'a>>
    where
        Self: 'a;

    /// Points
    fn points(&self) -> Self::PointIter<'_>;

    /// Point count
    fn point_count(&self) -> usize;

    /// volume
    fn volume(&self) -> usize;
}
