//! Grid geometry

use super::Grid;

/// The geometry of a grid
pub trait Geometry {
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
