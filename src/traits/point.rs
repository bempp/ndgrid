//! Point

use super::Grid;

/// A point
pub trait Point {
    /// The floating point type used for coordinates
    type Grid: Grid;

    /// Return the dimension of the point.
    fn dim(&self) -> usize {
        <Self::Grid as Grid>::geometry_dim()
    }

    /// Get the coordinates of the point.
    fn coords(&self, data: &mut [<Self::Grid as Grid>::T]);
}
