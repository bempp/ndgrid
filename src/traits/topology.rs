use super::Grid;

/// The topology of a grid

pub trait Topology {
    /// Grid type
    type Grid: Grid;

    /// Entity iterator
    type EntityIter<'a>: Iterator<Item = <Self::Grid as Grid>::Entity<'a>>
    where
        Self: 'a;

    /// Iterator over indices of connected entities
    fn connected_entities(&self, dim: usize) -> Self::EntityIter<'_>;
}
