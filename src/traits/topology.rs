use super::Grid;

/// The topology of an entity

pub trait Topology {
    /// Grid type
    type Grid: Grid;

    /// Entity iterator
    type EntityIter<'a>: Iterator<Item = <Self::Grid as Grid>::Entity<'a>>
    where
        Self: 'a;

    /// Iterator over indices of connected entities
    fn connected_entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// Iterator over sub-entities
    fn sub_entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// A sub-entity of this entity
    fn sub_entity(&self, dim: usize, index: usize) -> <Self::Grid as Grid>::Entity<'_>;

}
