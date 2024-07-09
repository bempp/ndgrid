//! The topology of an entity

/// The topology of an entity
pub trait Topology {
    /// Entity iterator
    type EntityIndexIter<'a>: Iterator<Item = usize>
    where
        Self: 'a;

    /// Entity iterator
    type ConnectedEntityIndexIter<'a>: Iterator<Item = usize>
    where
        Self: 'a;

    /// Iterator over indices of connected entities
    fn connected_entity_iter(&self, dim: usize) -> Self::ConnectedEntityIndexIter<'_>;

    /// Iterator over sub-entities
    fn sub_entity_iter(&self, dim: usize) -> Self::EntityIndexIter<'_>;

    /// An index of a sub-entity of this entity
    fn sub_entity(&self, dim: usize, index: usize) -> usize;
}
