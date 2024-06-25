//! The topology of an entity
use crate::traits::Entity;

/// The topology of an entity
pub trait Topology {
    /// Type used as identifier of different entity types
    type Entity<'a>: Entity<Topology<'a> = Self>
    where
        Self: 'a;

    /// Entity iterator
    type EntityIter<'a>: Iterator<Item = Self::Entity<'a>>
    where
        Self: 'a;

    /// Entity iterator
    type ConnectedEntityIter<'a>: Iterator<Item = Self::Entity<'a>>
    where
        Self: 'a;

    /// Iterator over indices of connected entities
    fn connected_entity_iter(&self, dim: usize) -> Self::ConnectedEntityIter<'_>;

    /// Iterator over sub-entities
    fn sub_entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// A sub-entity of this entity
    fn sub_entity(&self, dim: usize, index: usize) -> Self::Entity<'_>;
}
