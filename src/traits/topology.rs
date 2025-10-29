//! The topology of an entity
use std::{fmt::Debug, hash::Hash};

/// The topology of an entity
pub trait Topology {
    /// Entity type enum
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Entity iterator
    type EntityIndexIter<'a>: Iterator<Item = usize>
    where
        Self: 'a;

    /// Entity iterator
    type ConnectedEntityIndexIter<'a>: Iterator<Item = usize>
    where
        Self: 'a;

    /// Iterator over indices of connected entities
    fn connected_entity_iter(
        &self,
        entity_type: Self::EntityDescriptor,
    ) -> Self::ConnectedEntityIndexIter<'_>;

    /// Iterator over sub-entities
    fn sub_entity_iter(&self, entity_type: Self::EntityDescriptor) -> Self::EntityIndexIter<'_>;

    /// An index of a sub-entity of this entity
    fn sub_entity(&self, entity_type: Self::EntityDescriptor, index: usize) -> usize;

    /// A 32-bit integer that encodes the orientation differences between this entity and the corresponding reference entity
    fn orientation(&self) -> i32;
}
