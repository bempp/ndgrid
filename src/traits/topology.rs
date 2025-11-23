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

    /// Iterator over indices of connected entities.
    ///
    /// This iterator returns all entities of a given type connected to this entity.
    fn connected_entity_iter(
        &self,
        entity_type: Self::EntityDescriptor,
    ) -> Self::ConnectedEntityIndexIter<'_>;

    /// Iterator over sub-entities
    ///
    /// This iterator iterates over all subentities of a given type. It returns the index of the corresponding entity.
    /// This index can then be used with the [Grid::entity](crate::traits::Grid::entity) to return the corresponding actual entity.
    fn sub_entity_iter(&self, entity_type: Self::EntityDescriptor) -> Self::EntityIndexIter<'_>;

    /// An index of a sub-entity of this entity
    ///
    /// Assume that `e` is a triangle. Then calling `e.sub_entity(ReferenceCellType::Point, 0)` returns the entity index
    /// of the first vertex of that entity.
    fn sub_entity(&self, entity_type: Self::EntityDescriptor, index: usize) -> usize;

    /// A 32-bit integer that encodes the orientation differences between this entity and the corresponding reference entity
    fn orientation(&self) -> i32;
}
