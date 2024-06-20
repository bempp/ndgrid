//! Traits for a mesh entity
use super::Grid;
use crate::types::Ownership;
use std::iter::Iterator;

/// An entity
pub trait Entity {
    type Grid: Grid;

    /// Point Iterator
    type PointIterator<'a>: Iterator<Item = <Self::Grid as Grid>::Point<'a>>
    where
        Self: 'a;

    /// Iterator over sub-entities
    type EntityIter<'a>: Iterator<Item = Self>
    where
        Self: 'a;

    /// Iterator over connected entities
    type ConnectedEntityIter<'a>: Iterator<Item = Self>
    where
        Self: 'a;

    /// The entity type (eg triangle, quadrilateral) of this entity
    fn entity_type(&self) -> <Self::Grid as Grid>::EntityDescriptor;

    /// A sub-entity of this entity
    fn sub_entity(&self, dim: usize, index: usize) -> Self;

    /// Iterator over sub-entities
    fn sub_entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// Returns connected entities that are either on the same process or are ghost entities.
    fn connected_entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// The local index of this entity on the current process
    fn local_index(&self) -> usize;

    /// The global index of this entity on the current process
    fn global_index(&self) -> usize;

    /// The geometry of this entity
    fn geometry(
        &self,
        eval_points: &[<Self::Grid as Grid>::T],
    ) -> <Self::Grid as Grid>::Geometry<'_>;

    /// The ownership of this entity
    fn ownership(&self) -> Ownership;
}

/// An entity with an associated ID
pub trait EntityId: Entity {
    /// The id of this entity
    fn id(&self) -> usize;
}
