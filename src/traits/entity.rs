//! Traits for a mesh entity
use super::Grid;
use crate::types::Ownership;

/// An entity
pub trait Entity {
    /// Grid type
    type Grid: Grid;

    /// The entity type (eg triangle, quadrilateral) of this entity
    fn entity_type(&self) -> <Self::Grid as Grid>::EntityDescriptor;

    /// The local index of this entity on the current process
    fn local_index(&self) -> usize;

    /// The global index of this entity on the current process
    fn global_index(&self) -> usize;

    /// The geometry of this entity
    fn geometry(&self) -> <Self::Grid as Grid>::Geometry<'_>;

    /// The topology of this entity
    fn topology(&self) -> <Self::Grid as Grid>::Topology<'_>;

    /// The ownership of this entity
    fn ownership(&self) -> Ownership;
}

/// An entity with an associated ID
pub trait EntityId: Entity {
    /// The id of this entity
    fn id(&self) -> usize;
}
