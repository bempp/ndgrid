//! Traits for a mesh entity
use crate::traits::{Geometry, Topology};
use crate::types::{Ownership, Scalar};
use std::fmt::Debug;
use std::hash::Hash;

/// Definition of a grid entity
///
/// A grid entity can be a vertex, edge, face or cell. This trait
/// provides a unified interface to any type of entity.
pub trait Entity {
    /// Scalar type
    type T: Scalar;

    /// Type used as identifier of different entity types.
    /// In most cases this is given by [ReferenceCellType](ndelement::types::ReferenceCellType).
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Topology type
    type Topology<'a>: Topology<EntityDescriptor = Self::EntityDescriptor>
    where
        Self: 'a;

    /// Geometry type
    type Geometry<'a>: Geometry<T = Self::T>
    where
        Self: 'a;

    /// The entity type (eg triangle, quadrilateral) of this entity.
    fn entity_type(&self) -> Self::EntityDescriptor;

    /// The local index of this entity on the current process.
    fn local_index(&self) -> usize;

    /// The global index of this entity across all processes.
    fn global_index(&self) -> usize;

    /// The geometry of this entity.
    fn geometry(&self) -> Self::Geometry<'_>;

    /// The topology of this entity.
    fn topology(&self) -> Self::Topology<'_>;

    /// The ownership of this entity.
    fn ownership(&self) -> Ownership;

    /// Return true if the entity is owned.
    fn is_owned(&self) -> bool {
        matches!(self.ownership(), Ownership::Owned)
    }

    /// The insertion id of this entity.
    ///
    /// Return `None` if the entity has no insertion id.
    fn id(&self) -> Option<usize>;
}
