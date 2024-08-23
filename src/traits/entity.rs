//! Traits for a mesh entity
use crate::traits::{Geometry, Topology};
use crate::types::{Ownership, RealScalar};
use std::fmt::Debug;
use std::hash::Hash;

/// An entity
pub trait Entity {
    /// Scalar type
    type T: RealScalar;

    /// Type used as identifier of different entity types
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Topology type
    type Topology<'a>: Topology
    where
        Self: 'a;

    /// Geometry type
    type Geometry<'a>: Geometry<T = Self::T>
    where
        Self: 'a;

    /// The entity type (eg triangle, quadrilateral) of this entity
    fn entity_type(&self) -> Self::EntityDescriptor;

    /// The local index of this entity on the current process
    fn local_index(&self) -> usize;

    /// The global index of this entity on the current process
    fn global_index(&self) -> usize;

    /// The geometry of this entity
    fn geometry(&self) -> Self::Geometry<'_>;

    /// The topology of this entity
    fn topology(&self) -> Self::Topology<'_>;

    /// The ownership of this entity
    fn ownership(&self) -> Ownership;

    /// The insertion id of this entity
    fn id(&self) -> Option<usize>;
}
