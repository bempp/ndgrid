//! Traits for a mesh entity
use crate::types::Ownership;
use std::iter::Iterator;
use rlst::RlstScalar;
use super::Point;

/// The topology of an entity
pub trait EntityTopology {
    /// Entity index iterator
    type EntityIndexIter<'a>: Iterator<Item=usize> where Self: 'a;

    /// Topological dimension
    fn dim(&self) -> usize;
    /// Iterator over sub-entity indices
    fn sub_entity_indices(&self, dim: usize) -> Self::EntityIndexIter<'_>;
}
/// The geometry of an entity
pub trait EntityGeometry {
    /// Type used for coordinate values
    type T: RlstScalar;
    /// Point type
    type Point<'a>: Point<T=Self::T> where Self: 'a;
    /// Point iterator
    type PointIter<'a>: Iterator<Item=Self::Point<'a>> where Self:'a;

    /// Geometric dimension
    fn dim(&self) -> usize;
    /// Volume of the entity
    fn volume(&self) -> Self::T;
    /// Diameter of the entity
    fn diameter(&self) -> Self::T;
    /// The points of the entity
    fn points(&self) -> Self::PointIter<'_>;
}

/// An entity
pub trait Entity {
    /// Topology type
    type Topology<'a>: EntityTopology where Self: 'a;
    /// Geometry type
    type Geometry<'a>: EntityGeometry where Self: 'a;
    /// Type used as identifier of different entity types
    type EntityType;
    /// Iterator over sub-entities
    type SubentityIter<'a>: Iterator<Item=Self> where Self: 'a;

    /// The entity type (eg triangle, quadrilateral) of this entity
    fn entity_type(&self) -> Self::EntityType;
    /// A sub-entity of this entity
    fn sub_entity(&self, dim: usize, index: usize) -> Self;
    /// Iterator over sub-entities
    fn sub_entity_iter(&self, dim: usize) -> Self::SubentityIter<'_>;
    /// The local index of this entity on the current process
    fn local_index(&self) -> usize;
    /// The global index of this entity on the current process
    fn global_index(&self) -> usize;
    /// The topology of this entity
    fn topology(&self) -> Self::Topology<'_>;
    /// The geometry of this entity
    fn geometry(&self) -> Self::Geometry<'_>;
    /// The ownership of this entity
    fn ownership(&self) -> Ownership;
}

/// An entity with an associated ID
pub trait EntityId: Entity {
    /// The id of this entity
    fn id(&self) -> usize;
}
