//! Traits for a mesh entity
use super::{Entity, Point};
use rlst::RlstScalar;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::Iterator;

/// The topology of a grid
pub trait GridTopology {
    /// Entity index iterator
    type EntityIndexIter<'a>: Iterator<Item = usize>
    where
        Self: 'a;

    /// Topological dimension
    fn dim(&self) -> usize;
    /// Iterator over sub-entity indices
    fn entity_indices(&self, dim: usize) -> Self::EntityIndexIter<'_>;
    /// Entity count
    fn entity_count(&self, dim: usize) -> usize;
}
/// The geometry of a grid
pub trait GridGeometry {
    /// Type used for coordinate values
    type T: RlstScalar;
    /// Point type
    type Point<'a>: Point<T = Self::T>
    where
        Self: 'a;
    /// Point iterator
    type PointIter<'a>: Iterator<Item = Self::Point<'a>>
    where
        Self: 'a;

    /// Geometric dimension
    fn dim(&self) -> usize;
    /// Points
    fn points(&self) -> Self::PointIter<'_>;
    /// Point count
    fn point_count(&self) -> usize;
}

/// A grid
pub trait Grid {
    /// Topology type
    type Topology<'a>: GridTopology
    where
        Self: 'a;
    /// Geometry type
    type Geometry<'a>: GridGeometry
    where
        Self: 'a;
    /// Type used as identifier of different entity types
    type EntityType: Debug + PartialEq + Eq + Clone + Copy + Hash;
    /// Type used as identifier of different entity types
    type Entity: Entity;
    /// Iterator over sub-entities
    type EntityIter<'a>: Iterator<Item = Self::Entity>
    where
        Self: 'a;

    /// An entity in this grid
    fn entity(&self, dim: usize, index: usize) -> Self::Entity;
    /// Iterator over entities
    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;
    /// The topology of this entity
    fn topology(&self) -> Self::Topology<'_>;
    /// The geometry of this entity
    fn geometry(&self) -> Self::Geometry<'_>;
}
