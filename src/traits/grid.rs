//! Traits for a mesh entity
use super::{Entity, Geometry, Point, Topology};
use rlst::RlstScalar;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::Iterator;

/// A grid
pub trait Grid {
    /// World dimension
    const WORLD_DIM: usize;

    /// Grid dimension
    const GRID_DIM: usize;

    /// Scalar type
    type T: RlstScalar;

    /// Point Type
    type Point<'a>: Point
    where
        Self: 'a;

    /// Type used as identifier of different entity types
    type Entity<'a>: Entity
    where
        Self: 'a;

    /// Topology type
    type Topology<'a>: Topology
    where
        Self: 'a;

    /// Geometry type
    type Geometry<'a>: Geometry
    where
        Self: 'a;

    /// Type used as identifier of different entity types
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Iterator over sub-entities
    type EntityIter<'a>: Iterator<Item = Self::Entity<'a>>
    where
        Self: 'a;

    /// Iterator over points
    type PointIter<'a>: Iterator<Item = Self::Point<'a>>
    where
        Self: 'a;

    /// World dimension
    fn world_dim() -> usize {
        Self::WORLD_DIM
    }

    /// Grid dimension
    fn grid_dim() -> usize {
        Self::GRID_DIM
    }

    /// An entity in this grid
    fn entity(&self, dim: usize, local_index: usize) -> Self::Entity<'_>;

    /// Iterator over entities
    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;
}
