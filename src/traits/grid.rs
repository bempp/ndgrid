//! Traits for a mesh entity
use super::{Entity, Geometry, Point, Topology};
use crate::types::RealScalar;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::Iterator;

/// A grid
pub trait Grid {
    /// World dimension
    const GDIM: usize;

    /// Grid dimension
    const TDIM: usize;

    /// Scalar type
    type T: RealScalar;

    /// Point Type
    type Point<'a>: Point<T = Self::T>
    where
        Self: 'a;

    /// Type used as identifier of different entity types
    type Entity<'a>: Entity<
        EntityDescriptor = Self::EntityDescriptor,
        Topology<'a> = Self::Topology<'a>,
        Geometry<'a> = Self::Geometry<'a>,
    >
    where
        Self: 'a;

    /// Topology type
    type Topology<'a>: Topology
    where
        Self: 'a;

    /// Geometry type
    type Geometry<'a>: Geometry<Point<'a> = Self::Point<'a>>
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

    /// Dimension of the geometry of this grid
    fn geometry_dim() -> usize {
        Self::GDIM
    }

    /// Dimension of the topology of this grid
    fn topology_dim() -> usize {
        Self::TDIM
    }

    /// An entity in this grid
    fn entity(&self, dim: usize, local_index: usize) -> Self::Entity<'_>;

    /// Iterator over entities
    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// An entity in this grid from an insertion id
    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>>;
}
