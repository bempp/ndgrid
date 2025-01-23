//! Traits for a mesh entity
use super::{Entity, GeometryMap};
use crate::types::{Ownership, RealScalar};
use bempp_distributed_tools::IndexLayout;
use mpi::traits::Communicator;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::Iterator;
use std::rc::Rc;

/// A grid
pub trait Grid {
    /// Scalar type
    type T: RealScalar;

    /// Type used as identifier of different entity types
    type Entity<'a>: Entity<EntityDescriptor = Self::EntityDescriptor, T = Self::T>
    where
        Self: 'a;

    /// Geometry map type
    type GeometryMap<'a>: GeometryMap<T = Self::T>
    where
        Self: 'a;

    /// Type used as identifier of different entity types
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Iterator over sub-entities
    type EntityIter<'a>: Iterator<Item = Self::Entity<'a>>
    where
        Self: 'a;

    /// Dimension of the geometry of this grid
    fn geometry_dim(&self) -> usize;

    /// Dimension of the topology of this grid
    fn topology_dim(&self) -> usize;

    /// An entity in this grid
    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>>;

    /// The entity types of topological dimension `dim` contained in this grid
    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor];

    /// Number of entities
    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize;

    /// Number of cells
    fn cell_count(&self) -> usize {
        self.entity_types(self.topology_dim())
            .iter()
            .map(|&t| self.entity_count(t))
            .sum()
    }

    /// Owned cell count
    ///
    /// Note. The default implementation iterates through all grid to count the number of owned elements.
    /// Override this method if a more efficient implementation is available.
    fn owned_cell_count(&self) -> usize {
        self.entity_iter(self.topology_dim())
            .filter(|e| matches!(e.ownership(), Ownership::Owned))
            .count()
    }

    /// Iterator over entities
    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_>;

    /// An entity in this grid from an insertion id
    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>>;

    /// Geometry map from reference entity to physical entities at the given points
    ///
    /// `points` should have space [entity_topology_dim, npts] and use column-major ordering
    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_>;
}

pub trait ParallelGrid {
    //! MPI parallel grid

    /// Local grid type
    type LocalGrid: Grid;

    /// Communicator
    type C: Communicator;

    /// MPI communicator
    fn comm(&self) -> &Self::C;
    /// Local grid on the current process
    fn local_grid(&self) -> &Self::LocalGrid;

    /// Return the cell index layout that describes where each global cell lives
    fn cell_index_layout(&self) -> Rc<IndexLayout<'_, Self::C>>;
}
