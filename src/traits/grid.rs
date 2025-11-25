//! Traits for a mesh entity
use super::{Entity, GeometryMap};
#[cfg(feature = "mpi")]
use crate::types::GraphPartitioner;
use crate::types::{Ownership, Scalar};
#[cfg(feature = "mpi")]
use mpi::traits::Communicator;
#[cfg(feature = "mpi")]
use rlst::distributed_tools::IndexLayout;
use std::fmt::Debug;
use std::hash::Hash;
use std::iter::Iterator;
#[cfg(feature = "mpi")]
use std::rc::Rc;

/// A grid provides access to entities, their geometrical and their topological properties.
pub trait Grid {
    /// Scalar type
    type T: Scalar;

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
    fn entity(
        &self,
        entity_type: Self::EntityDescriptor,
        local_index: usize,
    ) -> Option<Self::Entity<'_>>;

    /// The entity types of topological dimension `dim` contained in this grid
    fn entity_types(&self, tdim: usize) -> &[Self::EntityDescriptor];

    /// Number of entities of type `entity_type`
    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize;

    /// Number of cells in the grid
    fn cell_count(&self) -> usize {
        self.entity_types(self.topology_dim())
            .iter()
            .map(|&t| self.entity_count(t))
            .sum()
    }

    /// Return the cell types in the grid
    fn cell_types(&self) -> &[Self::EntityDescriptor] {
        let tdim = self.topology_dim();
        self.entity_types(tdim)
    }

    /// Owned cell count
    ///
    /// Note. The default implementation iterates through all grid to count the number of owned elements.
    /// Override this method if a more efficient implementation is available.
    fn owned_cell_count(&self) -> usize {
        self.cell_types()
            .iter()
            .map(|t| {
                self.entity_iter(*t)
                    .filter(|e| matches!(e.ownership(), Ownership::Owned))
                    .count()
            })
            .sum()
    }

    /// Iterator over entities
    fn entity_iter(&self, entity_type: Self::EntityDescriptor) -> Self::EntityIter<'_>;

    /// An entity in this grid from an insertion id
    fn entity_from_id(
        &self,
        entity_type: Self::EntityDescriptor,
        id: usize,
    ) -> Option<Self::Entity<'_>>;

    /// Geometry map from reference entity to physical entities at the given points
    ///
    /// `points` should have space [entity_topology_dim, npts] and use column-major ordering
    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        geometry_degree: usize,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_>;
}

/// Definition of an MPI parallel grid
#[cfg(feature = "mpi")]
pub trait ParallelGrid {
    /// Type of the Grid
    type T: Scalar;

    /// Local grid type
    type LocalGrid: Grid<T = Self::T>;

    /// Communicator
    type C: Communicator;

    /// MPI communicator
    fn comm(&self) -> &Self::C;
    /// Local grid on the current process
    fn local_grid(&self) -> &Self::LocalGrid;

    /// Return the cell index layout that describes where each global cell lives
    fn cell_layout(&self) -> Rc<IndexLayout<'_, Self::C>>;

    /// Return the global number of cells
    fn global_cell_count(&self) -> usize {
        self.cell_layout().number_of_global_indices()
    }
}

/// A grid that can be be distributed across processes
#[cfg(feature = "mpi")]
pub trait DistributableGrid {
    /// Parallel grid type when distrubuted
    type ParallelGrid<'a, C: Communicator + 'a>: ParallelGrid<C = C>;

    /// Distribute this grid in parallel
    fn distribute<'a, C: Communicator>(
        &self,
        comm: &'a C,
        partitioner: GraphPartitioner,
    ) -> Self::ParallelGrid<'a, C>;
}
