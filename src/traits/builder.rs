//! Grid builder
use crate::{traits::Grid, types::Scalar};
#[cfg(feature = "mpi")]
use crate::{traits::ParallelGrid, types::GraphPartitioner};
#[cfg(feature = "mpi")]
use mpi::traits::Communicator;
use std::fmt::Debug;
use std::hash::Hash;

/// A builder is a factory that creates meshes.
///
/// After instantiation points and cells can be added.
/// To build the actual grid call [Builder::create_grid].
pub trait Builder {
    /// Type used as identifier of different entity types
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// The type of the grid that the builder creates
    type Grid: Grid<EntityDescriptor = Self::EntityDescriptor>;
    /// The floating point type used for coordinates
    type T: Scalar;
    /// The type of the data that is input to add a cell
    type CellData<'a>;

    /// Add a point to the grid
    fn add_point(&mut self, id: usize, data: &[Self::T]);

    /// Add a cell to the grid
    fn add_cell(&mut self, id: usize, cell_data: Self::CellData<'_>);

    /// Add a cell to the grid
    fn add_cell_from_nodes_and_type(
        &mut self,
        id: usize,
        nodes: &[usize],
        cell_type: Self::EntityDescriptor,
        cell_degree: usize,
    );

    /// Create the grid
    fn create_grid(&self) -> Self::Grid;

    /// Number of points
    fn point_count(&self) -> usize;

    /// Number of cells
    fn cell_count(&self) -> usize;

    /// Get the insertion ids of each point
    fn point_indices_to_ids(&self) -> &[usize];

    /// Get the insertion ids of each cell
    fn cell_indices_to_ids(&self) -> &[usize];

    /// Get the indices of the points of a cell
    fn cell_points(&self, index: usize) -> &[usize];

    /// Get the indices of the points of a cell
    fn cell_vertices(&self, index: usize) -> &[usize];

    /// Get the coordinates of a point
    fn point(&self, index: usize) -> &[Self::T];

    /// Get all points
    fn points(&self) -> &[Self::T];

    /// Get the type of a cell
    fn cell_type(&self, index: usize) -> Self::EntityDescriptor;

    /// Get the degree of a cell's geometry
    fn cell_degree(&self, index: usize) -> usize;

    /// Geometric dimension
    fn gdim(&self) -> usize;

    /// Topoligical dimension
    fn tdim(&self) -> usize;

    /// Number of points in a cell with the given type and degree
    fn npts(&self, cell_type: Self::EntityDescriptor, degree: usize) -> usize;
}

/// Trait for building a geometry
///
/// This trait is usually not called by the user. It provides
/// an interface to building the geometry information of the grid.
pub(crate) trait GeometryBuilder: Builder {
    /// Grid geometry type
    type GridGeometry;

    /// Create geometry
    fn create_geometry(
        &self,
        point_ids: &[usize],
        coordinates: &[Self::T],
        cell_points: &[usize],
        cell_types: &[Self::EntityDescriptor],
        cell_degrees: &[usize],
    ) -> Self::GridGeometry;
}

/// Trait for building a topology
///
/// This trait is usually not called by the user. It provides
/// an interface to building the topology information of the grid.
pub(crate) trait TopologyBuilder: Builder {
    /// Grid topology type
    type GridTopology;

    /// Create topology
    fn create_topology(
        &self,
        vertex_ids: Vec<usize>,
        cell_ids: Vec<usize>,
        cells: &[usize],
        cell_types: &[Self::EntityDescriptor],
    ) -> Self::GridTopology;

    /// Extract the cell vertices from the cell points
    fn extract_vertices(
        &self,
        cell_points: &[usize],
        cell_types: &[Self::EntityDescriptor],
        cell_degrees: &[usize],
    ) -> Vec<usize>;
}

/// Trait for building a grid from topology and geometry
///
/// This trait is usually not called by the user. It provides
/// an interface to building the grid from a given topology and Geometry.
pub(crate) trait GridBuilder: Builder + GeometryBuilder + TopologyBuilder {
    /// Create topology
    fn create_grid_from_topology_geometry(
        &self,
        topology: <Self as TopologyBuilder>::GridTopology,
        geometry: <Self as GeometryBuilder>::GridGeometry,
    ) -> <Self as Builder>::Grid;
}

/// MPI parallelized grid builder
#[cfg(feature = "mpi")]
pub trait ParallelBuilder: Builder {
    /// Parallel grid type
    type ParallelGrid<'a, C: Communicator + 'a>: ParallelGrid<C = C>
    where
        Self: 'a;

    /// Create a parallel grid (call from root)
    fn create_parallel_grid_root<'a, C: Communicator>(
        &self,
        comm: &'a C,
        partitioner: GraphPartitioner,
    ) -> Self::ParallelGrid<'a, C>;

    /// Create a parallel grid (call from other processes)
    fn create_parallel_grid<'a, C: Communicator>(
        &self,
        comm: &'a C,
        root_rank: i32,
    ) -> Self::ParallelGrid<'a, C>;
}
