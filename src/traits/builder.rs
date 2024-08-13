//! Grid builder
use crate::{traits::Grid, types::RealScalar};
use std::fmt::Debug;
use std::hash::Hash;

pub trait Builder {
    //! Object that can be used to build a mesh

    /// The type of the grid that the builder creates
    type Grid: Grid;
    /// The floating point type used for coordinates
    type T: RealScalar;
    /// The type of the data that is input to add a cell
    type CellData<'a>;
    /// Type used as identifier of different entity types
    type EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash;

    /// Add a point to the grid
    fn add_point(&mut self, id: usize, data: &[Self::T]);

    /// Add a cell to the grid
    fn add_cell(&mut self, id: usize, cell_data: Self::CellData<'_>);

    /// Create the grid
    fn create_grid(self) -> Self::Grid;

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

    /// Get the indices of the vertices of a cell
    fn cell_vertices(&self, index: usize) -> &[usize];

    /// Get the coordinates of a point
    fn point(&self, index: usize) -> &[Self::T];

    /// Get all points
    fn points(&self) -> &[Self::T];

    /// Get the type of a cell
    fn cell_type(&self, index: usize) -> Self::EntityDescriptor;

    /// Geometric dimension
    fn gdim(&self) -> usize;
}
