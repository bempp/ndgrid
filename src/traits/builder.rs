//! Grid builder
use crate::{traits::Grid, types::RealScalar};

pub trait Builder {
    //! Object that can be used to build a mesh

    /// The type of the grid that the builder creates
    type Grid: Grid;
    /// The floating point type used for coordinates
    type T: RealScalar;
    /// The type of the data that is input to add a cell
    type CellData<'a>;
    /// The type of the data that must be provided when initialising the builder
    type GridMetadata;

    /// Create a new grid builder
    fn new(gdim: usize, data: Self::GridMetadata) -> Self;

    /// Create a new grid builder with capacaty for a given number of points and cells
    fn new_with_capacity(
        gdim: usize,
        npoints: usize,
        ncells: usize,
        data: Self::GridMetadata,
    ) -> Self;

    /// Add a point to the grid
    fn add_point(&mut self, id: usize, data: &[Self::T]);

    /// Add a cell to the grid
    fn add_cell(&mut self, id: usize, cell_data: Self::CellData<'_>);

    /// Create the grid
    fn create_grid(self) -> Self::Grid;
}
