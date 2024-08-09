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

    /// Add a point to the grid
    fn add_point(&mut self, id: usize, data: &[Self::T]);

    /// Add a cell to the grid
    fn add_cell(&mut self, id: usize, cell_data: Self::CellData<'_>);

    /// Create the grid
    fn create_grid(self) -> Self::Grid;
}
