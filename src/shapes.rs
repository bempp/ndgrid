//! Functions to create simple example grids

mod cube;
mod regular_sphere;
mod screen;

pub use cube::{
    unit_cube, unit_cube_boundary, unit_cube_edges, unit_interval, unit_square,
    unit_square_boundary,
};
pub use regular_sphere::regular_sphere;
pub use screen::{screen_quadrilaterals, screen_triangles};
