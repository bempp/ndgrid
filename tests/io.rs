//! Test input/output
use ndelement::types::ReferenceCellType;
use ndgrid::{
    shapes::regular_sphere,
    traits::{Builder, GmshExport},
    SingleElementGridBuilder,
};

#[test]
fn test_regular_sphere_gmsh_io() {
    let g = regular_sphere::<f64>(2);
    g.export_as_gmsh(String::from("_test_io_sphere.msh"));
}

#[test]
fn test_gmsh_output_quads() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
    b.add_point(0, &[0.0, 0.0, 0.0]);
    b.add_point(1, &[1.0, 0.0, 0.0]);
    b.add_point(2, &[0.0, 1.0, 0.0]);
    b.add_point(3, &[1.0, 1.0, 0.0]);
    b.add_point(4, &[0.0, 0.0, 1.0]);
    b.add_point(5, &[1.0, 0.0, 1.0]);
    b.add_point(6, &[0.0, 1.0, 1.0]);
    b.add_point(7, &[1.0, 1.0, 1.0]);
    b.add_cell(0, &[0, 2, 1, 3]);
    b.add_cell(1, &[0, 1, 4, 5]);
    b.add_cell(2, &[0, 4, 2, 6]);
    b.add_cell(3, &[1, 3, 5, 7]);
    b.add_cell(4, &[2, 6, 3, 7]);
    b.add_cell(5, &[4, 5, 6, 7]);
    let g = b.create_grid();
    g.export_as_gmsh(String::from("_test_io_cube.msh"));
}
