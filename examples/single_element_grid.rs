use ndelement::types::ReferenceCellType;
use ndgrid::grid::serial::SingleElementGridBuilder;
use ndgrid::traits::{Builder, Entity, Geometry, Grid, Point, Topology};

/// Creating a single element grid
///
/// In a single element grid, the same finite element will be used to represent the geometry
/// of each cell. For example, a grid of bilinear quadrilaterals can be created by using a degree 1
/// element on a quadrilateral
fn main() {
    // When creating the grid builder, we give the physical/geometric dimension (3) and the cell type
    // and degree of the element
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
    // Add six points with ids 0 to 5
    b.add_point(0, &[0.0, 0.0, 0.0]);
    b.add_point(1, &[1.0, 0.0, 0.0]);
    b.add_point(2, &[2.0, 0.0, 0.2]);
    b.add_point(3, &[0.0, 1.0, 0.0]);
    b.add_point(4, &[1.0, 1.0, -0.2]);
    b.add_point(5, &[2.0, 1.0, 0.0]);
    // Add two cells
    b.add_cell(0, &[0, 1, 3, 4]);
    b.add_cell(1, &[1, 2, 4, 5]);
    // Create the grid
    let grid = b.create_grid();

    // Print the coordinates or each point in the mesh
    let mut coords = vec![0.0; grid.geometry_dim()];
    for point in grid.entity_iter(0) {
        point.geometry().points().collect::<Vec<_>>()[0].coords(coords.as_mut_slice());
        println!("point {}: {:#?}", point.local_index(), coords);
    }

    // Print the vertices of each cell
    for cell in grid.entity_iter(2) {
        println!(
            "cell {}: {:?} ",
            cell.local_index(),
            cell.topology().sub_entity_iter(0).collect::<Vec<_>>()
        );
    }
}
