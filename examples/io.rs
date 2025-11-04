use itertools::izip;
use ndelement::{ciarlet::CiarletElement, map::IdentityMap, types::ReferenceCellType};
use ndgrid::traits::{
    Builder, Entity, GmshExport, GmshImport, Grid, RONExport, RONImport, Topology,
};
use ndgrid::{SingleElementGrid, SingleElementGridBuilder, shapes};

/// Grid I/O
///
/// Demonstration of importing and exporting a grid in serial
///
/// Parallel I/O is demonstrated in the example `parallel_io.rs`
fn main() {
    // Create a grid using the shapes module: unit_cube_boundary will mesh the surface of a cube
    let g = shapes::unit_cube_boundary::<f64>(4, 5, 4, ReferenceCellType::Triangle);

    // If the serde option is used, the raw grid data can be exported in RON format
    g.export_as_ron("_unit_cube_boundary.ron");

    // A grid can be re-imported from raw RON data
    let g2 = SingleElementGrid::<f64, CiarletElement<f64, IdentityMap>>::import_from_ron(
        "_unit_cube_boundary.ron",
    );

    // Print the first 5 cells of each grid
    println!("The first 5 cells of the grids");
    for (cell, cell2) in izip!(
        g.entity_iter(ReferenceCellType::Triangle),
        g2.entity_iter(ReferenceCellType::Triangle)
    )
    .take(5)
    {
        println!(
            "{:?} {:?}",
            cell.topology()
                .sub_entity_iter(ReferenceCellType::Point)
                .collect::<Vec<_>>(),
            cell2
                .topology()
                .sub_entity_iter(ReferenceCellType::Point)
                .collect::<Vec<_>>(),
        );
    }

    println!();

    // Alternatively, grids can be exported and imported to/from Gmsh files

    // Export the grid as a Gmsh .msh file
    g.export_as_gmsh("_unit_cube_boundary.msh");

    // To import from a Gmsh .msh file, a builder is used
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh("_unit_cube_boundary.msh");
    let g3 = b.create_grid();

    // Print the first 5 cells of each grid
    println!("The first 5 cells of the grids");
    for (cell, cell3) in izip!(
        g.entity_iter(ReferenceCellType::Triangle),
        g3.entity_iter(ReferenceCellType::Triangle)
    )
    .take(5)
    {
        println!(
            "{:?} {:?}",
            cell.topology()
                .sub_entity_iter(ReferenceCellType::Point)
                .collect::<Vec<_>>(),
            cell3
                .topology()
                .sub_entity_iter(ReferenceCellType::Point)
                .collect::<Vec<_>>(),
        );
    }
}
