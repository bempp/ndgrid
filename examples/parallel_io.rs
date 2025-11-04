use itertools::izip;
use mpi::{collective::CommunicatorCollectives, environment::Universe, traits::Communicator};
use ndelement::{ciarlet::CiarletElement, map::IdentityMap, types::ReferenceCellType};
use ndgrid::traits::{
    DistributableGrid, Entity, Grid, ParallelBuilder, ParallelGrid, RONExportParallel,
    RONImportParallel, Topology,
};
use ndgrid::{
    ParallelGridImpl, SingleElementGrid, SingleElementGridBuilder, shapes, types::GraphPartitioner,
};

/// Grid I/O
///
/// Demonstration of importing and exporting a grid in parallel
///
/// Serial I/O is demonstrated in the example `io.rs`
fn main() {
    let universe: Universe = mpi::initialize().unwrap();
    let comm = universe.world();
    let rank = comm.rank();

    let g = if rank == 0 {
        // Create a grid using the shapes module: unit_cube_boundary will mesh the surface of a cube
        let serial_g = shapes::unit_cube_boundary::<f64>(4, 5, 4, ReferenceCellType::Triangle);

        // Distribute this grid across processes
        serial_g.distribute(&comm, GraphPartitioner::None)
    } else {
        let b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
        b.create_parallel_grid(&comm, 0)
    };

    // If the serde option is used, the raw grid data can be exported in RON format
    g.export_as_ron("_unit_cube_boundary_parallel.ron");

    // Wait for export to finish
    comm.barrier();

    // A grid can be re-imported from raw RON data. Note that it must be imported on the same number of processes as it was exported using
    let g2 = ParallelGridImpl::<'_, _, SingleElementGrid::<f64, CiarletElement<f64, IdentityMap>>>::import_from_ron(&comm, "_unit_cube_boundary_parallel.ron");

    // Print the first 5 cells of each grid on process 0
    if rank == 0 {
        println!("The first 5 cells of the grids");
        for (cell, cell2) in izip!(
            g.local_grid().entity_iter(ReferenceCellType::Triangle),
            g2.local_grid().entity_iter(ReferenceCellType::Triangle)
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
    }
}
