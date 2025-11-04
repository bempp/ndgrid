use mpi::{collective::CommunicatorCollectives, environment::Universe, traits::Communicator};
use ndelement::{ciarlet::CiarletElement, map::IdentityMap, types::ReferenceCellType};
use ndgrid::{
    SingleElementGrid, SingleElementGridBuilder,
    grid::ParallelGridImpl,
    traits::{Builder, Grid, ParallelBuilder, RONExportParallel, RONImportParallel},
    types::GraphPartitioner,
};

fn create_single_element_grid_data(b: &mut SingleElementGridBuilder<f64>, n: usize) {
    for y in 0..n {
        for x in 0..n {
            b.add_point(
                y * n + x,
                &[x as f64 / (n - 1) as f64, y as f64 / (n - 1) as f64, 0.0],
            );
        }
    }

    for i in 0..n - 1 {
        for j in 0..n - 1 {
            b.add_cell(
                i * (n - 1) + j,
                &[j * n + i, j * n + i + 1, j * n + i + n, j * n + i + n + 1],
            );
        }
    }
}

fn example_single_element_grid<C: Communicator>(
    comm: &C,
    n: usize,
) -> ParallelGridImpl<'_, C, SingleElementGrid<f64, CiarletElement<f64, IdentityMap>>> {
    let rank = comm.rank();

    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));

    if rank == 0 {
        create_single_element_grid_data(&mut b, n);
        b.create_parallel_grid_root(comm, GraphPartitioner::None)
    } else {
        b.create_parallel_grid(comm, 0)
    }
}

/// Test that grids can be exported as RON in parallel
fn test_parallel_export<C: Communicator>(comm: &C) {
    let size = comm.size();

    let n = 10;
    let grid = example_single_element_grid(comm, n);
    let filename = format!("_examples_parallel_io_{size}ranks.ron");
    grid.export_as_ron(&filename);
}

/// Test that grids can be imported from RON in parallel
fn test_parallel_import<C: Communicator>(comm: &C) {
    use ndgrid::traits::ParallelGrid;

    let size = comm.size();

    let filename = format!("_examples_parallel_io_{size}ranks.ron");
    let grid =
        ParallelGridImpl::<'_, C, SingleElementGrid<f64, CiarletElement<f64, IdentityMap>>>::import_from_ron(
            comm, &filename,
        );

    let n = 10;
    let grid2 = example_single_element_grid(comm, n);

    assert_eq!(
        grid.local_grid().entity_count(ReferenceCellType::Point),
        grid2.local_grid().entity_count(ReferenceCellType::Point)
    );
}

/// Run tests
fn main() {
    let universe: Universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();

    if rank == 0 {
        println!("Testing parallel grid export");
    }
    test_parallel_export(&world);

    world.barrier();

    if rank == 0 {
        println!("Testing parallel grid import");
    }
    test_parallel_import(&world);
}
