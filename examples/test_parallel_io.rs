//? mpirun -n {{NPROCESSES}} --features "mpi,serde"

#[cfg(feature = "mpi")]
use mpi::{collective::CommunicatorCollectives, environment::Universe, traits::Communicator};
#[cfg(feature = "mpi")]
use ndelement::{ciarlet::CiarletElement, types::ReferenceCellType};
#[cfg(feature = "mpi")]
use ndgrid::{
    grid::parallel::ParallelGrid,
    traits::{Builder, Grid, ParallelBuilder, RONExportParallel, RONImportParallel},
    SingleElementGrid, SingleElementGridBuilder,
};

#[cfg(feature = "mpi")]
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

#[cfg(feature = "mpi")]
fn example_single_element_grid<C: Communicator>(
    comm: &C,
    n: usize,
) -> ParallelGrid<'_, C, SingleElementGrid<f64, CiarletElement<f64>>> {
    let rank = comm.rank();

    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));

    if rank == 0 {
        create_single_element_grid_data(&mut b, n);
        b.create_parallel_grid(comm)
    } else {
        b.receive_parallel_grid(comm, 0)
    }
}

#[cfg(feature = "mpi")]
fn test_parallel_export<C: Communicator>(comm: &C) {
    let size = comm.size();

    let n = 10;
    let grid = example_single_element_grid(comm, n);
    let filename = format!("_examples_parallel_io_{}ranks.ron", size);
    grid.export_as_ron(&filename);
}

#[cfg(feature = "mpi")]
fn test_parallel_import<C: Communicator>(comm: &C) {
    let size = comm.size();

    let filename = format!("_examples_parallel_io_{}ranks.ron", size);
    let grid = ParallelGrid::<'_, C, SingleElementGrid<f64, CiarletElement<f64>>>::import_from_ron(
        comm, &filename,
    );

    let n = 10;
    let grid2 = example_single_element_grid(comm, n);

    assert_eq!(
        grid.entity_count(ReferenceCellType::Point),
        grid2.entity_count(ReferenceCellType::Point)
    );
}

#[cfg(feature = "mpi")]
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
#[cfg(not(feature = "mpi"))]
fn main() {}
