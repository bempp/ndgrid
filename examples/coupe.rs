//? run --features "mpi"

#[cfg(feature = "mpi")]
use mpi::{environment::Universe, topology::Communicator};
#[cfg(feature = "mpi")]
use ndelement::types::ReferenceCellType;
#[cfg(feature = "mpi")]
use ndgrid::{
    grid::serial::SingleElementGridBuilder,
    traits::{Builder, Grid, ParallelBuilder},
};

#[cfg(feature = "mpi")]
fn main() {
    let n = 8;

    let mut b = SingleElementGridBuilder::<f64>::new(2, (ReferenceCellType::Quadrilateral, 1));

    let mut i = 0;
    for y in 0..n {
        for x in 0..n {
            b.add_point(i, &[x as f64 / (n - 1) as f64, y as f64 / (n - 1) as f64]);
            i += 1;
        }
    }

    let mut i = 0;
    for y in 0..n - 1 {
        for x in 0..n - 1 {
            let sw = n * y + x;
            b.add_cell(i, &[sw, sw + 1, sw + n, sw + n + 1]);
            i += 1;
        }
    }

    let universe: Universe = mpi::initialize().unwrap();
    let comm = universe.world();
    let rank = comm.rank();
    let grid = if rank == 0 {
        b.create_parallel_grid(&comm)
    } else {
        b.receive_parallel_grid(&comm, 0)
    };

    println!(
        "MPI rank {rank} has {} vertices",
        grid.entity_count(ReferenceCellType::Point)
    );
    println!(
        "MPI rank {rank} has {} cells",
        grid.entity_count(ReferenceCellType::Quadrilateral)
    );
}

#[cfg(not(feature = "mpi"))]
fn main() {}
