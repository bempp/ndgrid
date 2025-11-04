use mpi::{environment::Universe, traits::Communicator};
use ndelement::types::ReferenceCellType;
use ndgrid::{
    SingleElementGridBuilder,
    traits::{Builder, Grid, ParallelBuilder, ParallelGrid},
    types::GraphPartitioner,
};

/// Test that using non-contiguous numbering does not cause panic
fn test_noncontiguous_numbering<C: Communicator>(comm: &C) {
    let rank = comm.rank();
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));

    let g = if rank == 0 {
        let n = 5;
        for y in 0..n {
            for x in 0..n {
                b.add_point(
                    2 * (y * n + x) + 5,
                    &[x as f64 / (n - 1) as f64, y as f64 / (n - 1) as f64, 0.0],
                );
            }
        }

        for i in 0..n - 1 {
            for j in 0..n - 1 {
                b.add_cell(
                    3 * (i * (n - 1) + j),
                    &[
                        2 * (j * n + i) + 5,
                        2 * (j * n + i + 1) + 5,
                        2 * (j * n + i + n) + 5,
                        2 * (j * n + i + n + 1) + 5,
                    ],
                );
            }
        }

        b.create_parallel_grid_root(comm, GraphPartitioner::None)
    } else {
        b.create_parallel_grid(comm, 0)
    };

    assert!(g.local_grid().entity_count(ReferenceCellType::Point) > 0);
}

/// Run tests
fn main() {
    let universe: Universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();

    if rank == 0 {
        println!("Testing non-contiguous numbering");
    }
    test_noncontiguous_numbering(&world);
}
