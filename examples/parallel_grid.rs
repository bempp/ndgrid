//? mpirun -n {{NPROCESSES}} --features "serde"

use mpi::{
    collective::SystemOperation, environment::Universe, topology::Communicator,
    traits::CommunicatorCollectives,
};
use ndelement::types::ReferenceCellType;
use ndgrid::{
    grid::local_grid::SingleElementGridBuilder,
    traits::{Builder, Entity, Grid, ParallelBuilder},
    types::Ownership,
};

fn main() {
    let n = 1000;

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
        b.create_parallel_grid_root(&comm)
    } else {
        b.create_parallel_grid(&comm, 0)
    };

    // Get the global indices.

    let global_vertices = grid
        .entity_iter(0)
        .filter(|e| matches!(e.ownership(), Ownership::Owned))
        .map(|e| e.global_index())
        .collect::<Vec<_>>();

    let nvertices = global_vertices.len();

    let global_cells = grid
        .entity_iter(2)
        .filter(|e| matches!(e.ownership(), Ownership::Owned))
        .map(|e| e.global_index())
        .collect::<Vec<_>>();

    let ncells = global_cells.len();

    let mut total_cells: usize = 0;
    let mut total_vertices: usize = 0;

    comm.all_reduce_into(&ncells, &mut total_cells, SystemOperation::sum());
    comm.all_reduce_into(&nvertices, &mut total_vertices, SystemOperation::sum());

    assert_eq!(total_cells, (n - 1) * (n - 1));
    assert_eq!(total_vertices, n * n);
}
