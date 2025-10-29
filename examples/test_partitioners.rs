use mpi::{
    collective::SystemOperation, environment::Universe, topology::Communicator,
    traits::CommunicatorCollectives,
};
use ndelement::types::ReferenceCellType;
use ndgrid::{
    grid::local_grid::SingleElementGridBuilder,
    traits::{Builder, Entity, Grid, ParallelBuilder, ParallelGrid},
    types::{GraphPartitioner, Ownership},
};

/// Test that a graph partitioner works
fn run_test<C: Communicator>(comm: &C, partitioner: GraphPartitioner) {
    let n = 10;

    let mut b = SingleElementGridBuilder::<f64>::new(2, (ReferenceCellType::Quadrilateral, 1));

    let rank = comm.rank();
    let grid = if rank == 0 {
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

        b.create_parallel_grid_root(comm, partitioner)
    } else {
        b.create_parallel_grid(comm, 0)
    };

    // Check that owned cells are sorted ahead of ghost cells

    let cell_count_owned = grid
        .local_grid()
        .entity_iter(ReferenceCellType::Quadrilateral)
        .filter(|entity| entity.is_owned())
        .count();

    // Now check that the first `cell_count_owned` entities are actually owned.
    for cell in grid
        .local_grid()
        .entity_iter(ReferenceCellType::Quadrilateral)
        .take(cell_count_owned)
    {
        assert!(cell.is_owned())
    }

    // Now make sure that the indices of the global cells are in consecutive order

    let mut cell_global_count = grid.cell_layout().local_range().0;

    for cell in grid
        .local_grid()
        .entity_iter(ReferenceCellType::Quadrilateral)
        .take(cell_count_owned)
    {
        assert_eq!(cell.global_index(), cell_global_count);
        cell_global_count += 1;
    }

    // Get the global indices.

    let global_vertices = grid
        .local_grid()
        .entity_iter(ReferenceCellType::Point)
        .filter(|e| matches!(e.ownership(), Ownership::Owned))
        .map(|e| e.global_index())
        .collect::<Vec<_>>();

    let nvertices = global_vertices.len();

    let global_cells = grid
        .local_grid()
        .entity_iter(ReferenceCellType::Quadrilateral)
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

/// Run tests
fn main() {
    let universe: Universe = mpi::initialize().unwrap();
    let comm = universe.world();

    if comm.rank() == 0 {
        println!("Testing GraphPartitioner::None");
    }
    run_test(&comm, GraphPartitioner::None);

    let mut p = vec![];
    for i in 0..81 {
        p.push(i % comm.size() as usize);
    }
    if comm.rank() == 0 {
        println!("Testing GraphPartitioner::Manual");
    }
    run_test(&comm, GraphPartitioner::Manual(p));

    if comm.rank() == 0 {
        println!("Testing GraphPartitioner::Coupe");
    }
    run_test(&comm, GraphPartitioner::Coupe);

    if comm.rank() == 0 {
        println!("Testing GraphPartitioner::Scotch");
    }
    run_test(&comm, GraphPartitioner::Scotch);
}
