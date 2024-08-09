//? mpirun -n {{NPROCESSES}} --features "mpi"

#[cfg(feature = "mpi")]
use approx::assert_relative_eq;
#[cfg(feature = "mpi")]
use ndgrid::{
    grid::{
        parallel_grid::ParallelGrid,
        single_element_grid::{SingleElementGrid, SingleElementGridBuilder},
    },
    traits::{Builder, Geometry, Grid, ParallelBuilder, ParallelGrid as PGTrait, Point, Entity},
    types::Ownership,
};
#[cfg(feature = "mpi")]
use mpi::{
    environment::Universe,
    request::WaitGuard,
    traits::{Communicator, Destination, Source},
};
#[cfg(feature = "mpi")]
use ndelement::{
    ciarlet::LagrangeElementFamily,
    types::{Continuity, ReferenceCellType},
};
#[cfg(feature = "mpi")]
use rlst::{CsrMatrix, Shape};
#[cfg(feature = "mpi")]
use std::collections::HashMap;

#[cfg(feature = "mpi")]
fn example_single_element_grid<C: Communicator>(
    comm: &C,
    n: usize,
) -> ParallelGrid<'_, C, SingleElementGrid<f64>> {
    let rank = comm.rank();
    let size = comm.size();

    let mut b = SingleElementGridBuilder::<3, f64>::new((ReferenceCellType::Quadrilateral, 1));

    if rank == 0 {
        for y in 0..n {
            for x in 0..n {
                b.add_point(
                    y * n + x,
                    [x as f64 / (n - 1) as f64, y as f64 / (n - 1) as f64, 0.0],
                );
            }
        }

        for i in 0..n - 1 {
            for j in 0..n - 1 {
                b.add_cell(
                    i * (n - 1) + j,
                    vec![j * n + i, j * n + i + 1, j * n + i + n, j * n + i + n + 1],
                );
            }
        }

        let ncells = (n - 1).pow(2);

        let mut owners = HashMap::new();
        let mut c = 0;
        for r in 0..size {
            let end = if r + 1 == size {
                ncells
            } else {
                (r + 1) as usize * ncells / size as usize
            };
            while c < end {
                owners.insert(c, r as usize);
                c += 1;
            }
        }
        b.create_parallel_grid(comm, &owners)
    } else {
        b.receive_parallel_grid(comm, 0)
    }
}

#[cfg(feature = "mpi")]
fn test_parallel_single_element_grid<C: Communicator>(comm: &C) {
    let rank = comm.rank();
    let size = comm.size();

    let n = 10;
    let grid = example_single_element_grid(comm, n);

    let mut area = 0.0;
    for cell in grid.iter_all_cells() {
        if cell.ownership() == Ownership::Owned {
            area += cell.geometry().volume();
        }
    }
    if rank != 0 {
        mpi::request::scope(|scope| {
            let _sreq2 = WaitGuard::from(comm.process_at_rank(0).immediate_send(scope, &area));
        });
    } else {
        for p in 1..size {
            let (a, _status) = comm.process_at_rank(p).receive::<f64>();
            area += a;
        }
        assert_relative_eq!(area, 1.0, max_relative = 1e-10);
    }

    let mut nvertices = 0;
    for v in 0..grid.number_of_vertices() {
        if grid.vertex_from_index(v).ownership() == Ownership::Owned {
            nvertices += 1
        }
    }
    if rank != 0 {
        mpi::request::scope(|scope| {
            let _sreq2 = WaitGuard::from(comm.process_at_rank(0).immediate_send(scope, &nvertices));
        });
    } else {
        for p in 1..size {
            let (nv, _status) = comm.process_at_rank(p).receive::<usize>();
            nvertices += nv;
        }
        assert_eq!(nvertices, n * n);
    }
}

#[cfg(feature = "mpi")]
fn main() {
    let universe: Universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();

    if rank == 0 {
        println!("Testing SingleElementGrid in parallel.");
    }
    test_parallel_single_element_grid(&world);
}
#[cfg(not(feature = "mpi"))]
fn main() {}
