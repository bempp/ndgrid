//? run --features "mpi"

#[cfg(feature = "mpi")]
use ndelement::types::ReferenceCellType;
#[cfg(feature = "mpi")]
use ndgrid::{
    grid::serial::SingleElementGridBuilder,
    traits::{Builder, ParallelBuilder},
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
    b.test();
}

#[cfg(not(feature = "mpi"))]
fn main() {}
