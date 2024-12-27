use ndelement::types::ReferenceCellType;
use ndgrid::{
    traits::{Builder, Grid},
    SingleElementGridBuilder,
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

    println!("Setting up grid.");

    let grid = b.create_grid();

    println!(
        "Grid created with {} vertices and {} cells",
        grid.entity_count(ReferenceCellType::Point),
        grid.entity_count(ReferenceCellType::Quadrilateral)
    );
}
