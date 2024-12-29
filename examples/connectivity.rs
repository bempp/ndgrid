//! Experiments with Element connectivity.

pub fn main() {
    use ndelement::reference_cell;

    let connectivity = reference_cell::connectivity(ndelement::types::ReferenceCellType::Triangle);

    // The zero-dimensional entity with index 1 is connected to the following zero-dimensional entities.
    println!("Connectivity: {:#?}", connectivity[0][1][0]);
}
