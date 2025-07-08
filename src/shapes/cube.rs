//! Cube grids

use crate::{
    grid::local_grid::{SingleElementGrid, SingleElementGridBuilder},
    traits::Builder,
    types::RealScalar,
};
use ndelement::{ciarlet::CiarletElement, map::IdentityMap, types::ReferenceCellType};

/// Create a unit interval grid
///
/// The unit interval is the interval between (0,) and (1,)
pub fn unit_interval<T: RealScalar>(
    nx: usize,
) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
    let mut b = SingleElementGridBuilder::new_with_capacity(
        1,
        nx + 1,
        nx,
        (ReferenceCellType::Interval, 1),
    );
    for i in 0..nx + 1 {
        b.add_point(i, &[T::from(i).unwrap() / T::from(nx).unwrap()]);
    }

    for i in 0..nx {
        b.add_cell(i, &[i, i + 1]);
    }

    b.create_grid()
}

/// Create a unit square grid
///
/// The unit square is the square with corners at (0,0), (1,0), (0,1) and (1,1)
pub fn unit_square<T: RealScalar>(
    nx: usize,
    ny: usize,
    cell_type: ReferenceCellType,
) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
    let mut b = SingleElementGridBuilder::new_with_capacity(
        2,
        (nx + 1) * (ny + 1),
        match cell_type {
            ReferenceCellType::Triangle => 2 * nx * ny,
            ReferenceCellType::Quadrilateral => 2 * nx * ny,
            _ => {
                panic!("Unsupported cell type: {cell_type:?}")
            }
        },
        (cell_type, 1),
    );
    for i in 0..nx + 1 {
        for j in 0..ny + 1 {
            b.add_point(
                i * (ny + 1) + j,
                &[
                    T::from(i).unwrap() / T::from(nx).unwrap(),
                    T::from(j).unwrap() / T::from(ny).unwrap(),
                ],
            );
        }
    }

    match cell_type {
        ReferenceCellType::Triangle => {
            for i in 0..nx {
                for j in 0..ny {
                    let dx = ny + 1;
                    let dy = 1;
                    let origin = i * dx + j * dy;
                    b.add_cell(2 * (j * nx + i), &[origin, origin + dx, origin + dx + dy]);
                    b.add_cell(
                        2 * (j * nx + i) + 1,
                        &[origin, origin + dx + dy, origin + dy],
                    );
                }
            }
        }
        ReferenceCellType::Quadrilateral => {
            for i in 0..nx {
                for j in 0..ny {
                    let dx = ny + 1;
                    let dy = 1;
                    let origin = i * dx + j * dy;
                    b.add_cell(
                        j * nx + i,
                        &[origin, origin + dx, origin + dy, origin + dx + dy],
                    );
                }
            }
        }
        _ => {
            panic!("Unsupported cell type: {cell_type:?}")
        }
    }

    b.create_grid()
}

/// Create a unit cube grid
///
/// The unit cube is the cube with corners at (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1),
/// (1,0,1), (0,1,1) and (1,1,1)
pub fn unit_cube<T: RealScalar>(
    nx: usize,
    ny: usize,
    nz: usize,
    cell_type: ReferenceCellType,
) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
    let mut b = SingleElementGridBuilder::new_with_capacity(
        3,
        (nx + 1) * (ny + 1) * (nz + 1),
        match cell_type {
            ReferenceCellType::Tetrahedron => 6 * nx * ny * nz,
            ReferenceCellType::Hexahedron => 2 * nx * ny * nz,
            _ => {
                panic!("Unsupported cell type: {cell_type:?}")
            }
        },
        (cell_type, 1),
    );
    for i in 0..nx + 1 {
        for j in 0..ny + 1 {
            for k in 0..nz + 1 {
                b.add_point(
                    (i * (ny + 1) + j) * (nz + 1) + k,
                    &[
                        T::from(i).unwrap() / T::from(nx).unwrap(),
                        T::from(j).unwrap() / T::from(ny).unwrap(),
                        T::from(k).unwrap() / T::from(nz).unwrap(),
                    ],
                );
            }
        }
    }

    match cell_type {
        ReferenceCellType::Tetrahedron => {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let dx = (ny + 1) * (nz + 1);
                        let dy = nz + 1;
                        let dz = 1;
                        let origin = i * dx + j * dy + k * dz;
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k),
                            &[origin, origin + dx, origin + dx + dy, origin + dx + dy + dz],
                        );
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k) + 1,
                            &[origin, origin + dy, origin + dx + dy, origin + dx + dy + dz],
                        );
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k) + 2,
                            &[origin, origin + dx, origin + dx + dz, origin + dx + dy + dz],
                        );
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k) + 3,
                            &[origin, origin + dz, origin + dx + dz, origin + dx + dy + dz],
                        );
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k) + 4,
                            &[origin, origin + dy, origin + dy + dz, origin + dx + dy + dz],
                        );
                        b.add_cell(
                            6 * ((j * nx + i) * ny + k) + 5,
                            &[origin, origin + dz, origin + dy + dz, origin + dx + dy + dz],
                        );
                    }
                }
            }
        }
        ReferenceCellType::Hexahedron => {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let dx = (ny + 1) * (nz + 1);
                        let dy = nz + 1;
                        let dz = 1;
                        let origin = i * dx + j * dy + k * dz;
                        b.add_cell(
                            (j * nx + i) * ny + k,
                            &[
                                origin,
                                origin + dx,
                                origin + dy,
                                origin + dy,
                                origin + dz,
                                origin + dx + dz,
                                origin + dy + dz,
                                origin + dy + dz,
                            ],
                        );
                    }
                }
            }
        }
        _ => {
            panic!("Unsupported cell type: {cell_type:?}")
        }
    }

    b.create_grid()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_unit_interval() {
        let _g = unit_interval::<f64>(1);
        let _g = unit_interval::<f64>(2);
        let _g = unit_interval::<f64>(4);
        let _g = unit_interval::<f64>(7);
    }

    #[test]
    fn test_unit_square_triangle() {
        let _g = unit_square::<f64>(1, 1, ReferenceCellType::Triangle);
        let _g = unit_square::<f64>(1, 1, ReferenceCellType::Triangle);
        let _g = unit_square::<f64>(4, 5, ReferenceCellType::Triangle);
        let _g = unit_square::<f64>(7, 6, ReferenceCellType::Triangle);
    }
    #[test]
    fn test_unit_square_quadrilateral() {
        let _g = unit_square::<f64>(1, 1, ReferenceCellType::Quadrilateral);
        let _g = unit_square::<f64>(2, 2, ReferenceCellType::Quadrilateral);
        let _g = unit_square::<f64>(4, 5, ReferenceCellType::Quadrilateral);
        let _g = unit_square::<f64>(7, 6, ReferenceCellType::Quadrilateral);
    }

    #[test]
    fn test_unit_cube_tetrahedron() {
        let _g = unit_cube::<f64>(1, 1, 1, ReferenceCellType::Tetrahedron);
        let _g = unit_cube::<f64>(2, 2, 2, ReferenceCellType::Tetrahedron);
        let _g = unit_cube::<f64>(4, 5, 5, ReferenceCellType::Tetrahedron);
        let _g = unit_cube::<f64>(7, 6, 4, ReferenceCellType::Tetrahedron);
    }
    #[test]
    fn test_unit_cube_hexahedron() {
        let _g = unit_cube::<f64>(1, 1, 1, ReferenceCellType::Hexahedron);
        let _g = unit_cube::<f64>(2, 2, 2, ReferenceCellType::Hexahedron);
        let _g = unit_cube::<f64>(4, 5, 5, ReferenceCellType::Hexahedron);
        let _g = unit_cube::<f64>(7, 6, 4, ReferenceCellType::Hexahedron);
    }
}
