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

/// Create a grid of the boundary of a unit square
///
/// The unit square is the square with corners at (0,0), (1,0), (0,1) and (1,1)
pub fn unit_square_boundary<T: RealScalar>(
    nx: usize,
    ny: usize,
) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
    let mut b = SingleElementGridBuilder::new_with_capacity(
        2,
        2 * (nx + ny),
        2 * (nx + ny),
        (ReferenceCellType::Interval, 1),
    );
    let dx = ny + 1;
    let dy = 1;

    println!("Starting: {nx},{ny}");

    for i in 0..nx + 1 {
        println!("i * dx = {}", i * dx);
        b.add_point(i * dx, &[T::from(i).unwrap() / T::from(nx).unwrap(), T::zero()]);
        println!("i * dx + ny * dy = {}", i * dx + ny * dy);
        b.add_point(i * dx + ny * dy, &[T::from(i).unwrap() / T::from(nx).unwrap(), T::one()]);
    }
    for j in 1..ny {
        println!("j * dy = {}", j * dy);
        b.add_point(j * dy, &[T::zero(), T::from(j).unwrap() / T::from(ny).unwrap()]);
        println!("nx * dx + j * dy = {}", nx * dx + j * dy);
        b.add_point(nx * dx + j * dy, &[T::one(), T::from(j).unwrap() / T::from(ny).unwrap()]);
    }

    println!();

    for i in 0..nx {
        let origin = i * dx;
        println!("[{} {}]", origin, origin + dx);
        b.add_cell(2 * i, &[origin, origin + dx]);
        let origin = i * dx + ny * dy;
        println!("[{} {}]", origin, origin + dx);
        b.add_cell(2 * i + 1, &[origin, origin + dx]);
    }
    for j in 0..ny {
        let origin = j * dy;
        println!("[{} {}]", origin, origin + dy);
        b.add_cell(2 * nx + 2 * j, &[origin, origin + dy]);
        let origin = nx * dx + j * dy;
        println!("[{} {}]", origin, origin + dy);
        b.add_cell(2 * nx + 2 * j + 1, &[origin, origin + dy]);
    }

    println!();

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

/// Create a grid of the boundary of a unit cube
///
/// The unit cube is the cube with corners at (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1),
/// (1,0,1), (0,1,1) and (1,1,1)
pub fn unit_cube_boundary<T: RealScalar>(
    nx: usize,
    ny: usize,
    nz: usize,
    cell_type: ReferenceCellType,
) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
    unimplemented!();
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
    use crate::traits::{Entity, Geometry, Grid, Point};
    use approx::*;
    use itertools::izip;

    fn max(values: &[f64]) -> f64 {
        let mut out = values[0];
        for i in &values[1..] {
            if *i > out {
                out = *i;
            }
        }
        out
    }

    fn check_volume(
        grid: &impl Grid<T = f64, EntityDescriptor = ReferenceCellType>,
        expected_volume: f64,
    ) {
        let mut volume = 0.0;
        let tdim = grid.topology_dim();
        let gdim = grid.geometry_dim();
        for cell in grid.entity_iter(tdim) {
            let g = cell.geometry();
            let mut point = vec![0.0; gdim];
            let mut min_p = vec![10.0; gdim];
            let mut max_p = vec![-10.0; gdim];
            for p in g.points() {
                p.coords(&mut point);
                for (j, v) in point.iter().enumerate() {
                    if *v < min_p[j] {
                        min_p[j] = *v;
                    }
                    if *v > max_p[j] {
                        max_p[j] = *v;
                    }
                }
            }
            volume += match cell.entity_type() {
                ReferenceCellType::Interval => {
                    max(&izip!(min_p, max_p).map(|(i, j)| j - i).collect::<Vec<_>>())
                }
                ReferenceCellType::Triangle => match gdim {
                    2 => (max_p[0] - min_p[0]) * (max_p[1] - min_p[1]) / 2.0,
                    3 => max(&[
                        (max_p[0] - min_p[0]) * (max_p[1] - min_p[1]) / 2.0,
                        (max_p[0] - min_p[0]) * (max_p[2] - min_p[2]) / 2.0,
                        (max_p[1] - min_p[1]) * (max_p[2] - min_p[2]) / 2.0,
                    ]),
                    _ => {
                        panic!("Unsupported dimension");
                    }
                },
                ReferenceCellType::Quadrilateral => (max_p[0] - min_p[0]) * (max_p[1] - min_p[1]),
                ReferenceCellType::Tetrahedron => {
                    (max_p[0] - min_p[0]) * (max_p[1] - min_p[1]) * (max_p[2] - min_p[2]) / 6.0
                }
                ReferenceCellType::Hexahedron => {
                    (max_p[0] - min_p[0]) * (max_p[1] - min_p[1]) * (max_p[2] - min_p[2])
                }
                _ => {
                    panic!("Unsupported cell");
                }
            };
        }
        assert_relative_eq!(volume, expected_volume, epsilon = 1e-10);
    }

    #[test]
    fn test_unit_interval() {
        check_volume(&unit_interval::<f64>(1), 1.0);
        check_volume(&unit_interval::<f64>(2), 1.0);
        check_volume(&unit_interval::<f64>(4), 1.0);
        check_volume(&unit_interval::<f64>(7), 1.0);
    }

    #[test]
    fn test_unit_square_triangle() {
        check_volume(&unit_square::<f64>(1, 1, ReferenceCellType::Triangle), 1.0);
        check_volume(&unit_square::<f64>(1, 1, ReferenceCellType::Triangle), 1.0);
        check_volume(&unit_square::<f64>(4, 5, ReferenceCellType::Triangle), 1.0);
        check_volume(&unit_square::<f64>(7, 6, ReferenceCellType::Triangle), 1.0);
    }

    #[test]
    fn test_unit_square_quadrilateral() {
        check_volume(
            &unit_square::<f64>(1, 1, ReferenceCellType::Quadrilateral),
            1.0,
        );
        check_volume(
            &unit_square::<f64>(2, 2, ReferenceCellType::Quadrilateral),
            1.0,
        );
        check_volume(
            &unit_square::<f64>(4, 5, ReferenceCellType::Quadrilateral),
            1.0,
        );
        check_volume(
            &unit_square::<f64>(7, 6, ReferenceCellType::Quadrilateral),
            1.0,
        );
    }

    #[test]
    fn test_unit_square_boundary() {
        check_volume(&unit_square_boundary::<f64>(1, 1), 4.0);
        check_volume(&unit_square_boundary::<f64>(1, 1), 4.0);
        check_volume(&unit_square_boundary::<f64>(4, 5), 4.0);
        check_volume(&unit_square_boundary::<f64>(7, 6), 4.0);
    }

    #[test]
    fn test_unit_cube_tetrahedron() {
        check_volume(
            &unit_cube::<f64>(1, 1, 1, ReferenceCellType::Tetrahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(2, 2, 2, ReferenceCellType::Tetrahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(4, 5, 5, ReferenceCellType::Tetrahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(7, 6, 4, ReferenceCellType::Tetrahedron),
            1.0,
        );
    }
    #[test]
    fn test_unit_cube_hexahedron() {
        check_volume(
            &unit_cube::<f64>(1, 1, 1, ReferenceCellType::Hexahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(2, 2, 2, ReferenceCellType::Hexahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(4, 5, 5, ReferenceCellType::Hexahedron),
            1.0,
        );
        check_volume(
            &unit_cube::<f64>(7, 6, 4, ReferenceCellType::Hexahedron),
            1.0,
        );
    }
}
