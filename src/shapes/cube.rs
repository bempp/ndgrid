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
    use crate::traits::{Entity, Geometry, Grid, Point};
    use approx::*;

    fn min(values: &[f64]) -> f64 {
        let mut out = values[0];
        for i in &values[1..] {
            if *i < out {
                out = *i;
            }
        }
        out
    }
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
        for cell in grid.entity_iter(grid.topology_dim()) {
            let g = cell.geometry();
            match cell.entity_type() {
                ReferenceCellType::Interval => {
                    let mut points = vec![0.0; 2];
                    for (i, p) in g.points().enumerate() {
                        p.coords(&mut points[i..i + 1]);
                    }
                    volume += points[1] - points[0];
                }
                ReferenceCellType::Triangle => {
                    let mut points = vec![0.0; 6];
                    for (i, p) in g.points().enumerate() {
                        p.coords(&mut points[2 * i..2 * i + 2]);
                    }
                    let x0 = min(&[points[0], points[2], points[4]]);
                    let y0 = min(&[points[1], points[3], points[5]]);
                    let x1 = max(&[points[0], points[2], points[4]]);
                    let y1 = max(&[points[1], points[3], points[5]]);
                    volume += (x1 - x0) * (y1 - y0) / 2.0;
                }
                ReferenceCellType::Quadrilateral => {
                    let mut points = vec![0.0; 8];
                    for (i, p) in g.points().enumerate() {
                        p.coords(&mut points[2 * i..2 * i + 2]);
                    }
                    let x0 = points[0];
                    let y0 = points[1];
                    let x1 = points[2];
                    let y1 = points[5];
                    volume += (x1 - x0) * (y1 - y0);
                }
                ReferenceCellType::Tetrahedron => {
                    let mut points = vec![0.0; 24];
                    for (i, p) in g.points().enumerate() {
                        p.coords(&mut points[3 * i..3 * i + 3]);
                    }
                    let x0 = min(&[points[0], points[3], points[6], points[9]]);
                    let y0 = min(&[points[1], points[4], points[7], points[10]]);
                    let z0 = min(&[points[2], points[5], points[8], points[11]]);
                    let x1 = max(&[points[0], points[3], points[6], points[9]]);
                    let y1 = max(&[points[1], points[4], points[7], points[10]]);
                    let z1 = max(&[points[2], points[5], points[8], points[11]]);
                    volume += (x1 - x0) * (y1 - y0) * (z1 - z0) / 6.0;
                }
                ReferenceCellType::Hexahedron => {
                    let mut points = vec![0.0; 24];
                    for (i, p) in g.points().enumerate() {
                        p.coords(&mut points[3 * i..3 * i + 3]);
                    }
                    let x0 = points[0];
                    let y0 = points[1];
                    let z0 = points[2];
                    let x1 = points[3];
                    let y1 = points[7];
                    let z1 = points[14];
                    volume += (x1 - x0) * (y1 - y0) * (z1 - z0);
                }
                _ => {
                    panic!("Unsupported cell");
                }
            }
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
