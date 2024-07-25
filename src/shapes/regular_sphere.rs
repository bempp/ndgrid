//! Regular sphere grid

use crate::{
    grid::serial::{SingleElementGrid, SingleElementGridBuilder},
    traits::Builder,
    types::RealScalar,
};
use ndelement::{ciarlet::CiarletElement, types::ReferenceCellType};
use num::Float;
use std::collections::{hash_map::Entry::Vacant, HashMap};

/// Create a regular sphere
///
/// A regular sphere is created by starting with a regular octahedron. The shape is then refined `refinement_level` times.
/// Each time the grid is refined, each triangle is split into four triangles (by adding lines connecting the midpoints of
/// each edge). The new points are then scaled so that they are a distance of 1 from the origin.
pub fn regular_sphere<T: RealScalar>(
    refinement_level: u32,
) -> SingleElementGrid<T, CiarletElement<T>> {
    let mut b = SingleElementGridBuilder::new_with_capacity(
        3,
        2 + usize::pow(4, refinement_level + 1),
        8 * usize::pow(4, refinement_level),
        (ReferenceCellType::Triangle, 1),
    );
    let zero = T::from(0.0).unwrap();
    let one = T::from(1.0).unwrap();
    let half = T::from(0.5).unwrap();
    b.add_point(0, &[zero, zero, one]);
    b.add_point(1, &[one, zero, zero]);
    b.add_point(2, &[zero, one, zero]);
    b.add_point(3, &[-one, zero, zero]);
    b.add_point(4, &[zero, -one, zero]);
    b.add_point(5, &[zero, zero, -one]);
    let mut point_n = 6;

    let mut cells = vec![
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 1],
        [5, 2, 1],
        [5, 3, 2],
        [5, 4, 3],
        [5, 1, 4],
    ];
    let mut v = [[zero, zero, zero], [zero, zero, zero], [zero, zero, zero]];

    for level in 0..refinement_level {
        let mut edge_points = HashMap::new();
        let mut new_cells = Vec::with_capacity(8 * usize::pow(6, level));
        for c in &cells {
            for i in 0..3 {
                for j in 0..3 {
                    v[i][j] = b.points[3 * c[i] + j];
                }
            }
            let edges = [[1, 2], [0, 2], [0, 1]]
                .iter()
                .map(|[i, j]| {
                    let mut pt_i = c[*i];
                    let mut pt_j = c[*j];
                    if pt_i > pt_j {
                        std::mem::swap(&mut pt_i, &mut pt_j);
                    }
                    if let Vacant(e) = edge_points.entry((pt_i, pt_j)) {
                        let v_i = v[*i];
                        let v_j = v[*j];
                        let mut new_pt = [
                            half * (v_i[0] + v_j[0]),
                            half * (v_i[1] + v_j[1]),
                            half * (v_i[2] + v_j[2]),
                        ];
                        let size =
                            Float::sqrt(new_pt.iter().map(|x| Float::powi(*x, 2)).sum::<T>());
                        for i in new_pt.iter_mut() {
                            *i /= size;
                        }
                        b.add_point(point_n, &new_pt);
                        e.insert(point_n);
                        point_n += 1;
                    }
                    edge_points[&(pt_i, pt_j)]
                })
                .collect::<Vec<_>>();
            new_cells.push([c[0], edges[2], edges[1]]);
            new_cells.push([c[1], edges[0], edges[2]]);
            new_cells.push([c[2], edges[1], edges[0]]);
            new_cells.push([edges[0], edges[1], edges[2]]);
        }
        cells = new_cells;
    }
    for (i, v) in cells.iter().enumerate() {
        b.add_cell(i, v);
    }

    b.create_grid()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::traits::{GeometryMap, Grid};
    use approx::assert_relative_eq;

    #[test]
    fn test_regular_sphere_0() {
        let _g = regular_sphere::<f64>(0);
    }

    #[test]
    fn test_regular_spheres() {
        let _g1 = regular_sphere::<f64>(1);
        let _g2 = regular_sphere::<f64>(2);
        let _g3 = regular_sphere::<f64>(3);
    }

    #[test]
    fn test_normal_is_outward() {
        for i in 0..3 {
            let g = regular_sphere::<f64>(i);
            let points = vec![1.0 / 3.0, 1.0 / 3.0];
            let map = g.geometry_map(ReferenceCellType::Triangle, &points);
            let mut mapped_pt = vec![0.0; 3];
            let mut j = vec![0.0; 6];
            let mut jdet = vec![0.0];
            let mut normal = vec![0.0; 3];
            for i in 0..g.entity_count(ReferenceCellType::Triangle) {
                map.points(i, &mut mapped_pt);
                map.jacobians_dets_normals(i, &mut j, &mut jdet, &mut normal);
                let dot = mapped_pt
                    .iter()
                    .zip(&normal)
                    .map(|(i, j)| i * j)
                    .sum::<f64>();
                assert!(dot > 0.0);
            }
        }
    }

    #[test]
    fn test_normal_is_unit() {
        for i in 0..3 {
            let g = regular_sphere::<f64>(i);
            let points = vec![1.0 / 3.0, 1.0 / 3.0];
            let map = g.geometry_map(ReferenceCellType::Triangle, &points);
            let mut j = vec![0.0; 6];
            let mut jdet = vec![0.0];
            let mut mapped_pt = vec![0.0; 3];
            let mut normal = vec![0.0; 3];
            for i in 0..g.entity_count(ReferenceCellType::Triangle) {
                map.points(i, &mut mapped_pt);
                map.jacobians_dets_normals(i, &mut j, &mut jdet, &mut normal);
                let dot = normal.iter().zip(&normal).map(|(i, j)| i * j).sum::<f64>();
                assert_relative_eq!(dot, 1.0, epsilon = 1e-10);
            }
        }
    }
}
