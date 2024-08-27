//! Geometry where each entity of a given dimension is represented by the same element
mod entity_geometry;
mod geometry;

pub use entity_geometry::SingleElementEntityGeometry;
#[cfg(feature = "serde")]
pub(crate) use geometry::SerializableGeometry;
pub use geometry::SingleElementGeometry;

#[cfg(test)]
mod test {
    use super::*;
    use crate::traits::{Geometry, Point};
    use approx::assert_relative_eq;
    use itertools::izip;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        types::Continuity,
    };
    use ndelement::{reference_cell, traits::FiniteElement, types::ReferenceCellType};
    use rlst::{rlst_dynamic_array2, Shape};
    use rlst::{DefaultIterator, RandomAccessMut};

    fn example_geometry_linear_interval1d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [1, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([0, 2]).unwrap() = 2.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Interval,
            points,
            &[0, 1, 1, 2],
            &family,
        )
    }

    fn example_geometry_linear_interval2d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [2, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 1.0;
        *points.get_mut([0, 2]).unwrap() = 2.0;
        *points.get_mut([1, 2]).unwrap() = 0.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Interval,
            points,
            &[0, 1, 1, 2],
            &family,
        )
    }

    fn example_geometry_linear_interval3d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [3, 3]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 1.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 2.0;
        *points.get_mut([1, 2]).unwrap() = 0.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Interval,
            points,
            &[0, 1, 1, 2],
            &family,
        )
    }

    fn example_geometry_quadratic_interval2d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [2, 5]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.5;
        *points.get_mut([1, 2]).unwrap() = 0.5;
        *points.get_mut([0, 3]).unwrap() = 2.0;
        *points.get_mut([1, 3]).unwrap() = 0.0;
        *points.get_mut([0, 4]).unwrap() = 1.5;
        *points.get_mut([1, 4]).unwrap() = 0.5;
        let family = LagrangeElementFamily::<f64>::new(2, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Interval,
            points,
            &[0, 1, 2, 1, 3, 4],
            &family,
        )
    }

    fn example_geometry_flat_triangle2d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [2, 4]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([0, 3]).unwrap() = 2.0;
        *points.get_mut([1, 3]).unwrap() = 1.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Triangle,
            points,
            &[0, 1, 2, 2, 1, 3],
            &family,
        )
    }

    fn example_geometry_flat_triangle3d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [3, 4]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([2, 0]).unwrap() = 1.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([2, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([2, 2]).unwrap() = 0.0;
        *points.get_mut([0, 3]).unwrap() = 2.0;
        *points.get_mut([1, 3]).unwrap() = 1.0;
        *points.get_mut([2, 3]).unwrap() = 0.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Triangle,
            points,
            &[0, 1, 2, 2, 1, 3],
            &family,
        )
    }

    fn example_geometry_quadratic_triangle2d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [2, 9]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([0, 3]).unwrap() = 0.5;
        *points.get_mut([1, 3]).unwrap() = 0.5;
        *points.get_mut([0, 4]).unwrap() = -0.2;
        *points.get_mut([1, 4]).unwrap() = 0.5;
        *points.get_mut([0, 5]).unwrap() = 0.5;
        *points.get_mut([1, 5]).unwrap() = 0.2;
        *points.get_mut([0, 6]).unwrap() = 1.0;
        *points.get_mut([1, 6]).unwrap() = 1.0;
        *points.get_mut([0, 7]).unwrap() = 1.2;
        *points.get_mut([1, 7]).unwrap() = 0.5;
        *points.get_mut([0, 8]).unwrap() = 0.5;
        *points.get_mut([1, 8]).unwrap() = 1.2;
        let family = LagrangeElementFamily::<f64>::new(2, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Triangle,
            points,
            &[0, 1, 2, 3, 4, 5, 2, 0, 6, 7, 8, 3],
            &family,
        )
    }

    fn example_geometry_quadrilateral2d() -> SingleElementGeometry<f64, CiarletElement<f64>> {
        let mut points = rlst_dynamic_array2!(f64, [2, 6]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 1.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 0.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([0, 3]).unwrap() = 1.0;
        *points.get_mut([1, 3]).unwrap() = 1.0;
        *points.get_mut([0, 4]).unwrap() = 2.0;
        *points.get_mut([1, 4]).unwrap() = 2.0;
        *points.get_mut([0, 5]).unwrap() = 3.0;
        *points.get_mut([1, 5]).unwrap() = 2.0;
        let family = LagrangeElementFamily::<f64>::new(1, Continuity::Standard);
        SingleElementGeometry::<f64, CiarletElement<f64>>::new(
            ReferenceCellType::Quadrilateral,
            points,
            &[0, 1, 2, 3, 1, 4, 3, 5],
            &family,
        )
    }

    macro_rules! make_tests {
        ($cellname:ident) => {
            paste::item! {
                #[test]
                fn [< test_points_ $cellname >]() {
                    //! Test that sub-entities have the correct points
                    let g = [< example_geometry_ $cellname >]();

                    let gdim = g.points().shape()[0];
                    let mut cell_vertices = vec![];
                    for cell in g.cells().col_iter() {
                        let mut pts = vec![];
                        for i in cell.iter() {
                            let mut pt = vec![];
                            for j in 0..gdim {
                                pt.push(g.points()[[j, i]]);
                            }
                            pts.push(pt)
                        }
                        cell_vertices.push(pts);
                    }
                    let conn = reference_cell::connectivity(g.element().cell_type());

                    let mut point = vec![0.0; gdim];
                    for (cell_i, vertices) in cell_vertices.iter().enumerate() {

                        for (dim, conn_dim) in conn.iter().enumerate() {
                            for (index, entity_vertices) in conn_dim.iter().enumerate() {
                                let entity = SingleElementEntityGeometry::<f64, CiarletElement<f64>>::new(
                                    &g, cell_i, dim, index,
                                );
                                for (v, pt) in izip!(&entity_vertices[0], entity.points()) {
                                    pt.coords(&mut point);
                                    for (i, j) in izip!(&point, &vertices[*v]) {
                                       assert_relative_eq!(*i, j);
                                    }
                                }


                            }
                        }
                    }
                }
            }
        };
    }

    make_tests!(linear_interval1d);
    make_tests!(linear_interval2d);
    make_tests!(linear_interval3d);
    make_tests!(quadratic_interval2d);
    make_tests!(flat_triangle2d);
    make_tests!(flat_triangle3d);
    make_tests!(quadratic_triangle2d);
    make_tests!(quadrilateral2d);
}
