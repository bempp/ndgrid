//! Geometry where each entity of a given dimension is represented by the same element
mod entity_geometry;
mod geometry;

pub use entity_geometry::MixedEntityGeometry;
pub use geometry::MixedGeometry;
#[cfg(feature = "serde")]
pub(crate) use geometry::SerializableGeometry;

#[cfg(test)]
mod test {
    use super::*;
    use crate::traits::{Geometry, Point};
    use approx::assert_relative_eq;
    use itertools::izip;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        map::IdentityMap,
        reference_cell,
        traits::FiniteElement,
        types::Continuity,
        types::ReferenceCellType,
    };
    use rlst::rlst_dynamic_array;

    fn example_geometry_triangles() -> MixedGeometry<f64, CiarletElement<f64, IdentityMap>> {
        let mut points = rlst_dynamic_array!(f64, [2, 7]);
        *points.get_mut([0, 0]).unwrap() = 0.0;
        *points.get_mut([1, 0]).unwrap() = 0.0;
        *points.get_mut([0, 1]).unwrap() = 2.0;
        *points.get_mut([1, 1]).unwrap() = 0.0;
        *points.get_mut([0, 2]).unwrap() = 1.0;
        *points.get_mut([1, 2]).unwrap() = 1.0;
        *points.get_mut([0, 3]).unwrap() = 2.5;
        *points.get_mut([1, 3]).unwrap() = 1.0;
        *points.get_mut([0, 4]).unwrap() = 0.0;
        *points.get_mut([1, 4]).unwrap() = 2.0;
        *points.get_mut([0, 5]).unwrap() = 1.0;
        *points.get_mut([1, 5]).unwrap() = 2.0;
        *points.get_mut([0, 6]).unwrap() = 2.0;
        *points.get_mut([1, 6]).unwrap() = 2.0;
        let families = vec![
            LagrangeElementFamily::<f64>::new(1, Continuity::Standard),
            LagrangeElementFamily::<f64>::new(2, Continuity::Standard),
        ];
        MixedGeometry::<f64, CiarletElement<f64, IdentityMap>>::new(
            &[ReferenceCellType::Triangle, ReferenceCellType::Triangle],
            points,
            &[0, 1, 4, 1, 4, 6, 5, 3, 2],
            &families,
            &[0, 1],
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
                    for n in 0..g.element_count() {
                        let mut cell_vertices = vec![];
                        for cell in g.cells(n).col_iter() {
                            let mut pts = vec![];
                            for i in cell.iter_ref() {
                                let mut pt = vec![];
                                for j in 0..gdim {
                                    pt.push(g.points()[[j, *i]]);
                                }
                                pts.push(pt)
                            }
                            cell_vertices.push(pts);
                        }
                        let conn = reference_cell::connectivity(g.element(n).cell_type());

                        let mut point = vec![0.0; gdim];
                        for (cell_i, vertices) in cell_vertices.iter().enumerate() {
                            for (dim, conn_dim) in conn.iter().enumerate() {
                                for (index, entity_vertices) in conn_dim.iter().enumerate() {
                                    let entity = MixedEntityGeometry::<f64, CiarletElement<f64, IdentityMap>>::new(
                                        &g, n, cell_i, dim, index,
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
            }
        };
    }

    make_tests!(triangles);
}
