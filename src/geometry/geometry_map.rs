//! Geometry map
use crate::{
    traits::GeometryMap as GeometryMapTrait,
    types::{Array2D, ArrayND, RealScalar},
};
use ndelement::{reference_cell, traits::FiniteElement, types::ReferenceCellType};
use rlst::{
    rlst_dynamic_array2, rlst_dynamic_array4, DefaultIterator, RandomAccessByRef, RlstScalar, Shape,
};

/// Single element geometry
#[derive(Debug)]
pub struct GeometryMap<'a, T: RealScalar> {
    geometry_points: &'a Array2D<T>,
    entities: &'a Array2D<usize>,
    tdim: usize,
    gdim: usize,
    table: ArrayND<4, T>,
}

fn det<T: RlstScalar>(m: &Array2D<T>) -> T {
    match m.shape()[0] {
        0 => T::from(1.0).unwrap(),
        1 => m[[0, 0]],
        2 => m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]],
        3 => {
            m[[0, 0]] * (m[[1, 1]] * m[[2, 2]] - m[[1, 2]] * m[[2, 1]])
                - m[[0, 1]] * (m[[1, 0]] * m[[2, 2]] - m[[1, 2]] * m[[2, 0]])
                + m[[0, 2]] * (m[[1, 0]] * m[[2, 1]] - m[[1, 1]] * m[[2, 0]])
        }
        _ => {
            unimplemented!();
        }
    }
}

fn norm<T: RlstScalar>(vector: &[T]) -> T {
    vector.iter().map(|&i| i * i).sum::<T>().sqrt()
}

impl<'a, T: RealScalar> GeometryMap<'a, T> {
    /// Create new
    pub fn new<A2D: RandomAccessByRef<2, Item = T> + Shape<2>>(
        element: &impl FiniteElement<CellType = ReferenceCellType, T = T>,
        points: &A2D,
        geometry_points: &'a Array2D<T>,
        entities: &'a Array2D<usize>,
    ) -> Self {
        let tdim = reference_cell::dim(element.cell_type());
        assert_eq!(points.shape()[0], tdim);
        let gdim = geometry_points.shape()[0];
        let npoints = points.shape()[1];

        let mut table = rlst_dynamic_array4!(T, element.tabulate_array_shape(1, npoints));
        element.tabulate(points, 1, &mut table);

        Self {
            geometry_points,
            entities,
            tdim,
            gdim,
            table,
        }
    }
}

impl<'a, T: RealScalar> GeometryMapTrait for GeometryMap<'a, T> {
    type T = T;

    fn entity_topology_dimension(&self) -> usize {
        self.tdim
    }
    fn geometry_dimension(&self) -> usize {
        self.gdim
    }
    fn point_count(&self) -> usize {
        self.table.shape()[1]
    }
    fn points(&self, entity_index: usize, points: &mut [T]) {
        let npts = self.table.shape()[1];
        assert_eq!(points.len(), self.gdim * npts);

        points.fill(T::from(0.0).unwrap());
        // TODO: can I use rlst better here?
        for i in 0..self.entities.shape()[0] {
            let v = self.entities[[i, entity_index]];
            for point_index in 0..npts {
                let t = *self.table.get([0, point_index, i, 0]).unwrap();
                for gd in 0..self.gdim {
                    points[gd + self.gdim * point_index] +=
                        *self.geometry_points.get([gd, v]).unwrap() * t;
                }
            }
        }
    }
    fn jacobians(&self, entity_index: usize, jacobians: &mut [T]) {
        let npts = self.table.shape()[1];
        assert_eq!(jacobians.len(), self.gdim * self.tdim * npts);

        jacobians.fill(T::from(0.0).unwrap());
        // TODO: can I use rlst better here?
        for i in 0..self.entities.shape()[0] {
            let v = self.entities[[i, entity_index]];
            for point_index in 0..npts {
                for gd in 0..self.gdim {
                    for td in 0..self.tdim {
                        jacobians[gd + self.gdim * td + self.gdim * self.tdim * point_index] +=
                            *self.geometry_points.get([gd, v]).unwrap()
                                * *self.table.get([1 + td, point_index, i, 0]).unwrap();
                    }
                }
            }
        }
    }

    fn jacobians_dets_normals(
        &self,
        entity_index: usize,
        jacobians: &mut [Self::T],
        jdets: &mut [Self::T],
        normals: &mut [Self::T],
    ) {
        if self.tdim + 1 != self.gdim {
            panic!("Can only compute normal for entities where tdim + 1 == gdim");
        }
        let npts = self.table.shape()[1];
        assert_eq!(jacobians.len(), self.gdim * self.tdim * npts);
        assert_eq!(jdets.len(), npts);
        assert_eq!(normals.len(), self.gdim * npts);

        // TODO: can we remove this memory assignment?
        let mut temp = rlst_dynamic_array2!(T, [self.tdim, self.tdim]);

        for point_index in 0..npts {
            for gd in 0..self.gdim {
                for td in 0..self.tdim {
                    jacobians[gd + self.gdim * td + self.gdim * self.tdim * point_index] = self
                        .entities
                        .view()
                        .slice(1, entity_index)
                        .iter()
                        .enumerate()
                        .map(|(i, v)| {
                            *self.geometry_points.get([gd, v]).unwrap()
                                * *self.table.get([1 + td, point_index, i, 0]).unwrap()
                        })
                        .sum::<T>();
                }
            }
            for gd in 1..self.gdim {
                for td in 0..self.tdim {
                    temp[[gd - 1, td]] =
                        jacobians[gd + self.gdim * td + self.gdim * self.tdim * point_index];
                }
            }
            for gd in 0..self.gdim {
                normals[gd + self.gdim * point_index] =
                    if gd % 2 == 0 { det(&temp) } else { -det(&temp) };
                if gd < self.tdim {
                    for td in 0..self.tdim {
                        temp[[gd, td]] =
                            jacobians[gd + self.gdim * td + self.gdim * self.tdim * point_index];
                    }
                }
            }
            jdets[point_index] =
                norm(&normals[self.gdim * point_index..self.gdim * (point_index + 1)]);
            for gd in 0..self.gdim {
                normals[gd + self.gdim * point_index] /= jdets[point_index];
            }
        }
    }
}
