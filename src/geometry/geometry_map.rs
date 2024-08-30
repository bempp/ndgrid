//! Geometry map
use crate::{
    traits::GeometryMap as GeometryMapTrait,
    types::{Array2D, ArrayND, RealScalar},
};
use ndelement::{reference_cell, traits::FiniteElement, types::ReferenceCellType};
use rlst::UnsafeRandomAccessByRef;
use rlst::{rlst_dynamic_array4, RandomAccessByRef, RlstScalar, Shape};

/// Single element geometry
#[derive(Debug)]
pub struct GeometryMap<'a, T: RealScalar> {
    geometry_points: &'a Array2D<T>,
    entities: &'a Array2D<usize>,
    tdim: usize,
    gdim: usize,
    table: ArrayND<4, T>,
}

fn norm<T: RlstScalar>(vector: &[T]) -> T {
    vector.iter().map(|&i| i * i).sum::<T>().sqrt()
}

fn cross<T: RlstScalar>(mat: &[T], result: &mut [T]) {
    match mat.len() {
        0 => {}
        2 => {
            debug_assert!(result.len() == 2);
            unsafe {
                *result.get_unchecked_mut(0) = *mat.get_unchecked(1);
                *result.get_unchecked_mut(1) = -*mat.get_unchecked(0);
            }
        }
        6 => {
            debug_assert!(result.len() == 3);
            unsafe {
                *result.get_unchecked_mut(0) = *mat.get_unchecked(1) * *mat.get_unchecked(5)
                    - *mat.get_unchecked(2) * *mat.get_unchecked(4);
                *result.get_unchecked_mut(1) = *mat.get_unchecked(2) * *mat.get_unchecked(3)
                    - *mat.get_unchecked(0) * *mat.get_unchecked(5);
                *result.get_unchecked_mut(2) = *mat.get_unchecked(0) * *mat.get_unchecked(4)
                    - *mat.get_unchecked(1) * *mat.get_unchecked(3);
            }
        }
        _ => {
            unimplemented!();
        }
    }
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
        debug_assert!(points.shape()[0] == tdim);
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
        debug_assert!(points.len() == self.gdim * npts);

        points.fill(T::zero());
        for i in 0..self.entities.shape()[0] {
            let v = unsafe { *self.entities.get_unchecked([i, entity_index]) };
            for point_index in 0..npts {
                let t = unsafe { *self.table.get_unchecked([0, point_index, i, 0]) };
                for gd in 0..self.gdim {
                    unsafe {
                        *points.get_unchecked_mut(gd + self.gdim * point_index) +=
                            *self.geometry_points.get_unchecked([gd, v]) * t
                    };
                }
            }
        }
    }
    fn jacobians(&self, entity_index: usize, jacobians: &mut [T]) {
        let npts = self.table.shape()[1];
        debug_assert!(jacobians.len() == self.gdim * self.tdim * npts);

        jacobians.fill(T::zero());
        for i in 0..self.entities.shape()[0] {
            let v = unsafe { *self.entities.get_unchecked([i, entity_index]) };
            for point_index in 0..npts {
                for td in 0..self.tdim {
                    let t = unsafe { *self.table.get_unchecked([1 + td, point_index, i, 0]) };
                    for gd in 0..self.gdim {
                        unsafe {
                            *jacobians.get_unchecked_mut(
                                gd + self.gdim * td + self.gdim * self.tdim * point_index,
                            ) += *self.geometry_points.get_unchecked([gd, v]) * t
                        };
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
        debug_assert!(jacobians.len() == self.gdim * self.tdim * npts);
        debug_assert!(jdets.len() == npts);
        debug_assert!(normals.len() == self.gdim * npts);

        self.jacobians(entity_index, jacobians);

        for point_index in 0..npts {
            let j = unsafe {
                &jacobians.get_unchecked(
                    self.gdim * self.tdim * point_index..self.gdim * self.tdim * (point_index + 1),
                )
            };
            let n = unsafe {
                &mut normals
                    .get_unchecked_mut(self.gdim * point_index..self.gdim * (point_index + 1))
            };
            let jd = unsafe { jdets.get_unchecked_mut(point_index) };
            cross(j, n);
            *jd = norm(n);
            for n_i in n.iter_mut() {
                *n_i /= *jd;
            }
        }
    }
}
