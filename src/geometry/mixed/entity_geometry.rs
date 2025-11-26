//! Geometry where each entity of a given dimension is represented by the same element
use super::MixedGeometry;
use crate::{
    geometry::{Point, PointIter},
    traits::Geometry,
    types::Scalar,
};
use ndelement::traits::MappedFiniteElement;

/// Geometry of an entity
#[derive(Debug)]
pub struct MixedEntityGeometry<'a, T: Scalar, E: MappedFiniteElement> {
    geometry: &'a MixedGeometry<T, E>,
    element_index: usize,
    cell_index: usize,
    sub_entity_dim: usize,
    sub_entity_index: usize,
}

impl<'a, T: Scalar, E: MappedFiniteElement> MixedEntityGeometry<'a, T, E> {
    /// Create new
    pub fn new(
        geometry: &'a MixedGeometry<T, E>,
        element_index: usize,
        cell_index: usize,
        sub_entity_dim: usize,
        sub_entity_index: usize,
    ) -> Self {
        Self {
            geometry,
            element_index,
            cell_index,
            sub_entity_dim,
            sub_entity_index,
        }
    }
}

impl<T: Scalar, E: MappedFiniteElement> Geometry for MixedEntityGeometry<'_, T, E> {
    type T = T;
    type Point<'a>
        = Point<'a, T>
    where
        Self: 'a;
    type PointIter<'a>
        = PointIter<'a, T>
    where
        Self: 'a;

    fn points(&self) -> PointIter<'_, T> {
        let gdim = self.geometry.dim();
        let mut pts = vec![];
        for index in self
            .geometry
            .element(self.element_index)
            .entity_closure_dofs(self.sub_entity_dim, self.sub_entity_index)
            .unwrap()
        {
            let i = self.geometry.cells(self.element_index)[[*index, self.cell_index]];
            pts.push((
                i,
                &self.geometry.points().data().unwrap()[i * gdim..(i + 1) * gdim],
            ))
        }

        PointIter::new(pts)
    }

    fn point_count(&self) -> usize {
        self.geometry.cells(self.element_index).shape()[0]
    }
    fn degree(&self) -> usize {
        self.geometry
            .element(self.element_index)
            .lagrange_superdegree()
    }
}
