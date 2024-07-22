//! Geometry where each entity of a given dimension is represented by the same element
use crate::{
    traits::{Geometry, Point as PointTrait},
    types::{Array2D, RealScalar},
};
use ndelement::{
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::ReferenceCellType,
};
use rlst::{rlst_dynamic_array2, RawAccessMut};

use std::iter::Copied;

/// Single element geometry
pub struct SingleElementGeometry<T: RealScalar, E: FiniteElement> {
    points: Array2D<T>,
    cells: Array2D<usize>,
    elements: Vec<E>,
}

impl<T: RealScalar, E: FiniteElement> SingleElementGeometry<T, E> {
    /// Create single element geometry
    pub fn new(
        cell_type: ReferenceCellType,
        points: Array2D<T>,
        cells_input: &[usize],
        element_family: &impl ElementFamily<T = T, CellType = ReferenceCellType, FiniteElement = E>,
    ) -> Self {
        let mut elements = vec![];
        for et in reference_cell::entity_types(cell_type)
            .iter()
            .skip(1)
            .filter(|i| i.len() > 0)
        {
            for c in et.iter().skip(1) {
                if *c != et[0] {
                    panic!("Unsupported cell type for SingleElementGeometry: {cell_type:?}");
                }
            }
            println!("{et:?}");
            elements.push(element_family.element(et[0]));
        }
        let points_per_cell = elements[elements.len() - 1].dim();
        let mut cells = rlst_dynamic_array2!(
            usize,
            [points_per_cell, cells_input.len() / points_per_cell]
        );
        cells.data_mut().copy_from_slice(cells_input);
        Self {
            points,
            cells,
            elements,
        }
    }
    /// Points
    pub fn points(&self) -> &Array2D<T> {
        &self.points
    }
    /// Cells
    pub fn cells(&self) -> &Array2D<usize> {
        &self.cells
    }
    /// Element
    pub fn element(&self, tdim: usize) -> &E {
        &self.elements[tdim]
    }
}

#[derive(Clone, Copy)]
// TODO: move this
/// A points
pub struct Point<'a, T: RealScalar> {
    coordinates: &'a [T],
}

impl<'a, T: RealScalar> PointTrait for Point<'a, T> {
    type T = T;

    fn dim(&self) -> usize {
        self.coordinates.len()
    }

    fn coords(&self, data: &mut [T]) {
        data.copy_from_slice(self.coordinates);
    }
}

/// Geometry of a cell
pub struct SingleElementEntityGeometry<'a, T: RealScalar, E: FiniteElement> {
    geometry: &'a SingleElementGeometry<T, E>,
    cell_index: usize,
    sub_entity_dimension: usize,
    sub_entity_index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement> SingleElementEntityGeometry<'a, T, E> {
    /// Create new
    pub fn new(
        geometry: &'a SingleElementGeometry<T, E>,
        cell_index: usize,
        sub_entity_dimension: usize,
        sub_entity_index: usize,
    ) -> Self {
        Self {
            geometry,
            cell_index,
            sub_entity_dimension,
            sub_entity_index,
        }
    }
}

impl<'g, T: RealScalar, E: FiniteElement> Geometry for SingleElementEntityGeometry<'g, T, E> {
    type Point<'a> = Point<'a, T> where Self: 'a;

    type PointIter<'a> = Copied<std::slice::Iter<'a, Point<'a, T>>> where Self: 'a;

    fn points(&self) -> Self::PointIter<'_> {
        unimplemented!();
    }

    fn point_count(&self) -> usize {
        unimplemented!();
    }

    fn volume(&self) -> usize {
        unimplemented!();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        types::Continuity,
    };
    use rlst::RandomAccessMut;

    fn example_geometry_triangle() -> SingleElementGeometry<f64, CiarletElement<f64>> {
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

    #[test]
    fn test_points_triangle() {
        let g = example_geometry_triangle();
        let cell0 = SingleElementEntityGeometry::<f64, CiarletElement<f64>>::new(&g, 0, 2, 0);
        let cell1 = SingleElementEntityGeometry::<f64, CiarletElement<f64>>::new(&g, 0, 2, 0);
    }
}
