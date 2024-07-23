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
use rlst::{rlst_dynamic_array2, RawAccess, RawAccessMut, Shape};

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
    /// Element for a sub-entity
    pub fn entity_element(&self, tdim: usize) -> &E {
        &self.elements[tdim - 1]
    }
    /// Element for a cell
    pub fn element(&self) -> &E {
        &self.elements[self.elements.len() - 1]
    }
}

#[derive(Clone, Copy)]
// TODO: move this
/// A points
pub struct Point<'a, T: RealScalar> {
    coordinates: &'a [T],
}

impl<'a, T: RealScalar> Point<'a, T> {
    /// Create new
    pub fn new(coordinates: &'a [T]) -> Self {
        Self { coordinates }
    }
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

/// Iterator over points
pub struct PointIter<'a, T: RealScalar> {
    points: Vec<&'a [T]>,
    index: usize,
}
impl<'a, T: RealScalar> PointIter<'a, T> {
    /// Create new
    pub fn new(points: Vec<&'a [T]>) -> Self {
        Self { points, index: 0 }
    }
}
impl<'a, T: RealScalar> Iterator for PointIter<'a, T> {
    type Item = Point<'a, T>;

    fn next(&mut self) -> Option<Point<'a, T>> {
        self.index += 1;
        if self.index <= self.points.len() {
            Some(Point::new(self.points[self.index - 1]))
        } else {
            None
        }
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

    type PointIter<'a> = PointIter<'a, T> where Self: 'a;

    fn points(&self) -> PointIter<'_, T> {
        let gdim = self.geometry.points().shape()[0];
        let mut pts = vec![];
        for index in self
            .geometry
            .element()
            .entity_closure_dofs(self.sub_entity_dimension, self.sub_entity_index)
            .unwrap()
        {
            let i = self.geometry.cells()[[*index, self.cell_index]];
            pts.push(&self.geometry.points().data()[i * gdim..(i + 1) * gdim])
        }

        PointIter::new(pts)
    }

    fn point_count(&self) -> usize {
        self.geometry.cells().shape()[0]
    }

    fn volume(&self) -> usize {
        unimplemented!();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;
    use itertools::izip;
    use ndelement::{
        ciarlet::{CiarletElement, LagrangeElementFamily},
        types::Continuity,
    };
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
}
