//! Geometry where each entity of a given dimension is represented by the same element
use crate::types::{Array2D, Array2DBorrowed, RealScalar};
use rlst::{rlst_dynamic_array2, RawAccessMut, Shape};
use std::fmt::{Debug, Formatter};
use ndelement::traits::FiniteElement;

/// Single element geometry
pub struct MixedGeometry<T: RealScalar, E: FiniteElement> {
    points: Array2D<T>,
    cells: Array2D<usize>,
    elements: Vec<E>,
}

impl<T: RealScalar, E: FiniteElement> Debug for MixedGeometry<T, E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        f.debug_struct("MixedGeometry")
            .field("points", &self.points)
            .field("cells", &self.cells)
            .finish()
    }
}

// TODO: Serializable

impl<T: RealScalar, E: FiniteElement> MixedGeometry<T, E> {
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
            .filter(|i| !i.is_empty())
        {
            for c in et.iter().skip(1) {
                if *c != et[0] {
                    panic!("Unsupported cell type for MixedGeometry: {cell_type:?}");
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
    /// Geometric dimension
    pub fn dim(&self) -> usize {
        self.points().shape()[0]
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

/// Single element geometry with borrowed data
pub struct MixedGeometryBorrowed<'a, T: RealScalar, E: FiniteElement> {
    points: Array2DBorrowed<'a, T>,
    cells: Array2DBorrowed<'a, usize>,
    elements: Vec<E>,
}

impl<T: RealScalar, E: FiniteElement> Debug for MixedGeometryBorrowed<'_, T, E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        f.debug_struct("MixedGeometryBorrowed")
            .field("points", &self.points)
            .field("cells", &self.cells)
            .finish()
    }
}

impl<'a, T: RealScalar, E: FiniteElement> MixedGeometryBorrowed<'a, T, E> {
    /// Create single element geometry
    pub fn new(
        points: Array2DBorrowed<'a, T>,
        cells: Array2DBorrowed<'a, usize>,
        elements: Vec<E>,
    ) -> Self {
        Self {
            points,
            cells,
            elements,
        }
    }

    /// Geometric dimension
    pub fn dim(&self) -> usize {
        self.points().shape()[0]
    }
    /// Points
    pub fn points(&self) -> &Array2DBorrowed<'_, T> {
        &self.points
    }
    /// Cells
    pub fn cells(&self) -> &Array2DBorrowed<'_, usize> {
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
