//! Geometry where each entity of a given dimension is represented by the same element
#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::types::{Array2D, RealScalar};
#[cfg(feature = "serde")]
use ndelement::{
    ciarlet::{lagrange, CiarletElement},
    types::Continuity,
};
use ndelement::{
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::ReferenceCellType,
};
#[cfg(feature = "serde")]
use rlst::RawAccess;
use rlst::{rlst_dynamic_array2, RawAccessMut, Shape};
use std::fmt::{Debug, Formatter};

/// Single element geometry
pub struct SingleElementGeometry<T: RealScalar, E: FiniteElement> {
    points: Array2D<T>,
    cells: Array2D<usize>,
    elements: Vec<E>,
}

impl<T: RealScalar, E: FiniteElement> Debug for SingleElementGeometry<T, E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        f.debug_struct("SingleElementGeometry")
            .field("points", &self.points)
            .field("cells", &self.cells)
            .finish()
    }
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
#[serde(bound = "for<'de2> T: serde::Deserialize<'de2>")]
pub struct SerializableGeometry<T: RealScalar + serde::Serialize>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    points: (Vec<T>, [usize; 2]),
    cells: (Vec<usize>, [usize; 2]),
    elements: Vec<(ReferenceCellType, usize)>,
}

#[cfg(feature = "serde")]
impl<T: RealScalar + serde::Serialize> ConvertToSerializable
    for SingleElementGeometry<T, CiarletElement<T>>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    type SerializableType = SerializableGeometry<T>;
    fn to_serializable(&self) -> SerializableGeometry<T> {
        SerializableGeometry {
            points: (self.points.data().to_vec(), self.points.shape()),
            cells: (self.cells.data().to_vec(), self.cells.shape()),
            elements: self
                .elements
                .iter()
                .map(|e| (e.cell_type(), e.embedded_superdegree()))
                .collect::<Vec<_>>(),
        }
    }
    fn from_serializable(s: SerializableGeometry<T>) -> Self {
        Self {
            points: {
                let (data, shape) = &s.points;
                let mut p = rlst_dynamic_array2!(T, *shape);
                p.data_mut().copy_from_slice(data.as_slice());
                p
            },
            cells: {
                let (data, shape) = &s.cells;
                let mut p = rlst_dynamic_array2!(usize, *shape);
                p.data_mut().copy_from_slice(data.as_slice());
                p
            },
            elements: s
                .elements
                .iter()
                .map(|(cell, degree)| lagrange::create(*cell, *degree, Continuity::Standard))
                .collect::<Vec<_>>(),
        }
    }
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
            .filter(|i| !i.is_empty())
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
