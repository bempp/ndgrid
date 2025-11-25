//! Geometry where each entity of a given dimension is represented by the same element
#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::types::Scalar;
#[cfg(feature = "serde")]
use ndelement::{
    ciarlet::{CiarletElement, lagrange},
    map::IdentityMap,
    traits::FiniteElement,
    types::Continuity,
};
use ndelement::{
    reference_cell,
    traits::{ElementFamily, MappedFiniteElement},
    types::ReferenceCellType,
};
use rlst::{DynArray, rlst_dynamic_array};
use std::fmt::{Debug, Formatter};

/// Single element geometry
pub struct SingleElementGeometry<T: Scalar, E: MappedFiniteElement> {
    points: DynArray<T, 2>,
    cells: DynArray<usize, 2>,
    elements: Vec<E>,
}

impl<T: Scalar, E: MappedFiniteElement> Debug for SingleElementGeometry<T, E> {
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
pub struct SerializableGeometry<T: Scalar + serde::Serialize>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    points: (Vec<T>, [usize; 2]),
    cells: (Vec<usize>, [usize; 2]),
    elements: Vec<(ReferenceCellType, usize)>,
}

#[cfg(feature = "serde")]
impl<T: Scalar + serde::Serialize> ConvertToSerializable
    for SingleElementGeometry<T, CiarletElement<T, IdentityMap, T>>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    type SerializableType = SerializableGeometry<T>;
    fn to_serializable(&self) -> SerializableGeometry<T> {
        SerializableGeometry {
            points: (self.points.data().unwrap().to_vec(), self.points.shape()),
            cells: (self.cells.data().unwrap().to_vec(), self.cells.shape()),
            elements: self
                .elements
                .iter()
                .map(|e| (e.cell_type(), e.lagrange_superdegree()))
                .collect::<Vec<_>>(),
        }
    }
    fn from_serializable(s: SerializableGeometry<T>) -> Self {
        Self {
            points: {
                let (data, shape) = &s.points;
                let mut p = DynArray::<T, 2>::from_shape(*shape);
                p.data_mut().unwrap().copy_from_slice(data.as_slice());
                p
            },
            cells: {
                let (data, shape) = &s.cells;
                let mut c = DynArray::<usize, 2>::from_shape(*shape);
                c.data_mut().unwrap().copy_from_slice(data.as_slice());
                c
            },
            elements: s
                .elements
                .iter()
                .map(|(cell, degree)| lagrange::create(*cell, *degree, Continuity::Standard))
                .collect::<Vec<_>>(),
        }
    }
}

impl<T: Scalar, E: MappedFiniteElement> SingleElementGeometry<T, E> {
    /// Create single element geometry
    pub fn new(
        cell_type: ReferenceCellType,
        points: DynArray<T, 2>,
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
        let mut cells = rlst_dynamic_array!(
            usize,
            [points_per_cell, cells_input.len() / points_per_cell]
        );
        cells.data_mut().unwrap().copy_from_slice(cells_input);
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
    pub fn points(&self) -> &DynArray<T, 2> {
        &self.points
    }
    /// Cells
    pub fn cells(&self) -> &DynArray<usize, 2> {
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
