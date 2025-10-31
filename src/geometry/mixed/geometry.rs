//! Geometry where each entity of a given dimension is represented by the same element
use crate::types::{Array2D, RealScalar};
use itertools::izip;
use ndelement::{
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::ReferenceCellType,
};
use rlst::{rlst_dynamic_array2, RawAccessMut, Shape};
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};

/// Single element geometry
pub struct MixedGeometry<T: RealScalar, E: FiniteElement> {
    points: Array2D<T>,
    cells: Vec<Array2D<usize>>,
    cell_types: Vec<ReferenceCellType>,
    elements: Vec<HashMap<ReferenceCellType, E>>,
    element_insertion_indices: Vec<Vec<usize>>,
    pub(crate) insertion_indices_to_element_indices: Vec<usize>,
    pub(crate) insertion_indices_to_cell_indices: Vec<usize>,
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
        cell_types_in: &[ReferenceCellType],
        points: Array2D<T>,
        cells_input: &[usize],
        element_families: &[impl ElementFamily<
            T = T,
            CellType = ReferenceCellType,
            FiniteElement = E,
        >],
        cell_families: &[usize],
    ) -> Self {
        let mut element_info = HashMap::new();
        let mut elements = vec![];
        let mut cells = vec![];
        let mut cell_types = vec![];
        let mut points_per_cell = vec![];
        let mut element_insertion_indices = vec![];
        let mut insertion_indices_to_element_indices = vec![];
        let mut insertion_indices_to_cell_indices = vec![];

        let mut start = 0;

        for (i, (findex, cell_type)) in izip!(cell_families, cell_types_in).enumerate() {
            let eindex = *element_info.entry((findex, cell_type)).or_insert_with(|| {
                let n = elements.len();

                let mut new_e = HashMap::new();
                for et in reference_cell::entity_types(*cell_type)
                    .iter()
                    .skip(1)
                    .filter(|i| !i.is_empty())
                {
                    for e in et {
                        new_e
                            .entry(*e)
                            .or_insert(element_families[*findex].element(*e));
                    }
                }
                points_per_cell.push(new_e[cell_type].dim());
                elements.push(new_e);
                cells.push(vec![]);
                element_insertion_indices.push(vec![]);
                n
            });

            insertion_indices_to_element_indices.push(eindex);
            insertion_indices_to_cell_indices.push(element_insertion_indices[eindex].len());
            element_insertion_indices[eindex].push(i);
            cell_types.push(*cell_type);
            for i in 0..points_per_cell[eindex] {
                cells[eindex].push(cells_input[start + i]);
            }
            start += points_per_cell[eindex];
        }
        let cells = izip!(points_per_cell, cells)
            .map(|(n, c_in)| {
                let mut c = rlst_dynamic_array2!(usize, [n, c_in.len() / n]);
                c.data_mut().copy_from_slice(&c_in);
                c
            })
            .collect::<Vec<_>>();
        Self {
            points,
            cells,
            cell_types,
            elements,
            element_insertion_indices,
            insertion_indices_to_element_indices,
            insertion_indices_to_cell_indices,
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
    pub fn cells(&self, element_index: usize) -> &Array2D<usize> {
        &self.cells[element_index]
    }
    /// Element for a cell
    pub fn element(&self, element_index: usize) -> &E {
        &self.elements[element_index][&self.cell_types[element_index]]
    }
    /// Number of elments
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
}
