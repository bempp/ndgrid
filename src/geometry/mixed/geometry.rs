//! Geometry where each entity of a given dimension is represented by the same element
use crate::types::RealScalar;
use itertools::izip;
use ndelement::{
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::ReferenceCellType,
};
use rlst::{DynArray, rlst_dynamic_array};
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};

/// Single element geometry
pub struct MixedGeometry<T: RealScalar, E: FiniteElement> {
    points: DynArray<T, 2>,
    cells: Vec<DynArray<usize, 2>>,
    elements: Vec<HashMap<ReferenceCellType, E>>,
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
        points: DynArray<T, 2>,
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
        let mut cell_counts = vec![];
        let mut points_per_cell = vec![];
        let mut insertion_indices_to_element_indices = vec![];
        let mut insertion_indices_to_cell_indices = vec![];

        let mut start = 0;

        for (findex, cell_type) in izip!(cell_families, cell_types_in) {
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
                cell_counts.push(0);
                n
            });

            insertion_indices_to_element_indices.push(eindex);
            insertion_indices_to_cell_indices.push(cell_counts[eindex]);
            for i in 0..points_per_cell[eindex] {
                cells[eindex].push(cells_input[start + i]);
            }
            cell_counts[eindex] += 1;
            start += points_per_cell[eindex];
        }
        let cells = izip!(cell_counts, points_per_cell, cells)
            .map(|(ncells, npts, c_in)| {
                let mut c = rlst_dynamic_array!(usize, [npts, ncells]);
                c.data_mut().copy_from_slice(&c_in);
                c
            })
            .collect::<Vec<_>>();
        Self {
            points,
            cells,
            elements,
            insertion_indices_to_element_indices,
            insertion_indices_to_cell_indices,
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
    pub fn cells(&self, element_index: usize) -> &DynArray<usize, 2> {
        &self.cells[element_index]
    }
    /// Element for a cell
    pub fn element(&self, element_index: usize) -> &E {
        for (ct, e) in &self.elements[element_index] {
            if reference_cell::dim(*ct) == self.dim() {
                return e;
            }
        }
        panic!("Could not find element");
    }
    /// Number of elments
    pub fn element_count(&self) -> usize {
        self.elements.len()
    }
}
