//! Geometry where each entity of a given dimension is represented by the same element
#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::types::Scalar;
use itertools::izip;
#[cfg(feature = "serde")]
use ndelement::{
    ciarlet::{CiarletElement, lagrange},
    map::IdentityMap,
    types::Continuity,
};
use ndelement::{
    reference_cell,
    traits::{ElementFamily, MappedFiniteElement},
    types::ReferenceCellType,
};
use rlst::{DynArray, rlst_dynamic_array};
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};

/// Single element geometry
pub struct MixedGeometry<T: Scalar, E: MappedFiniteElement> {
    points: DynArray<T, 2>,
    cells: Vec<DynArray<usize, 2>>,
    elements: Vec<HashMap<ReferenceCellType, E>>,
    pub(crate) insertion_indices_to_element_indices: Vec<usize>,
    pub(crate) insertion_indices_to_cell_indices: Vec<usize>,
}

impl<T: Scalar, E: MappedFiniteElement> Debug for MixedGeometry<T, E> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        f.debug_struct("MixedGeometry")
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
    cells: Vec<(Vec<usize>, [usize; 2])>,
    elements: Vec<HashMap<ReferenceCellType, usize>>,
    insertion_indices_to_element_indices: Vec<usize>,
    insertion_indices_to_cell_indices: Vec<usize>,
}

#[cfg(feature = "serde")]
impl<T: Scalar + serde::Serialize> ConvertToSerializable
    for MixedGeometry<T, CiarletElement<T, IdentityMap, T>>
where
    for<'de2> T: serde::Deserialize<'de2>,
{
    type SerializableType = SerializableGeometry<T>;
    fn to_serializable(&self) -> SerializableGeometry<T> {
        SerializableGeometry {
            points: (self.points.data().unwrap().to_vec(), self.points.shape()),
            cells: self
                .cells
                .iter()
                .map(|c| (c.data().unwrap().to_vec(), c.shape()))
                .collect::<Vec<_>>(),
            elements: self
                .elements
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|(b, c)| (*b, c.lagrange_superdegree()))
                        .collect::<HashMap<_, _>>()
                })
                .collect::<Vec<_>>(),
            insertion_indices_to_element_indices: self.insertion_indices_to_element_indices.clone(),
            insertion_indices_to_cell_indices: self.insertion_indices_to_cell_indices.clone(),
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
            cells: s
                .cells
                .iter()
                .map(|(data, shape)| {
                    let mut c = DynArray::<usize, 2>::from_shape(*shape);
                    c.data_mut().unwrap().copy_from_slice(data.as_slice());
                    c
                })
                .collect::<Vec<_>>(),
            elements: s
                .elements
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|(cell, degree)| {
                            (
                                *cell,
                                lagrange::create(*cell, *degree, Continuity::Standard),
                            )
                        })
                        .collect::<HashMap<_, _>>()
                })
                .collect::<Vec<_>>(),
            insertion_indices_to_element_indices: s.insertion_indices_to_element_indices,
            insertion_indices_to_cell_indices: s.insertion_indices_to_cell_indices,
        }
    }
}

impl<T: Scalar, E: MappedFiniteElement> MixedGeometry<T, E> {
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
                c.data_mut().unwrap().copy_from_slice(&c_in);
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
