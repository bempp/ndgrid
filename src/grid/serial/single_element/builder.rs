//! Grid builder

use super::SingleElementGrid;
use crate::{
    geometry::SingleElementGeometry, topology::serial::SingleTypeTopology, traits::Builder,
    types::RealScalar,
};
use ndelement::{
    ciarlet::{lagrange, CiarletElement, LagrangeElementFamily},
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::{Continuity, ReferenceCellType},
};
use rlst::{rlst_dynamic_array2, RawAccessMut};
use std::collections::HashMap;

/// Grid builder for a single element grid
pub struct SingleElementGridBuilder<T: RealScalar> {
    gdim: usize,
    element_data: (ReferenceCellType, usize),
    points_per_cell: usize,
    pub(crate) points: Vec<T>,
    cells: Vec<usize>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: Vec<usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
}

impl<T: RealScalar> SingleElementGridBuilder<T> {
    /// Create a new grid builder
    pub fn new(gdim: usize, data: (ReferenceCellType, usize)) -> Self {
        let points_per_cell = lagrange::create::<T>(data.0, data.1, Continuity::Standard).dim();

        Self {
            gdim,
            element_data: data,
            points_per_cell,
            points: vec![],
            cells: vec![],
            point_indices_to_ids: vec![],
            point_ids_to_indices: HashMap::new(),
            cell_indices_to_ids: vec![],
            cell_ids_to_indices: HashMap::new(),
        }
    }

    /// Create a new grid builder with capacaty for a given number of points and cells
    pub fn new_with_capacity(
        gdim: usize,
        npoints: usize,
        ncells: usize,
        data: (ReferenceCellType, usize),
    ) -> Self {
        let points_per_cell = lagrange::create::<T>(data.0, data.1, Continuity::Standard).dim();
        Self {
            gdim,
            element_data: data,
            points_per_cell,
            points: Vec::with_capacity(npoints * gdim),
            cells: Vec::with_capacity(ncells * points_per_cell),
            point_indices_to_ids: Vec::with_capacity(npoints),
            point_ids_to_indices: HashMap::new(),
            cell_indices_to_ids: Vec::with_capacity(ncells),
            cell_ids_to_indices: HashMap::new(),
        }
    }
}

impl<T: RealScalar> Builder for SingleElementGridBuilder<T> {
    type Grid = SingleElementGrid<T, CiarletElement<T>>;
    type T = T;
    type CellData<'a> = &'a [usize];

    fn add_point(&mut self, id: usize, data: &[T]) {
        if data.len() != self.gdim {
            panic!("Point has wrong number of coordinates");
        }
        if self.point_indices_to_ids.contains(&id) {
            panic!("Cannot add point with duplicate id.");
        }
        self.point_ids_to_indices
            .insert(id, self.point_indices_to_ids.len());
        self.point_indices_to_ids.push(id);
        self.points.extend_from_slice(data);
    }

    fn add_cell(&mut self, id: usize, cell_data: &[usize]) {
        if self.cell_indices_to_ids.contains(&id) {
            panic!("Cannot add cell with duplicate id.");
        }
        assert_eq!(cell_data.len(), self.points_per_cell);
        self.cell_ids_to_indices
            .insert(id, self.cell_indices_to_ids.len());
        self.cell_indices_to_ids.push(id);
        for id in cell_data {
            self.cells.push(self.point_ids_to_indices[id]);
        }
    }

    fn create_grid(self) -> SingleElementGrid<T, CiarletElement<T>> {
        let family = LagrangeElementFamily::<T>::new(self.element_data.1, Continuity::Standard);
        let element = family.element(self.element_data.0);
        let mut vertices = vec![];
        for v in 0..reference_cell::entity_counts(self.element_data.0)[0] {
            for d in element.entity_dofs(0, v).unwrap() {
                vertices.push(*d);
            }
        }

        let mut tcells = vec![];
        let mut start = 0;
        let npoints = element.dim();
        while start < self.cells.len() {
            for v in &vertices {
                tcells.push(self.cells[start + v])
            }
            start += npoints;
        }

        let topology = SingleTypeTopology::new(&tcells, self.element_data.0);

        let npts = self.point_indices_to_ids.len();
        let mut points = rlst_dynamic_array2!(T, [3, npts]);
        points.data_mut().copy_from_slice(&self.points);

        let geometry = SingleElementGeometry::<T, CiarletElement<T>>::new(
            self.element_data.0,
            points,
            &self.cells,
            &family,
        );

        SingleElementGrid::new(topology, geometry)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[should_panic]
    fn test_duplicate_point_id() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(2, &[0.0, 0.0, 0.0]);
        b.add_point(0, &[1.0, 0.0, 0.0]);
        b.add_point(1, &[0.0, 1.0, 0.0]);
        b.add_point(2, &[1.0, 1.0, 0.0]);
    }

    #[test]
    #[should_panic]
    fn test_duplicate_cell_id() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);

        b.add_cell(0, &[0, 1, 2]);
        b.add_cell(0, &[1, 2, 3]);
    }
}
