//! Grid builder

use super::SingleElementGrid;
use crate::{
    geometry::SingleElementGeometry,
    topology::single_type::SingleTypeTopology,
    traits::{Builder, GeometryBuilder, GridBuilder, TopologyBuilder},
    types::RealScalar,
};
use ndelement::{
    ciarlet::{lagrange, CiarletElement, LagrangeElementFamily},
    map::IdentityMap,
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::{Continuity, ReferenceCellType},
};
use rlst::{rlst_dynamic_array2, RawAccessMut};
use std::collections::{HashMap, HashSet};

/// Grid builder for a single element grid
#[derive(Debug)]
pub struct SingleElementGridBuilder<T: RealScalar> {
    gdim: usize,
    element_data: (ReferenceCellType, usize),
    element: CiarletElement<T, IdentityMap>,
    points_per_cell: usize,
    pub(crate) points: Vec<T>,
    cells: Vec<usize>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: Vec<usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
    point_indices: HashSet<usize>,
    cell_indices: HashSet<usize>,
}

impl<T: RealScalar> SingleElementGridBuilder<T> {
    /// Create a new grid builder
    pub fn new(gdim: usize, data: (ReferenceCellType, usize)) -> Self {
        Self::new_with_capacity(gdim, 0, 0, data)
    }

    /// Create a new grid builder with capacaty for a given number of points and cells
    pub fn new_with_capacity(
        gdim: usize,
        npoints: usize,
        ncells: usize,
        data: (ReferenceCellType, usize),
    ) -> Self {
        let points_per_cell = lagrange::create::<T>(data.0, data.1, Continuity::Standard).dim();
        let element = LagrangeElementFamily::<T>::new(data.1, Continuity::Standard).element(data.0);
        Self {
            gdim,
            element_data: data,
            element,
            points_per_cell,
            points: Vec::with_capacity(npoints * gdim),
            cells: Vec::with_capacity(ncells * points_per_cell),
            point_indices_to_ids: Vec::with_capacity(npoints),
            point_ids_to_indices: HashMap::new(),
            cell_indices_to_ids: Vec::with_capacity(ncells),
            cell_ids_to_indices: HashMap::new(),
            point_indices: HashSet::new(),
            cell_indices: HashSet::new(),
        }
    }
}

impl<T: RealScalar> Builder for SingleElementGridBuilder<T> {
    type Grid = SingleElementGrid<T, CiarletElement<T, IdentityMap>>;
    type T = T;
    type CellData<'a> = &'a [usize];
    type EntityDescriptor = ReferenceCellType;

    fn add_point(&mut self, id: usize, data: &[T]) {
        if data.len() != self.gdim {
            panic!("Point has wrong number of coordinates");
        }
        if self.point_indices.contains(&id) {
            panic!("Cannot add point with duplicate id.");
        }
        self.point_ids_to_indices
            .insert(id, self.point_indices_to_ids.len());
        self.point_indices.insert(id);
        self.point_indices_to_ids.push(id);
        self.points.extend_from_slice(data);
    }

    fn add_cell(&mut self, id: usize, cell_data: &[usize]) {
        if self.cell_indices.contains(&id) {
            panic!("Cannot add cell with duplicate id.");
        }
        assert_eq!(cell_data.len(), self.points_per_cell);
        self.cell_ids_to_indices
            .insert(id, self.cell_indices_to_ids.len());
        self.cell_indices.insert(id);
        self.cell_indices_to_ids.push(id);
        for id in cell_data {
            self.cells.push(self.point_ids_to_indices[id]);
        }
    }

    fn create_grid(&self) -> SingleElementGrid<T, CiarletElement<T, IdentityMap>> {
        let cell_vertices =
            self.extract_vertices(&self.cells, &[self.element_data.0], &[self.element_data.1]);

        // Add the vertex ids. But need to make sure that we don't have duplicates.
        // So first generate a set of vertex ids, then convert it to a vector.
        // We finally sort the vector to make sure that grid generation is reproducible across runs
        // as sets do not have a stable ordering.
        let mut vertex_ids = HashSet::<usize>::new();
        for v in &cell_vertices {
            vertex_ids.insert(*v);
        }

        let vertex_ids = {
            let mut tmp = Vec::<usize>::with_capacity(vertex_ids.len());
            tmp.extend(vertex_ids.iter());
            tmp.sort();
            tmp
        };

        let geometry = self.create_geometry(
            &self.point_indices_to_ids,
            &self.points,
            &self.cells,
            &[self.element_data.0],
            &[self.element_data.1],
        );

        let topology = self.create_topology(
            vertex_ids,
            (0..self.cell_count()).collect::<Vec<_>>(),
            &cell_vertices,
            &[self.element_data.0],
        );

        self.create_grid_from_topology_geometry(topology, geometry)
    }

    fn point_count(&self) -> usize {
        self.point_indices_to_ids.len()
    }
    fn cell_count(&self) -> usize {
        self.cell_indices_to_ids.len()
    }
    fn point_indices_to_ids(&self) -> &[usize] {
        &self.point_indices_to_ids
    }
    fn cell_indices_to_ids(&self) -> &[usize] {
        &self.cell_indices_to_ids
    }
    fn cell_points(&self, index: usize) -> &[usize] {
        &self.cells[self.points_per_cell * index..self.points_per_cell * (index + 1)]
    }
    fn cell_vertices(&self, index: usize) -> &[usize] {
        &self.cells[self.points_per_cell * index
            ..self.points_per_cell * index + reference_cell::entity_counts(self.element_data.0)[0]]
    }
    fn point(&self, index: usize) -> &[T] {
        &self.points[self.gdim * index..self.gdim * (index + 1)]
    }
    fn points(&self) -> &[T] {
        &self.points
    }
    fn cell_type(&self, _index: usize) -> ReferenceCellType {
        self.element_data.0
    }
    fn cell_degree(&self, _index: usize) -> usize {
        self.element_data.1
    }
    fn gdim(&self) -> usize {
        self.gdim
    }
    fn tdim(&self) -> usize {
        reference_cell::dim(self.element_data.0)
    }

    fn npts(&self, _cell_type: Self::EntityDescriptor, _degree: usize) -> usize {
        self.points_per_cell
    }
}

impl<T: RealScalar> GeometryBuilder for SingleElementGridBuilder<T> {
    type GridGeometry = SingleElementGeometry<T, CiarletElement<T, IdentityMap>>;
    fn create_geometry(
        &self,
        point_ids: &[usize],
        coordinates: &[Self::T],
        cell_points: &[usize],
        _cell_types: &[ReferenceCellType],
        _cell_degrees: &[usize],
    ) -> SingleElementGeometry<T, CiarletElement<T, IdentityMap>> {
        let npts = point_ids.len();
        let mut points = rlst_dynamic_array2!(T, [self.gdim(), npts]);
        points.data_mut().copy_from_slice(coordinates);

        let family = LagrangeElementFamily::<T>::new(self.element_data.1, Continuity::Standard);

        SingleElementGeometry::<T, CiarletElement<T, IdentityMap>>::new(
            self.element_data.0,
            points,
            cell_points,
            &family,
        )
    }
}

impl<T: RealScalar> TopologyBuilder for SingleElementGridBuilder<T> {
    type GridTopology = SingleTypeTopology;
    fn create_topology(
        &self,
        vertex_ids: Vec<usize>,
        cell_ids: Vec<usize>,
        cells: &[usize],
        _cell_types: &[ReferenceCellType],
    ) -> SingleTypeTopology {
        // Create a map from point ids to the corresponding positions in the points array.
        let vertex_ids_to_pos = {
            let mut tmp = HashMap::<usize, usize>::new();
            for (i, id) in vertex_ids.iter().enumerate() {
                tmp.insert(*id, i);
            }
            tmp
        };

        let cells = {
            let mut new_cells = Vec::<usize>::with_capacity(cells.len());
            for id in cells {
                new_cells.push(vertex_ids_to_pos[id]);
            }
            new_cells
        };

        // let cells = cells
        //     .iter()
        //     .map(|v| vertex_ids.iter().position(|i| *i == *v).unwrap())
        //     .collect::<Vec<_>>();

        SingleTypeTopology::new(
            &cells,
            self.element_data.0,
            Some(vertex_ids),
            Some(cell_ids),
        )
    }

    fn extract_vertices(
        &self,
        cell_points: &[usize],
        _cell_types: &[Self::EntityDescriptor],
        _cell_degrees: &[usize],
    ) -> Vec<usize> {
        let mut vertices = vec![];
        for v in 0..reference_cell::entity_counts(self.element_data.0)[0] {
            for d in self.element.entity_dofs(0, v).unwrap() {
                vertices.push(*d);
            }
        }
        let mut cell_vertices = vec![];
        let mut start = 0;
        let npoints = self.element.dim();
        while start < cell_points.len() {
            for v in &vertices {
                cell_vertices.push(cell_points[start + v])
            }
            start += npoints;
        }
        cell_vertices
    }
}

impl<T: RealScalar> GridBuilder for SingleElementGridBuilder<T> {
    fn create_grid_from_topology_geometry(
        &self,
        topology: <Self as TopologyBuilder>::GridTopology,
        geometry: <Self as GeometryBuilder>::GridGeometry,
    ) -> <Self as Builder>::Grid {
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

    #[test]
    fn test_non_contiguous_ids() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(4, &[1.0, 1.0, 0.0]);

        b.add_cell(0, &[0, 1, 2]);
        b.add_cell(2, &[1, 2, 4]);

        b.create_grid();
    }
}
